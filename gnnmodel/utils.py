"Module for utils like plotting data."
import csv
import json
import os.path as osp
import re

import numpy as np
import onnxruntime as ort
from django.conf import settings
from gnnepcsaft.data.ogb_utils import smiles2graph
from gnnepcsaft.data.rdkit_util import assoc_number, inchitosmiles, mw, smilestoinchi
from gnnepcsaft.epcsaft.utils import pure_den_feos, pure_vp_feos
from rdkit.Chem import AllChem as Chem

from .models import GnnepcsaftPara, ThermoMLDenData, ThermoMLVPData

# lazy import
# import polars as pl

ort.set_default_logger_severity(3)

file_dir = osp.dirname(__file__)
dataset_dir = osp.join(file_dir, "data")


def make_dataset():
    "Make dict dataset for inference."
    import polars as pl  # # pylint: disable = C0415

    data = pl.read_parquet(osp.join(dataset_dir, "thermoml/raw/pure.parquet"))
    inchis = data.unique("inchi1")["inchi1"].to_list()
    _tml_data = {}
    for inchi in inchis:
        try:
            molw = mw(inchi)
        except (TypeError, ValueError) as e:
            print(f"Error for InChI:\n {inchi}", e, sep="\n\n", end="\n\n")
            continue
        vp = (
            data.filter(pl.col("inchi1") == inchi, pl.col("tp") == 3)
            .select("TK", "PPa", "phase", "tp", "m")
            .to_numpy()
        )
        rho = (
            data.filter(pl.col("inchi1") == inchi, pl.col("tp") == 1)
            .select("TK", "PPa", "phase", "tp", "m")
            .to_numpy()
        )

        rho[:, -1] *= 1000 / molw  # convert to mol/ m³

        _tml_data[inchi] = (rho, vp)
    return _tml_data


# pylint: disable=R0914
def plotdata(para: np.ndarray, inchi: str) -> tuple[list, list]:
    "Organize data for plotting."
    tml_data = make_dataset()
    plotden, plotvp = None, None
    if inchi in tml_data:
        rho, vp = tml_data[inchi]
        pred_rho, pred_vp = rhovp_data(para, rho, vp)
        # plot rho data
        if rho.shape[0] > 2:
            idx_p = abs(rho[:, 1] - 101325) < 1_000
            rho: np.ndarray = rho[idx_p]
            pred_rho = pred_rho[idx_p]

            if rho.shape[0] > 2:
                idx = np.argsort(rho[:, 0], 0)
                plotden = {
                    "T": rho[idx, 0].tolist(),
                    "TML": rho[idx, -1].tolist(),
                    "GNN": pred_rho[idx].tolist(),
                }
        # plot vp data
        if vp.shape[0] > 2:
            idx = np.argsort(vp[:, 0], 0)
            plotvp = {
                "T": vp[idx, 0].tolist(),
                "TML": (vp[idx, -1] / 1000).tolist(),  # Pressure in kPa
                "GNN": (pred_vp[idx] / 1000).tolist(),
            }

    return plotden, plotvp


def plotmol(inchi: str) -> str:
    "Make Mol block for 3Dmol."

    mol = Chem.MolFromInchi(inchi)
    mol = Chem.AddHs(mol)
    params = Chem.ETKDGv3()
    params.randomSeed = 0xF00D
    result = Chem.EmbedMolecule(mol, params)
    if result == 0:
        Chem.MMFFOptimizeMolecule(
            mol, maxIters=1000, nonBondedThresh=100, ignoreInterfragInteractions=False
        )
    # mol = Chem.RemoveHs(mol, implicitOnly=False)
    imgmol = Chem.MolToV3KMolBlock(mol)
    return imgmol


def para_update_database(app, schema_editor):  # pylint: disable=W0613
    "fn to update database with epcsaft parameters"
    tml_data = make_dataset()

    data = []
    for inchi in tml_data:
        try:
            smiles = inchitosmiles(inchi, False, False)
            para, _, _ = prediction(smiles)
        except ValueError as e:
            print(e, inchi)
            continue
        if para[0] is None:
            continue
        new_comp = GnnepcsaftPara(
            inchi=inchi,
            smiles=smiles,
            m=para[0],
            sigma=para[1],
            e=para[2],
            k_ab=para[3],
            e_ab=para[4],
            mu=para[5],
            na=para[6],
            nb=para[7],
        )
        data.append(
            [
                para[0],
                para[1],
                para[2],
                para[3],
                para[4],
                para[5],
                para[6],
                para[7],
                inchi,
                smiles,
            ]
        )
        new_comp.save()
    with open(
        settings.BASE_DIR / "gnnmodel/static/mydata.csv", "w", encoding="UTF-8"
    ) as f:
        writer = csv.writer(f, delimiter="|")

        writer.writerow(
            ["m", "sigma", "e", "k_ab", "e_ab", "mu", "na", "nb", "inchi", "smiles"]
        )
        writer.writerows(data)
    print("Updated database with epcsaft parameters")


def thermo_update_database(app, schema_editor):  # pylint: disable=W0613
    "fn to update database with plotden, plotvp"
    tml_data = make_dataset()

    for inchi in tml_data:
        comp = GnnepcsaftPara.objects.filter(inchi=inchi)  # pylint: disable=E1101
        if len(comp) == 0:
            continue
        comp = comp[0]
        para = np.asarray(
            [
                comp.m,
                comp.sigma,
                comp.e,
                comp.k_ab,
                comp.e_ab,
                comp.mu,
                comp.na,
                comp.nb,
            ]
        )
        plotden, plotvp = plotdata(para, inchi)
        if plotden:
            new_comp = ThermoMLDenData(inchi=inchi, den=json.dumps(plotden))
            new_comp.save()
        if plotvp:
            new_comp = ThermoMLVPData(inchi=inchi, vp=json.dumps(plotvp))
            new_comp.save()
    print("Updated database with plotden, plotvp")


def prediction(smiles: str) -> tuple[np.ndarray, bool, str]:
    "Predict ePC-SAFT parameters."
    msigmae_onnx = ort.InferenceSession(settings.STATIC_ROOT / "msigmae_7.onnx")
    assoc_onnx = ort.InferenceSession(settings.STATIC_ROOT / "assoc_8.onnx")
    inchi = checking_inchi(smiles)
    try:
        graph = smiles2graph(smiles)
        na, nb = assoc_number(inchi)
        x, edge_index, edge_attr = (
            graph["node_feat"],
            graph["edge_index"],
            graph["edge_feat"],
        )

        assoc = (
            10
            ** (
                assoc_onnx.run(
                    None,
                    {
                        "x": x,
                        "edge_index": edge_index,
                        "edge_attr": edge_attr,
                    },
                )[0][0]
                * np.asarray([-1.0, 1.0])
            )
        ).round(decimals=4)
        if na == 0 and nb == 0:
            assoc *= 0
        msigmae = msigmae_onnx.run(
            None,
            {
                "x": x,
                "edge_index": edge_index,
                "edge_attr": edge_attr,
            },
        )[0].round(decimals=4)[0]
        munanb = np.asarray([0.0, na, nb])
        pred = np.hstack([msigmae, assoc, munanb])
        output = True
    except (ValueError, TypeError, AttributeError, IndexError) as e:
        print("\n\n", e)
        print("error for query:", smiles)
        pred = [None]
        output = False

    return pred, output, inchi


def checking_inchi(query: str) -> str:
    "Check if query is inchi and return an inchi."
    inchi_check = re.search("^InChI=", query)
    inchi = query
    if not inchi_check:
        try:
            inchi = smilestoinchi(query)
        except ValueError as e:
            print(e)
            print("error for query:", query)
    return inchi


def rhovp_data(parameters: np.ndarray, rho: np.ndarray, vp: np.ndarray):
    """Calculates density and vapor pressure with ePC-SAFT"""
    parameters = np.abs(parameters)
    den = []
    if rho.shape[0] > 0:
        for state in rho:
            den += [pure_den_feos(parameters, state)]
    den = np.asarray(den)

    vpl = []
    if vp.shape[0] > 0:
        for state in vp:
            try:
                vpl += [pure_vp_feos(parameters, state)]
            except (AssertionError, RuntimeError):
                vpl += [np.nan]
    vp = np.asarray(vpl)

    return den, vp
