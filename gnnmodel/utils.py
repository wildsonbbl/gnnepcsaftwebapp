"Module for utils like plotting data."

import csv
import json
import os.path as osp

import numpy as np
import onnxruntime as ort
from django.conf import settings
from gnnepcsaft.data.ogb_utils import smiles2graph
from gnnepcsaft.data.rdkit_util import assoc_number, inchitosmiles, mw, smilestoinchi
from gnnepcsaft.epcsaft.epcsaft_feos import (
    pure_den_feos,
    pure_h_lv_feos,
    pure_s_lv_feos,
    pure_vp_feos,
)
from rdkit.Chem import AllChem as Chem

from .models import GnnepcsaftPara, ThermoMLDenData, ThermoMLVPData

# lazy import
# import polars as pl

ort.set_default_logger_severity(3)

file_dir = osp.dirname(__file__)
dataset_dir = osp.join(file_dir, "data")

msigmae_onnx = ort.InferenceSession(settings.STATIC_ROOT / "msigmae_7.onnx")
assoc_onnx = ort.InferenceSession(settings.STATIC_ROOT / "assoc_8.onnx")


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
def plotdata(para: np.ndarray, inchi: str, tml_data: dict) -> tuple[list, list]:
    "Organize data for plotting."
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
    from tqdm import tqdm  # pylint: disable=C0415

    tml_data = make_dataset()

    data = []
    for inchi in tqdm(tml_data):
        try:
            smiles = inchitosmiles(inchi, False, False)
            para = prediction(smiles)
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
        writer = csv.writer(f, delimiter="|", lineterminator="\n")

        writer.writerow(
            ["m", "sigma", "e", "k_ab", "e_ab", "mu", "na", "nb", "inchi", "smiles"]
        )
        writer.writerows(data)
    print("Updated database with epcsaft parameters")


def thermo_update_database(app, schema_editor):  # pylint: disable=W0613
    "fn to update database with plotden, plotvp"
    from tqdm import tqdm  # pylint: disable=C0415

    tml_data = make_dataset()

    for inchi in tqdm(tml_data):
        molecule = GnnepcsaftPara.objects.filter(inchi=inchi)  # pylint: disable=E1101
        if len(molecule) == 0:
            continue
        molecule = molecule[0]
        para = np.asarray(
            [
                molecule.m,
                molecule.sigma,
                molecule.e,
                molecule.k_ab,
                molecule.e_ab,
                molecule.mu,
                molecule.na,
                molecule.nb,
            ]
        )
        plotden, plotvp = plotdata(para, inchi, tml_data)
        if plotden:
            new_comp = ThermoMLDenData(inchi=inchi, den=json.dumps(plotden))
            new_comp.save()
        if plotvp:
            new_comp = ThermoMLVPData(inchi=inchi, vp=json.dumps(plotvp))
            new_comp.save()
    print("Updated database with plotden, plotvp")


def prediction(smiles: str) -> np.ndarray:
    "Predict ePC-SAFT parameters."
    lower_bounds = np.asarray([1.0, 1.9, 50.0, 0.0, 0.0, 0, 0, 0])
    upper_bounds = np.asarray([25.0, 4.5, 550.0, 0.9, 5000.0, np.inf, np.inf, np.inf])

    inchi = smilestoinchi(smiles)

    graph = smiles2graph(smiles)
    na, nb = assoc_number(inchi)
    x, edge_index, edge_attr = (
        graph["node_feat"],
        graph["edge_index"],
        graph["edge_feat"],
    )

    assoc = 10 ** (
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
    if na == 0 and nb == 0:
        assoc *= 0
    msigmae = msigmae_onnx.run(
        None,
        {
            "x": x,
            "edge_index": edge_index,
            "edge_attr": edge_attr,
        },
    )[0][0]
    munanb = np.asarray([0.0, na, nb])
    pred = np.hstack([msigmae, assoc, munanb])
    pred = np.clip(pred, lower_bounds, upper_bounds)

    return pred


def rhovp_data(parameters: np.ndarray, rho: np.ndarray, vp: np.ndarray):
    """Calculates density and vapor pressure with ePC-SAFT"""

    all_pred_den = []
    if rho.shape[0] > 0:
        for state in rho:
            try:
                all_pred_den += [pure_den_feos(parameters, state)]
            except (AssertionError, RuntimeError):
                all_pred_den += [np.nan]
    all_pred_den = np.asarray(all_pred_den)

    all_pred_vp = []
    if vp.shape[0] > 0:
        for state in vp:
            try:
                all_pred_vp += [pure_vp_feos(parameters, state)]
            except (AssertionError, RuntimeError):
                all_pred_vp += [np.nan]
    all_pred_vp = np.asarray(all_pred_vp)

    return all_pred_den, all_pred_vp


def custom_plot(
    parameters: list,
    temp_min: float,
    temp_max: float,
    pressure: float,
    checkboxes: list,
):
    """
    Custom plot function for ePC-SAFT parameters."
    args:
    parameters: list
    list with ePC-SAFT parameters
    temp_min: float
    minimum temperature in Kelvin
    temp_max: float
    maximum temperature in Kelvin
    pressure: float
    pressure in Pa
    checkboxes: list
    list with checks to plot in the order [density, vapor pressure, enthalpy, entropy]
    """
    temp_range = np.linspace(temp_min, temp_max, 100, dtype=np.float64)
    p_range = np.asarray([pressure] * 100, dtype=np.float64)
    states = np.stack((temp_range, p_range), 1)
    prop_fns = [
        pure_den_feos,
        pure_vp_feos,
        pure_h_lv_feos,
        pure_s_lv_feos,
    ]  # more prop later
    prop_ids = ["den_plot", "vp_plot", "h_lv_plot", "s_lv_plot"]
    prop_names = [
        "Liquid Density (mol / m³)",
        "Vapor pressure (Pa)",
        "Enthalpy of vaporization (kJ/mol)",
        "Entropy of vaporization (J/mol/K)",
    ]
    xlegendpos = [0, 0, 0, 0]
    all_plots = []

    for prop_fn, prop_id, prop_name, xpos, checkbox in zip(
        prop_fns,
        prop_ids,
        prop_names,
        xlegendpos,
        checkboxes,
    ):
        if checkbox:
            plot_data = {"T": [], "GNN": [], "TML": []}
            for state in states:
                try:

                    prop_for_state = prop_fn(
                        np.asarray(parameters, dtype=np.float64), state
                    )
                    plot_data["T"].append(state[0])
                    plot_data["GNN"].append(prop_for_state)
                except (AssertionError, RuntimeError) as e:
                    print(e)
            all_plots.append((json.dumps(plot_data), xpos, prop_name, prop_id))
    return all_plots
