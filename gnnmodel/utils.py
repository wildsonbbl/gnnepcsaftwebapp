"Module for utils like plotting data."
import csv
import os.path as osp
import re
import sqlite3 as lite

import numpy as np
import onnxruntime as ort
import polars as pl
from gnnepcsaft.data.ogb_utils import smiles2graph
from gnnepcsaft.data.rdkit_util import assoc_number, inchitosmiles, mw, smilestoinchi
from gnnepcsaft.epcsaft.utils import pure_den_feos, pure_vp_feos
from rdkit.Chem import AllChem as Chem

ort.set_default_logger_severity(3)

file_dir = osp.dirname(__file__)
dataset_dir = osp.join(file_dir, "data")


def make_dataset():
    "Make dict dataset for inference."
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


tml_data = make_dataset()


# pylint: disable=R0914
def plotdata(para: np.ndarray, inchi: str) -> tuple[list, list]:
    "Organize data for plotting."
    plotden, plotvp = [], []
    if inchi in tml_data:
        rho, vp = tml_data[inchi]
        pred_rho, pred_vp = rhovp_data(para, rho, vp)
        # plot rho data
        if rho.shape[0] > 2:
            idx_p = abs(rho[:, 1] - 101325) < 1_000
            rho = rho[idx_p]
            pred_rho = pred_rho[idx_p]

            if rho.shape[0] > 1:
                idx = np.argsort(rho[:, 0], 0)
                plotden = np.stack(
                    [
                        rho[idx, 0],
                        rho[idx, -1],
                        pred_rho[idx],
                    ],
                    -1,
                )
        # plot vp data
        if vp.shape[0] > 2:
            idx = np.argsort(vp[:, 0], 0)
            plotvp = np.stack(
                [
                    vp[idx, 0],
                    vp[idx, -1] / 1000,  # Pressure in kPa
                    pred_vp[idx] / 1000,
                ],
                -1,
            )

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


def update_database():
    "fn to update database with epcsaft parameters, plotden, plotvp and plotmol"
    workdir = "/workspaces/webapp/gnnepcsaftwebapp"  # set mannualy
    con = lite.connect(osp.join(workdir, "mydatabase"))

    data = []
    for inchi in tml_data:
        with con:
            cur = con.cursor()
            cur.execute(
                "SELECT \
                  ROUND(m, 4), \
                  ROUND(sigma, 4), \
                  ROUND(e, 4), \
                  ROUND(k_ab, 4), \
                  ROUND(e_ab, 4), \
                  ROUND(mu, 4), \
                  ROUND(na, 4), \
                  ROUND(nb, 4), \
                  inchi, \
                  smiles \
                FROM gnnmodel_gnnepcsaftpara WHERE inchi=?",
                (inchi,),
            )
            smiles = inchitosmiles(inchi, False, False)
            result = cur.fetchall()
            if len(result) == 0:
                para, _, _ = prediction(inchi)
                plotden, plotvp = plotdata(para, inchi)
                para = para.tolist()
                cur.execute(
                    """
                    INSERT INTO gnnmodel_gnnepcsaftpara 
                    (m, sigma, e, k_ab, e_ab, mu, na, nb, inchi, smiles) 
                    VALUES (?,?,?,?,?,?,?,?,?,?)
                    """,
                    (
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
                    ),
                )
                if len(plotden) > 0:
                    for row in plotden:
                        cur.execute(
                            """
                            INSERT INTO gnnmodel_thermomldendata 
                            (inchi, T, den_tml, den_gnn) 
                            VALUES (?,?,?,?)
                            """,
                            (inchi, row[0], row[1], row[2]),
                        )
                if len(plotvp) > 0:
                    for row in plotvp:
                        cur.execute(
                            """
                            INSERT INTO gnnmodel_thermomlvpdata 
                            (inchi, T, vp_tml, vp_gnn) 
                            VALUES (?,?,?,?)
                            """,
                            (inchi, row[0], row[1], row[2]),
                        )
            else:
                data.append(result[0])
    with open("./static/mydata.csv", "w", encoding="UTF-8") as f:
        writer = csv.writer(f, delimiter="|")

        writer.writerow(
            ["m", "sigma", "e", "k_ab", "e_ab", "mu", "na", "nb", "inchi", "smiles"]
        )
        writer.writerows(data)


def prediction(query: str) -> tuple[np.ndarray, bool, str]:
    "Predict ePC-SAFT parameters."
    msigmae_onnx = ort.InferenceSession(osp.join(dataset_dir, "msigmae_6.onnx"))
    assoc_onnx = ort.InferenceSession(osp.join(dataset_dir, "assoc.onnx"))
    inchi = checking_inchi(query)
    try:
        graph = smiles2graph(query)
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
        print(e)
        print("error for query:", query)
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
        except (ValueError, TypeError, AttributeError, IndexError) as e:
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
