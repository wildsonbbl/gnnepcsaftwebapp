"Module for utils like plotting data."
import csv
import os
import os.path as osp
import re
import sqlite3 as lite
import textwrap
from urllib.parse import quote
from urllib.request import HTTPError, urlopen

import numpy as np
import torch
from gnnepcsaft.configs.default import get_config
from gnnepcsaft.data.graph import from_InChI, inchitosmiles, smilestoinchi
from gnnepcsaft.data.graphdataset import Ramirez, ThermoMLDataset
from gnnepcsaft.train.models import PnaconvsParams, PNApcsaftL, ReadoutMLPParams
from gnnepcsaft.train.utils import calc_deg, rhovp_data
from langchain_core.prompts import ChatPromptTemplate
from langchain_google_genai import ChatGoogleGenerativeAI
from rdkit import Chem
from rdkit.Chem import AllChem

file_dir = osp.dirname(__file__)
dataset_dir = osp.join(file_dir, "data")
deg = calc_deg("ramirez", file_dir)
device = torch.device("cpu")


def make_datasets(_dataset_dir):
    "make datasets."
    dt = Ramirez(_dataset_dir + "/ramirez2022")
    _ra_data = {}
    for gh in dt:
        _ra_data[gh.InChI] = gh.para
    dt = ThermoMLDataset(_dataset_dir + "/thermoml")
    _tml_data = {}
    for gh in dt:
        _tml_data[gh.InChI] = [prop.numpy() for prop in [gh.rho, gh.vp]]
    return _ra_data, _tml_data


ra_data, tml_data = make_datasets(dataset_dir)


# pylint: disable=R0914
def plotdata(para: np.ndarray, inchi: str) -> tuple[list, list, list]:
    "Organize data for plotting."
    plotden, plotvp = [], []
    if inchi in tml_data:
        rho, vp = tml_data[inchi]
        ra_rho, ra_vp = np.ndarray, np.ndarray
        pred_rho, pred_vp = rhovp_data(para, rho, vp)
        ra = False
        if inchi in ra_data:
            ra = True
            para = ra_data[inchi].squeeze().numpy()
            ra_rho, ra_vp = rhovp_data(para, rho, vp)
        # plot rho data
        if ~np.all(rho == np.zeros_like(rho)):
            idx_p = abs(rho[:, 1] - 101325) < 15_000
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
                if ra:
                    plotden = np.concatenate(
                        [
                            plotden,
                            ra_rho[idx_p][idx][:, np.newaxis],
                        ],
                        -1,
                    )
                else:
                    plotden = np.concatenate(
                        [
                            plotden,
                            np.full(
                                (plotden.shape[0], 1),
                                np.nan,
                            ),
                        ],
                        -1,
                    )
        # plot vp data
        if ~np.all(vp == np.zeros_like(vp)) and vp.shape[0] > 1:
            idx = np.argsort(vp[:, 0], 0)
            plotvp = np.stack(
                [
                    vp[idx, 0],
                    vp[idx, -1] / 1000,  # Pressure in kPa
                    pred_vp[idx] / 1000,
                ],
                -1,
            )
            if ra:
                plotvp = np.concatenate(
                    [
                        plotvp,
                        ra_vp[idx][:, np.newaxis] / 1000,
                    ],
                    -1,
                )
            else:
                plotvp = np.concatenate(
                    [
                        plotvp,
                        np.full(
                            (plotvp.shape[0], 1),
                            np.nan,
                        ),
                    ],
                    -1,
                )

    return plotden, plotvp


def plotmol(inchi: str) -> str:
    "Make Mol block for 3Dmol."

    mol = Chem.MolFromInchi(inchi)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xF00D
    AllChem.EmbedMolecule(mol, params)
    mol = Chem.RemoveHs(mol, implicitOnly=False)
    imgmol = Chem.MolToMolBlock(mol)
    return imgmol


def resume_mol(inchi: str):
    "Describe the molecule with google's gemini."

    llm = ChatGoogleGenerativeAI(
        model="gemini-pro", google_api_key=os.environ["API_KEY"]
    )

    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/description/json?inchi="
        + quote(inchi, safe="")
    )
    try:
        with urlopen(url) as ans:
            ans = ans.read().decode("utf8").rstrip()
    except (TypeError, HTTPError, ValueError):
        ans = "no data available."

    query = """
            You are a chemistry expert who is given this InChI {inchi} to analyse.
            Make sure to answer each one of the bellow questions. 
            To be able to do that you are gonna need to take into account all the 
            organic groups known in chemistry, 
            the difference between a Lewis acid and base, the concept of a 
            hydrogen bond, a hydrogen bond donor
            and hydrogen bond acceptor. Once you have all this 
            information gathered, you will be able to answer. You can use
            the passage bellow as reference.

            QUESTIONS: '
            First, describe the molecule with this InChI in detail.
            
            Then, answer the following questions about this molecule with 
            a detailed explanation of the answer:
            
               - Is it a Lewis acid or base or both? 
               - Can it do hydrogen bonds? 
               - Is it a hydrogen bond donor or acceptor?'
               
            PASSAGE: '{ans}'

            """
    query = textwrap.dedent(query)

    prompt = ChatPromptTemplate.from_template(query)

    # pylint: disable=E1131
    chain = prompt | llm
    # pylint: enable=E1131
    response = chain.invoke({"inchi": inchi, "ans": ans})

    return response.content


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
                  ROUND(m, 2), \
                  ROUND(sigma,2), \
                  ROUND(e, 2), \
                  inchi\
                FROM gnnmodel_gnnepcsaftpara WHERE inchi=?",
                (inchi,),
            )
            smiles = inchitosmiles(inchi, False, False)
            result = cur.fetchall()
            if len(result) >= 0:
                para, _, _ = prediction(inchi)
                plotden, plotvp = plotdata(para, inchi)
                para = para.tolist()
                cur.execute(
                    """
                    INSERT INTO gnnmodel_gnnepcsaftpara 
                    (m, sigma, e, inchi, smiles) 
                    VALUES (?,?,?,?,?)
                    """,
                    (para[0], para[1], para[2], inchi, smiles),
                )
                if len(plotden) > 0:
                    for row in plotden:
                        print(row)
                        cur.execute(
                            """
                            INSERT INTO gnnmodel_thermomldendata 
                            (inchi, T, den_tml, den_gnn, den_ra) 
                            VALUES (?,?,?,?,?)
                            """,
                            (inchi, row[0], row[1], row[2], row[3]),
                        )
                if len(plotvp) > 0:
                    for row in plotvp:
                        print(row)
                        cur.execute(
                            """
                            INSERT INTO gnnmodel_thermomlvpdata 
                            (inchi, T, vp_tml, vp_gnn, vp_ra) 
                            VALUES (?,?,?,?,?)
                            """,
                            (inchi, row[0], row[1], row[2], row[3]),
                        )
            else:
                data.append(result[0])
    with open("./static/mydata.csv", "w", encoding="UTF-8") as f:
        writer = csv.writer(f, delimiter="|")

        writer.writerow(["m", "sigma", "e", "inchi", "smiles"])
        writer.writerows(data)


def prediction(query: str) -> tuple[torch.Tensor, bool, str]:
    "Predict ePC-SAFT parameters."
    config = get_config()
    config.hidden_dim = 128
    model = PNApcsaftL(
        pna_params=PnaconvsParams(
            propagation_depth=2,
            pre_layers=1,
            post_layers=3,
            deg=deg,
        ),
        mlp_params=ReadoutMLPParams(num_mlp_layers=1, num_para=3),
        config=config,
    )
    model.to("cpu")

    checkpoint = torch.load(
        osp.join(dataset_dir, "model.ckpt"), map_location="cpu", weights_only=False
    )
    model.load_state_dict(checkpoint["state_dict"])

    model.eval()

    inchi = checking_inchi(query)

    try:
        graph = from_InChI(inchi).to(device)
        with torch.no_grad():
            pred = model.forward(graph)[0]
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
