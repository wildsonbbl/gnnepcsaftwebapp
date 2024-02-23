"Module for utils like plotting data."
import base64
import os
import os.path as osp
import textwrap
from io import BytesIO
from urllib.parse import quote
from urllib.request import HTTPError, urlopen

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from gnnepcsaft.data.graphdataset import Ramirez, ThermoMLDataset
from gnnepcsaft.train.utils import rhovp_data
from langchain_core.prompts import ChatPromptTemplate
from langchain_google_genai import ChatGoogleGenerativeAI
from rdkit import Chem
from rdkit.Chem import Draw

file_dir = osp.dirname(__file__)
dataset_dir = osp.join(file_dir, "static/data")


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


sns.set_theme(style="ticks")


def pltcustom(ra, scale="linear", ylabel=""):
    """
    Add legend and lables for `plotdata`.
    """
    plt.xlabel("T (K)")
    plt.ylabel(ylabel)
    plt.title("")
    legend = ["ThermoML Archive", "GNNePCSAFT"]
    if ra:
        legend += ["Ramírez-Vélez et al. (2022)"]
    plt.legend(legend)
    plt.grid(False)
    plt.yscale(scale)


def pltline(x, y):
    "Line plot."
    return plt.plot(x, y, linewidth=0.5)


def pltscatter(x, y):
    "Scatter plot."
    return plt.scatter(x, y, marker="x", s=10, c="black")


# pylint: disable=R0914
def plotdata(para: np.ndarray, inchi: str) -> tuple[str, str]:
    "plot den and vp data if available."
    plotden = ""
    plotvp = ""
    if inchi in tml_data:
        rho, vp = tml_data[inchi]
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
                x = rho[idx, 0]
                y = rho[idx, -1]
                pltscatter(x, y)
                y = pred_rho[idx]
                pltline(x, y)
                if ra:
                    y = ra_rho[idx_p][idx]
                    pltline(x, y)
                pltcustom(ra, ylabel="Density (mol / m³)")
                sns.despine(trim=True)
                imgbio = BytesIO()
                plt.savefig(
                    imgbio,
                    dpi=300,
                    format="png",
                    bbox_inches="tight",
                    transparent=True,
                )
                imgbio.seek(0)
                plotden = base64.b64encode(imgbio.read()).decode("ascii")

                plt.close()
        # plot vp data
        if ~np.all(vp == np.zeros_like(vp)) and vp.shape[0] > 1:
            idx = np.argsort(vp[:, 0], 0)
            x = vp[idx, 0]
            y = vp[idx, -1] / 1000
            pltscatter(x, y)
            y = pred_vp[idx] / 1000
            pltline(x, y)
            if ra:
                y = ra_vp[idx] / 1000
                pltline(x, y)
            pltcustom(ra, ylabel="Vapor pressure (kPa)")
            sns.despine(trim=True)
            imgbio = BytesIO()
            plt.savefig(
                imgbio,
                dpi=300,
                format="png",
                bbox_inches="tight",
                transparent=True,
            )
            imgbio.seek(0)
            plotvp = base64.b64encode(imgbio.read()).decode("ascii")

            plt.close()
    return plotden, plotvp


def plotmol(inchi: str):
    "Plot molecule."

    mol = Chem.MolFromInchi(inchi)

    options = Draw.DrawingOptions()
    options.bgColor = None
    options.atomLabelMinFontSize = 4
    options.useFraction = 1

    img = Draw.MolToMPL(mol, options=options)
    imgbio = BytesIO()
    plt.axis("off")
    img.savefig(
        imgbio,
        dpi=100,
        format="png",
        bbox_inches="tight",
        transparent=True,
    )
    plt.close()
    imgbio.seek(0)
    return base64.b64encode(imgbio.read()).decode("ascii")


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
