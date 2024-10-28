"Module for utils like plotting data."
import base64
import datetime
import os
import os.path as osp
import random
import re
import string
import textwrap
from io import BytesIO
from urllib.parse import quote
from urllib.request import HTTPError, urlopen

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import torch
from django.conf import settings
from gnnepcsaft.configs.default import get_config
from gnnepcsaft.data.graph import from_InChI, smilestoinchi
from gnnepcsaft.data.graphdataset import Ramirez, ThermoMLDataset
from gnnepcsaft.train.models import PnaconvsParams, PNApcsaftL, ReadoutMLPParams
from gnnepcsaft.train.utils import calc_deg, rhovp_data
from langchain_core.prompts import ChatPromptTemplate
from langchain_google_genai import ChatGoogleGenerativeAI
from rdkit import Chem
from rdkit.Chem import Draw

from .models import GnnepcsaftPara

file_dir = osp.dirname(__file__)
dataset_dir = osp.join(settings.MEDIA_ROOT, "data")
workdir = settings.MEDIA_ROOT
images_dir = osp.join(settings.MEDIA_ROOT, "images")

deg = calc_deg("ramirez", workdir)
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
    "plot and save fig of den and vp data if available."
    plotden, plotvp = "", ""
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
                basename = "".join(
                    random.choice(string.ascii_lowercase) for i in range(16)
                )
                suffix = datetime.datetime.now().strftime("%y%m%d_%H%M%S")
                plotden = "_".join([basename, suffix, ".png"])
                img_path = osp.join(images_dir, plotden)
                plt.savefig(
                    img_path,
                    dpi=300,
                    format="png",
                    bbox_inches="tight",
                    transparent=True,
                )
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
            basename = "".join(random.choice(string.ascii_lowercase) for i in range(16))
            suffix = datetime.datetime.now().strftime("%y%m%d_%H%M%S")
            plotvp = "_".join([basename, suffix, ".png"])
            img_path = osp.join(images_dir, plotvp)
            plt.savefig(
                img_path,
                dpi=300,
                format="png",
                bbox_inches="tight",
                transparent=True,
            )
            plt.close()
    return plotden, plotvp


def plotmol(inchi: str) -> str:
    "Plot and save fig of molecule."

    mol = Chem.MolFromInchi(inchi)

    options = Draw.DrawingOptions()
    options.bgColor = None
    options.atomLabelMinFontSize = 4
    options.useFraction = 1

    img = Draw.MolToMPL(mol, options=options)
    basename = "".join(random.choice(string.ascii_lowercase) for i in range(16))
    suffix = datetime.datetime.now().strftime("%y%m%d_%H%M%S")
    imgmol = "_".join([basename, suffix, ".png"])
    img_path = osp.join(images_dir, imgmol)
    plt.axis("off")
    img.savefig(
        img_path,
        dpi=100,
        format="png",
        bbox_inches="tight",
        transparent=True,
    )
    plt.close()
    return imgmol


def plotmol_temp(inchi: str) -> str:
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


def update_database():
    "fn to update database with epcsaft parameters, plotden, plotvp and plotmol"

    for inchi in tml_data:
        # pylint: disable=E1101
        comp = GnnepcsaftPara.objects.filter(inchi=inchi).all()
        # pylint: enable=E1101
        if len(comp) == 0:
            para, _, _ = prediction(inchi)
            plotden, plotvp = plotdata(para, inchi)
            pltmol = plotmol(inchi)
            db_update(para, inchi, comp, [plotden, plotvp, pltmol])


def db_update(pred, inchi, comp, plots):
    "Updates the gnnepcsaft db."
    if len(comp) == 0:
        new_comp = GnnepcsaftPara(
            inchi=inchi,
            m=pred[0],
            sigma=pred[1],
            e=pred[2],
            plot_den=plots[0],
            plot_vp=plots[1],
            plot_mol=plots[2],
        )
        new_comp.save()
        return [new_comp]


def prediction(query: str) -> tuple[torch.Tensor, bool, str]:
    "Predict ePC-SAFT parameters."

    model = PNApcsaftL(
        pna_params=PnaconvsParams(
            propagation_depth=2,
            pre_layers=1,
            post_layers=3,
            deg=deg,
        ),
        mlp_params=ReadoutMLPParams(num_mlp_layers=1, num_para=3),
        config=get_config(),
    )
    model.to("cpu")

    checkpoint = torch.load(osp.join(workdir, "model.ckpt"), map_location="cpu")
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
