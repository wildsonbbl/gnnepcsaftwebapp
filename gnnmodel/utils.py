import base64
import os.path as osp
from io import BytesIO

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import torch

from model.data.graphdataset import Ramirez, ThermoMLDataset
from model.demo.utils import pltline, pltscatter
from model.train.utils import rhovp_data

file_dir = osp.dirname(__file__)
dt = Ramirez("./model/data/ramirez2022")
ra_data = {}
for gh in dt:
    ra_data[gh.InChI] = gh.para
dt = ThermoMLDataset("./model/data/thermoml")
tml_data = {}
for gh in dt:
    tml_data[gh.InChI] = [prop.numpy() for prop in [gh.rho, gh.vp]]


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

        if ~np.any(rho == np.zeros_like(rho)):
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
            if ~np.any(vp == np.zeros_like(vp)) and vp.shape[0] > 1:
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
