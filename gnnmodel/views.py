"request handler."
import json
import os.path as osp

import torch
from django.conf import settings
from django.shortcuts import render
from markdown import markdown

from .forms import InChIorSMILESinput
from .models import GnnepcsaftPara, ThermoMLDenData, ThermoMLVPData, db_update
from .utils import checking_inchi, plotmol, prediction

file_dir = osp.dirname(__file__)
images_dir = osp.join(settings.MEDIA_ROOT, "images")

available_params = [
    "Segment number",
    "Segment diameter (Å)",
    "Dispersion energy (K)",
]


# Create your views here.
def estimator(request):
    "handle request"

    pred = None
    query = ""
    output = False
    plotden, plotvp, molimg = "", "", ""
    if request.method == "POST":
        form = InChIorSMILESinput(request.POST)

        if form.is_valid():
            query = form.cleaned_data["query"]
            inchi = checking_inchi(query)
            # pylint: disable=E1101
            comp = GnnepcsaftPara.objects.filter(inchi=inchi).all()
            if len(comp) == 0:
                pred, output, inchi = prediction(query)
                comp = db_update(pred, inchi)
            comp = comp[0]
            pred = torch.tensor([comp.m, comp.sigma, comp.e])
            alldata = ThermoMLVPData.objects.filter(inchi=inchi).all()
            if len(alldata) > 0:
                plotvp = {"T": [], "TML": [], "GNN": [], "RA": []}
                for row in alldata:
                    plotvp["T"].append(row.T)
                    plotvp["TML"].append(row.vp_tml)
                    plotvp["GNN"].append(row.vp_gnn)
                    plotvp["RA"].append(row.vp_ra)
                plotvp = json.dumps(plotvp)

            alldata = ThermoMLDenData.objects.filter(inchi=inchi).all()
            if len(alldata) > 0:
                plotden = {"T": [], "TML": [], "GNN": [], "RA": []}
                for row in alldata:
                    plotden["T"].append(row.T)
                    plotden["TML"].append(row.den_tml)
                    plotden["GNN"].append(row.den_gnn)
                    plotden["RA"].append(row.den_ra)
                plotden = json.dumps(plotden)

            molimg = plotmol(inchi)
            output = True

            with open(
                osp.join(file_dir, "templates/description.txt"), "w", encoding="utf-8"
            ) as file:
                file.write(inchi)
    else:
        form = InChIorSMILESinput()

    context = {
        "form": form,
        "predicted_para": (
            [
                (paraname, round(para.item(), 2))
                for para, paraname in zip(pred, available_params)
            ]
            if output
            else [(None, None)]
        ),
        "query": query,
        "output": output,
        "plotden": plotden != "",
        "plotvp": plotvp != "",
        "den_data": plotden,
        "vp_data": plotvp,
        "mol_data": molimg,
    }

    return render(request, "pred.html", context)


def homepage(request):
    "handle request"
    return render(request, "homepage.html")


def authorpage(request):
    "handle request"
    return render(request, "author.html")


def description(request):
    "handle request"

    html_output = ""

    with open(
        osp.join(file_dir, "templates/description.txt"), "r", encoding="utf-8"
    ) as file:
        inchi = file.readline()

    if inchi != "":
        html_output = markdown("Temporary disabled")

        with open(
            osp.join(file_dir, "templates/description.txt"), "w", encoding="utf-8"
        ) as file:
            file.write("")

    return render(request, "description.html", {"output": html_output})
