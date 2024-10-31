"request handler."
import os.path as osp

import torch
from django.conf import settings
from django.shortcuts import render
from markdown import markdown

from .forms import InChIorSMILESinput
from .models import GnnepcsaftPara, db_update
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
            # pylint: enable=E1101
            if len(comp) == 0:
                pred, output, inchi = prediction(query)
                molimg = plotmol(inchi, images_dir)
                comp = db_update(pred, inchi, comp, plots=[plotden, plotvp, molimg])
            comp = comp[0]
            pred = torch.tensor([comp.m, comp.sigma, comp.e])
            plotden = comp.plot_den
            plotvp = comp.plot_vp
            molimg = comp.plot_mol
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
        "den_uri": plotden,
        "vp_uri": plotvp,
        "mol_uri": molimg,
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
