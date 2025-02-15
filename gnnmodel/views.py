"request handler."
import os.path as osp

from django.conf import settings
from django.shortcuts import render

from .forms import InChIorSMILESinput
from .models import GnnepcsaftPara, ThermoMLDenData, ThermoMLVPData
from .utils import checking_inchi, plotmol, prediction

file_dir = osp.dirname(__file__)
images_dir = osp.join(settings.MEDIA_ROOT, "images")

available_params = [
    "Segment number",
    "Segment diameter (Å)",
    "Dispersion energy (K)",
    "Association volume",
    "Association energy (K)",
    "Dipole moment (D)*",
    "Number of association site A",
    "Number of association site B",
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
                pred = pred.tolist()
            else:
                comp = comp[0]
                pred = [
                    comp.m,
                    comp.sigma,
                    comp.e,
                    comp.k_ab,
                    comp.e_ab,
                    comp.mu,
                    comp.na,
                    comp.nb,
                ]
            alldata = ThermoMLVPData.objects.filter(inchi=inchi).all()
            if len(alldata) > 0:
                plotvp = alldata[0].vp

            alldata = ThermoMLDenData.objects.filter(inchi=inchi).all()
            if len(alldata) > 0:
                plotden = alldata[0].den

            molimg = plotmol(inchi)
            output = True

    else:
        form = InChIorSMILESinput()

    context = {
        "form": form,
        "predicted_para": (
            [
                (paraname, round(para, 4))
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
