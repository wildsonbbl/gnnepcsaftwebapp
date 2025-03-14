"request handler."

import os.path as osp

from django.conf import settings
from django.shortcuts import render

from .forms import (
    CustomPlotCheckForm,
    CustomPlotConfigForm,
    InChIorSMILESinput,
    RhoCheckForm,
    VPCheckForm,
)
from .models import GnnepcsaftPara, ThermoMLDenData, ThermoMLVPData
from .utils import checking_inchi, custom_plot, plotmol, prediction

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
    custom_plots = []
    if request.method == "POST":
        form = InChIorSMILESinput(request.POST)
        plot_config = CustomPlotConfigForm(request.POST)
        plot_checkbox = CustomPlotCheckForm(request.POST)
        rho_checkbox = RhoCheckForm(request.POST)
        vp_checkbox = VPCheckForm(request.POST)

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
            plot_config.full_clean()
            plot_checkbox.full_clean()
            rho_checkbox.full_clean()
            vp_checkbox.full_clean()
            if plot_checkbox.cleaned_data["custom_plot_checkbox"]:
                custom_plots = custom_plot(
                    pred,
                    plot_config.cleaned_data["temp_min"],
                    plot_config.cleaned_data["temp_max"],
                    plot_config.cleaned_data["pressure"],
                    [
                        rho_checkbox.cleaned_data["rho_checkbox"],
                        vp_checkbox.cleaned_data["vp_checkbox"],
                    ],
                )

    else:
        form = InChIorSMILESinput()
        plot_config = CustomPlotConfigForm()
        plot_checkbox = CustomPlotCheckForm()
        rho_checkbox = RhoCheckForm()
        vp_checkbox = VPCheckForm()

    context = {
        "form": form,
        "plot_config": plot_config,
        "plot_checkboxes": [plot_checkbox, rho_checkbox, vp_checkbox],
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
        "custom_plots": custom_plots,
    }

    return render(request, "pred.html", context)


def homepage(request):
    "handle request"
    return render(request, "homepage.html")


def authorpage(request):
    "handle request"
    return render(request, "author.html")


def about(request):
    "handle request"
    return render(request, "about.html")
