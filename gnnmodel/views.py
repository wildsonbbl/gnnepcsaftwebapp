"request handler."

import os.path as osp

from django.shortcuts import render
from gnnepcsaft.epcsaft.epcsaft_feos import critical_points_feos, phase_diagram_feos

from .forms import (
    CustomPlotCheckForm,
    CustomPlotConfigForm,
    HlvCheckForm,
    InChIorSMILESinput,
    PhaseDiagramCheckForm,
    RhoCheckForm,
    SlvCheckForm,
    VPCheckForm,
)
from .models import GnnepcsaftPara, ThermoMLDenData, ThermoMLVPData
from .utils import get_custom_plots_data, get_forms, get_main_plots_data, get_pred

file_dir = osp.dirname(__file__)

available_params = [
    "Segment number",
    "Segment diameter (Å)",
    "Dispersion energy (K)",
    "Association volume",
    "Association energy (K)",
    "Dipole moment (D)*",
    "Number of association site A",
    "Number of association site B",
    "Critical temperature (K)",
    "Critical pressure (Bar)",
]


# Create your views here.
def estimator(request):  # pylint: disable=R0914
    "handle request"

    pred = []
    output = False
    plotden, plotvp, molimg = ("", "", "")
    custom_plots, phase_diagrams = [], []
    inchi = ""
    smiles = ""
    if request.method == "POST":
        (
            form,
            plot_config,
            plot_checkbox,
            rho_checkbox,
            vp_checkbox,
            h_lv_checkbox,
            s_lv_checkbox,
            phase_diagram_checkbox,
        ) = get_forms(request)

        if form.is_valid():
            smiles, inchi = form.cleaned_data["query"]

            pred = get_pred(smiles, inchi)
            plotden, plotvp, molimg = get_main_plots_data(inchi)
            output = True

            plot_checkbox.full_clean()
            if plot_checkbox.cleaned_data["custom_plot_checkbox"]:
                phase_diagrams, custom_plots = get_custom_plots_data(
                    pred[:-2],
                    plot_config,
                    (
                        rho_checkbox,
                        vp_checkbox,
                        h_lv_checkbox,
                        s_lv_checkbox,
                        phase_diagram_checkbox,
                    ),
                )

    else:
        form = InChIorSMILESinput()
        plot_config = CustomPlotConfigForm()
        plot_checkbox = CustomPlotCheckForm()
        rho_checkbox = RhoCheckForm()
        vp_checkbox = VPCheckForm()
        h_lv_checkbox = HlvCheckForm()
        s_lv_checkbox = SlvCheckForm()
        phase_diagram_checkbox = PhaseDiagramCheckForm()

    context = {
        "form": form,
        "plot_config": plot_config,
        "plot_checkboxes": [
            plot_checkbox,
            rho_checkbox,
            vp_checkbox,
            h_lv_checkbox,
            s_lv_checkbox,
            phase_diagram_checkbox,
        ],
        "predicted_para": (
            [
                (paraname, round(para, 4))
                for para, paraname in zip(pred, available_params)
            ]
            if output
            else [(None, None)]
        ),
        "mol_identifiers": (
            [("InChI", inchi), ("SMILES", smiles)] if output else [(None, None)]
        ),
        "output": output,
        "plotden": plotden != "",
        "plotvp": plotvp != "",
        "den_data": plotden,
        "vp_data": plotvp,
        "mol_data": molimg,
        "custom_plots": custom_plots,
        "phase_diagrams": phase_diagrams,
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
