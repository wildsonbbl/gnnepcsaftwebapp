"request handler."

import os.path as osp

from django.conf import settings
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
    "Critical temperature (K)",
    "Critical pressure (Bar)",
]


# Create your views here.
def estimator(request):
    "handle request"

    pred = None
    query = ""
    output = False
    plotden, plotvp, molimg, phase_diagrams = "", "", "", ""
    custom_plots = []
    if request.method == "POST":
        form = InChIorSMILESinput(request.POST)
        plot_config = CustomPlotConfigForm(request.POST)
        plot_checkbox = CustomPlotCheckForm(request.POST)
        rho_checkbox = RhoCheckForm(request.POST)
        vp_checkbox = VPCheckForm(request.POST)
        h_lv_checkbox = HlvCheckForm(request.POST)
        s_lv_checkbox = SlvCheckForm(request.POST)
        phase_diagram_checkbox = PhaseDiagramCheckForm(request.POST)

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
            try:
                critical_points = critical_points_feos(pred)
            except RuntimeError:
                critical_points = [0.0, 0.0]
            pred.append(critical_points[0])
            pred.append(critical_points[1] * 0.00001)  # convert from Pa to Bar
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

            if plot_checkbox.cleaned_data["custom_plot_checkbox"]:
                rho_checkbox.full_clean()
                vp_checkbox.full_clean()
                h_lv_checkbox.full_clean()
                s_lv_checkbox.full_clean()
                phase_diagram_checkbox.full_clean()
                custom_plots = custom_plot(
                    pred,
                    plot_config.cleaned_data["temp_min"],
                    plot_config.cleaned_data["temp_max"],
                    plot_config.cleaned_data["pressure"],
                    [
                        rho_checkbox.cleaned_data["rho_checkbox"],
                        vp_checkbox.cleaned_data["vp_checkbox"],
                        h_lv_checkbox.cleaned_data["h_lv_checkbox"],
                        s_lv_checkbox.cleaned_data["s_lv_checkbox"],
                    ],
                )
                if phase_diagram_checkbox.cleaned_data["phase_diagram_checkbox"]:
                    try:
                        phase_diagrams_all_data = phase_diagram_feos(
                            pred, [plot_config.cleaned_data["temp_min"]]
                        )
                        phase_diagrams = [
                            phase_diagrams_all_data.get("temperature", [0]),
                            phase_diagrams_all_data.get(
                                "pressure",
                                phase_diagrams_all_data.get("pressure vapor", [0]),
                            ),
                            phase_diagrams_all_data["density liquid"],
                            phase_diagrams_all_data["density vapor"],
                        ]
                    except RuntimeError:
                        phase_diagrams = ""

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
        "query": query,
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
