"request handler."

import os.path as osp

from django.shortcuts import render
from gnnepcsaft.data.rdkit_util import mw

from .forms import (
    CustomPlotCheckForm,
    CustomPlotConfigForm,
    GoogleAPIKeyForm,
    HlvCheckForm,
    InChIorSMILESareaInput,
    InChIorSMILESareaInputforMixture,
    InChIorSMILESinput,
    PhaseDiagramCheckForm,
    RhoCheckForm,
    SlvCheckForm,
    STCheckForm,
    VPCheckForm,
)
from .utils import (
    get_custom_plots_data,
    get_forms,
    get_main_plots_data,
    get_mixture_plots_data,
    get_pred,
)
from .utils_llm import resume_mol

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
def pure(request):  # pylint: disable=R0914
    "handle request for pure substance"

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
            st_checkbox,
            google_api_key_form,
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
                        st_checkbox,
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
        st_checkbox = STCheckForm()
        google_api_key_form = GoogleAPIKeyForm()

    context = {
        "form": form,
        "plot_config": plot_config,
        "plot_checkboxes": [
            plot_checkbox,
            rho_checkbox,
            vp_checkbox,
            h_lv_checkbox,
            s_lv_checkbox,
            st_checkbox,
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
        "google_api_key_form": google_api_key_form,
    }

    return render(request, "pure.html", context)


def batch(request):
    "handle request for batch of substances"
    pred_list = []
    output = False
    if request.method == "POST":
        form = InChIorSMILESareaInput(request.POST)
        if form.is_valid():
            inchi_list, smiles_list = form.cleaned_data["text_area"]

            for smiles, inchi in zip(smiles_list, inchi_list):
                pred_list.append([round(para, 5) for para in get_pred(smiles, inchi)])
            output = True
    else:
        form = InChIorSMILESareaInput()

    context = {
        "form": form,
        "available_params": available_params,
        "parameters_list": pred_list,
        "output": output,
    }

    return render(request, "batch.html", context)


def mixture(request):
    "handle request for mixture"
    para_pred_list = []
    output = False
    mixture_plots = ([], [])
    mole_fractions_list = []
    para_pred_for_plot = []
    if request.method == "POST":
        form = InChIorSMILESareaInputforMixture(request.POST)
        plot_config = CustomPlotConfigForm(request.POST)
        if form.is_valid():
            inchi_list, smiles_list, mole_fractions_list = form.cleaned_data[
                "text_area"
            ]
            for smiles, inchi in zip(smiles_list, inchi_list):
                para_pred_list.append(
                    [round(para, 5) for para in get_pred(smiles, inchi)]
                )
                para_pred_for_plot.append(para_pred_list[-1] + [mw(inchi)])
            mixture_plots = get_mixture_plots_data(
                para_pred_for_plot, mole_fractions_list, plot_config
            )
            output = True
    else:
        form = InChIorSMILESareaInputforMixture()
        plot_config = CustomPlotConfigForm()

    context = {
        "form": form,
        "plot_config": plot_config,
        "available_params": available_params,
        "parameters_molefractions_list": list(zip(para_pred_list, mole_fractions_list)),
        "mixture_plots": mixture_plots[0],
        "vp_plots": mixture_plots[1],
        "output": output,
    }

    return render(request, "mixture.html", context)


def homepage(request):
    "handle request"
    return render(request, "homepage.html")


def authorpage(request):
    "handle request"
    return render(request, "author.html")


def about(request):
    "handle request for about page"
    return render(request, "about.html")


def description(request):
    "handle request for molecule description"

    html_output = ""

    if request.method == "POST":
        form = InChIorSMILESinput(request.POST)
        google_api_key_form = GoogleAPIKeyForm(request.POST)
        if google_api_key_form.is_valid() and form.is_valid():
            smiles, inchi = form.cleaned_data["query"]
            html_output = resume_mol(
                inchi, smiles, google_api_key_form.cleaned_data["google_api_key"]
            )
    else:
        form = InChIorSMILESinput()
        google_api_key_form = GoogleAPIKeyForm()

    return render(
        request,
        "description.html",
        {
            "output": html_output,
            "google_api_key_form": google_api_key_form,
        },
    )
