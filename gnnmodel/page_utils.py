"""Helpers for page forms, post processing, and context assembly."""

from __future__ import annotations

from typing import Any, Dict, Optional

from .forms import (
    BinaryLLECheckForm,
    BinaryVLECheckForm,
    BinaryVLEpxyCheckForm,
    CustomPlotConfigForm,
    HlvCheckForm,
    InChIorSMILESareaInputforMixture,
    InChIorSMILESinput,
    KijCheckForm,
    MixtureForms,
    PhaseDiagramCheckForm,
    PureForms,
    RhoCheckForm,
    SlvCheckForm,
    STCheckForm,
    TernaryLLECheckForm,
    TernaryVLEpxCheckForm,
    VPCheckForm,
)
from .utils import (
    available_params,
    get_mixture_plots_data,
    get_pred,
    get_pure_plots_data,
)
from .utils_data import (
    retrieve_available_data_binary,
    retrieve_available_data_pure,
    retrieve_available_data_ternary,
    retrieve_vle_for_kij,
)
from .utils_mix import optimize_binary_kij_for_vle


def init_pure_forms(post_data=None):
    """Initialize pure-page forms."""

    if post_data:
        return PureForms(
            InChIorSMILESinput(post_data),
            CustomPlotConfigForm(post_data),
            RhoCheckForm(post_data),
            VPCheckForm(post_data),
            HlvCheckForm(post_data),
            SlvCheckForm(post_data),
            PhaseDiagramCheckForm(post_data),
            STCheckForm(post_data),
        )
    return PureForms(
        InChIorSMILESinput(),
        CustomPlotConfigForm(),
        RhoCheckForm(),
        VPCheckForm(),
        HlvCheckForm(),
        SlvCheckForm(),
        PhaseDiagramCheckForm(),
        STCheckForm(),
    )


def process_pure_post(forms: PureForms) -> Optional[Dict[str, Any]]:
    """Process POST data from the pure page."""

    if not forms.form.is_valid():
        return None

    smiles, inchi = forms.form.cleaned_data["query"]
    phase_diagrams, custom_plots = get_pure_plots_data(smiles=smiles, forms=forms)

    return {
        "smiles": smiles,
        "inchi": inchi,
        "pred": get_pred(smiles),
        "output": True,
        "pure_plots": custom_plots,
        "phase_diagrams": phase_diagrams,
        "available_exp_data": retrieve_available_data_pure(smiles=smiles),
    }


def build_pure_context(forms: PureForms, post_data=None):
    """Build context for the pure component page."""

    data = {
        "form": forms.form,
        "plot_config": forms.plot_config,
        "plot_checkboxes": [
            forms.rho_checkbox,
            forms.vp_checkbox,
            forms.h_lv_checkbox,
            forms.s_lv_checkbox,
            forms.st_checkbox,
            forms.phase_diagram_checkbox,
        ],
        "predicted_para": [(None, None)],
        "mol_identifiers": [(None, None)],
        "output": False,
        "pure_plots": [],
        "phase_diagrams": [],
    }

    if post_data:
        data.update(
            {
                "predicted_para": [
                    (paraname, round(para, 4))
                    for para, paraname in zip(post_data["pred"], available_params)
                ],
                "mol_identifiers": [
                    ("InChI", post_data["inchi"]),
                    ("SMILES", post_data["smiles"]),
                ],
                "output": post_data["output"],
                "pure_plots": post_data["pure_plots"],
                "phase_diagrams": post_data["phase_diagrams"],
                **post_data["available_exp_data"],
            }
        )
    return data


def _build_kij_matrix(
    smiles_list: list[str], kij_values: list[float]
) -> list[list[float]]:
    kij_matrix = [
        [0.0 for _ in range(len(smiles_list))] for _ in range(len(smiles_list))
    ]
    k_idx = 0
    if kij_values:
        for i in range(len(smiles_list)):
            for j in range(i + 1, len(smiles_list)):
                kij_matrix[i][j] = kij_values[k_idx]
                kij_matrix[j][i] = kij_values[k_idx]
                k_idx += 1
    return kij_matrix


def _available_mixture_exp_data(smiles_list: list[str]):
    if len(smiles_list) == 2:
        return retrieve_available_data_binary(smiles_list=smiles_list)
    if len(smiles_list) == 3:
        return retrieve_available_data_ternary(smiles_list=smiles_list)
    return {}


def init_mixture_forms(post_data=None):
    """Initialize mixture-page forms."""

    if post_data:
        return MixtureForms(
            InChIorSMILESareaInputforMixture(post_data),
            CustomPlotConfigForm(post_data),
            RhoCheckForm(post_data),
            VPCheckForm(post_data),
            BinaryVLECheckForm(post_data),
            BinaryLLECheckForm(post_data),
            TernaryLLECheckForm(post_data),
            BinaryVLEpxyCheckForm(post_data),
            TernaryVLEpxCheckForm(post_data),
            KijCheckForm(post_data),
        )
    return MixtureForms(
        InChIorSMILESareaInputforMixture(),
        CustomPlotConfigForm(),
        RhoCheckForm(),
        VPCheckForm(),
        BinaryVLECheckForm(),
        BinaryLLECheckForm(),
        TernaryLLECheckForm(),
        BinaryVLEpxyCheckForm(),
        TernaryVLEpxCheckForm(),
        KijCheckForm(),
    )


def process_mixture_post(forms: MixtureForms):
    """Process POST data from the mixture page."""

    para_pred_list = []
    mole_fractions_list = []
    mixture_plots_ = ([], "", "", "")
    available_exp_data = {}
    kij_value = None

    if forms.form.is_valid():
        _, smiles_list = forms.form.cleaned_data["smiles_inchi_list"]
        mole_fractions_list = forms.form.cleaned_data["mole_fractions"]
        kij = forms.form.cleaned_data["kijs"]
        forms.kij_checkform.full_clean()
        if forms.kij_checkform.cleaned_data["kij_checkbox"]:
            vle = retrieve_vle_for_kij(smiles_list=smiles_list)
            if vle is not None:
                forms.plot_config.full_clean()
                kij_value = optimize_binary_kij_for_vle(
                    smiles_list=smiles_list,
                    initial_kij=kij[0] if kij else 0.05,
                    vle=vle,
                    npoints=forms.plot_config.cleaned_data["npoints"] or 5,
                )
                kij = [kij_value]
        kij_matrix = _build_kij_matrix(smiles_list, kij)
        para_pred_list = [
            [round(para, 5) for para in get_pred(smiles)] for smiles in smiles_list
        ]
        mixture_plots_ = get_mixture_plots_data(
            smiles_list=smiles_list,
            mole_fractions_list=mole_fractions_list,
            forms=forms,
            kij_matrix=kij_matrix,
        )
        available_exp_data = _available_mixture_exp_data(smiles_list)

    return {
        "form": forms.form,
        "plot_config": forms.plot_config,
        "plot_checkboxes": [
            forms.rho_checkform,
            forms.vp_checkform,
            forms.binary_vle_checkform,
            forms.binary_vlepxy_checkform,
            forms.binary_lle_checkform,
            forms.ternary_lle_checkform,
            forms.ternary_vlepx_checkform,
            forms.kij_checkform,
        ],
        "para_pred_list": para_pred_list,
        "mole_fractions_list": mole_fractions_list,
        "mixture_plots": mixture_plots_,
        "output": bool(para_pred_list),
        "kij_value": kij_value,
        "available_exp_data": available_exp_data,
    }


def build_mixture_context(post_data=None):
    """Build context for the mixture page."""

    if post_data:
        return {
            "form": post_data["form"],
            "plot_config": post_data["plot_config"],
            "plot_checkboxes": post_data["plot_checkboxes"],
            "available_params": available_params,
            "parameters_molefractions_list": list(
                zip(post_data["para_pred_list"], post_data["mole_fractions_list"])
            ),
            "mixture_plots": post_data["mixture_plots"][0],
            "binary_lle_phase_diagram_data": post_data["mixture_plots"][1],
            "vle_phase_diagram_data": post_data["mixture_plots"][2],
            "ternary_lle_phase_diagram_data": post_data["mixture_plots"][3],
            "output": post_data["output"],
            "kij_value": post_data["kij_value"],
            **post_data["available_exp_data"],
        }

    return {
        "form": InChIorSMILESareaInputforMixture(),
        "plot_config": CustomPlotConfigForm(),
        "plot_checkboxes": [
            RhoCheckForm(),
            VPCheckForm(),
            BinaryVLECheckForm(),
            BinaryVLEpxyCheckForm(),
            BinaryLLECheckForm(),
            TernaryLLECheckForm(),
            TernaryVLEpxCheckForm(),
        ],
        "available_params": available_params,
        "parameters_molefractions_list": [],
        "mixture_plots": [],
        "binary_lle_phase_diagram_data": "",
        "vle_phase_diagram_data": "",
        "ternary_lle_phase_diagram_data": "",
        "output": False,
        "kij_value": None,
    }
