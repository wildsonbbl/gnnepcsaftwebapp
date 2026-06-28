"Module for utils like plotting data."

import json
import os.path as osp
from typing import Any, Dict, List, Optional, Tuple, Union

import onnxruntime as ort
from django.conf import settings
from gnnepcsaft.pcsaft.pcsaft_feos import critical_points_feos
from gnnepcsaft_mcp_server.utils import predict_pcsaft_parameters
from rdkit.Chem import AllChem as Chem

from . import logger
from .forms import (
    BinaryLLECheckForm,
    BinaryVLECheckForm,
    CustomPlotCheckForm,
    CustomPlotConfigForm,
    HlvCheckForm,
    InChIorSMILESareaInputforMixture,
    InChIorSMILESinput,
    PhaseDiagramCheckForm,
    RhoCheckForm,
    SlvCheckForm,
    STCheckForm,
    TernaryLLECheckForm,
    VPCheckForm,
)
from .utils_data import (
    retrieve_available_data_binary,
    retrieve_available_data_pure,
    retrieve_available_data_ternary,
    retrieve_bubble_pressure_data,
    retrieve_lle_binary_data,
    retrieve_lle_ternary_data,
    retrieve_rho_binary_data,
    retrieve_rho_pure_data,
    retrieve_rho_ternary_data,
    retrieve_st_pure_data,
    retrieve_vle_binary_data,
    retrieve_vle_pxy_binary_data,
    retrieve_vle_ternary_data,
    retrieve_vle_ternary_tx_fixed_data,
    retrieve_vp_pure_data,
)
from .utils_mix import (
    MixDenParams,
    MixLLEParams,
    MixVpParams,
    mix_den,
    mix_lle,
    mix_ternary_lle,
    mix_vle,
    mix_vp,
)
from .utils_pure import (
    pure_den,
    pure_h_lv,
    pure_phase_diagram,
    pure_s_lv,
    pure_surface_tension,
    pure_vp,
)

# lazy import
# import polars as pl

ort.set_default_logger_severity(3)

available_params = [
    "Segment number",
    "Segment diameter (Å)",
    "Dispersion energy (K)",
    "Association volume",
    "Association energy (K)",
    "Dipole moment (D)*",
    "Number of association site A",
    "Number of association site B",
    "Molecular weight (g/mol)",
    "Critical temperature (K)",
    "Critical pressure (Bar)",
]

file_dir = osp.dirname(__file__)
dataset_dir = osp.join(file_dir, "data")


def plotmol(inchi: str) -> str:
    "Make Mol block for 3Dmol."

    mol = Chem.MolFromInchi(inchi)
    mol = Chem.AddHs(mol)  # type: ignore
    params = Chem.ETKDGv3()  # type: ignore
    params.randomSeed = 0xF00D
    result = Chem.EmbedMolecule(mol, params)  # type: ignore
    if result == 0:
        Chem.MMFFOptimizeMolecule(  # type: ignore
            mol, maxIters=1000, nonBondedThresh=100, ignoreInterfragInteractions=False
        )
    # mol = Chem.RemoveHs(mol, implicitOnly=False)
    imgmol = Chem.MolToV3KMolBlock(mol)  # type: ignore
    return imgmol


def pure_plots(
    smiles: str,
    temp_min: float,
    temp_max: float,
    pressure: float,
    selected_checkboxes: Optional[List[str]] = None,
    npoints: int = 10,
) -> Union[List[tuple[str, str, str, str, str]], List]:
    """
    Pure plots.

    args
    ---------
    smiles: str
      SMILES
    temp_min: float
      minimum temperature in Kelvin
    temp_max: float
      maximum temperature in Kelvin
    pressure: float
      pressure in Pa
    checkboxes: list
      list with ids to plot
    npoints: int
      number of data points
    """

    if selected_checkboxes is None:
        selected_checkboxes = [
            "den_plot",
            "vp_plot",
            "h_lv_plot",
            "s_lv_plot",
            "st_plot",
        ]

    all_plots = []

    if (
        "den_plot" in selected_checkboxes
        and temp_min is not None
        and temp_max is not None
        and pressure is not None
    ):
        plot_data = {}
        try:
            plot_data["GNN"] = pure_den(
                smiles=smiles,
                min_temp=temp_min,
                max_temp=temp_max,
                pressure=pressure,
                npoints=npoints,
            )
        except (AssertionError, RuntimeError) as e:
            logger.debug(e)
        try:
            exp_data = retrieve_rho_pure_data(smiles=smiles, pressure=pressure / 1000)
            if exp_data is not None:
                plot_data["TML"] = exp_data.T.tolist()
            else:
                plot_data["TML"] = ([], [])
        except (AssertionError, RuntimeError) as e:
            logger.debug(e)
        all_plots.append(
            (
                json.dumps(plot_data),
                "Temperature (K)",
                "Density (mol / m³)",
                f"Density at {pressure} Pa for {smiles}",
                "den_plot",
            )
        )

    if (
        "vp_plot" in selected_checkboxes
        and temp_min is not None
        and temp_max is not None
    ):
        plot_data = {}
        try:
            plot_data["GNN"] = pure_vp(
                smiles=smiles,
                min_temp=temp_min,
                max_temp=temp_max,
                npoints=npoints,
            )
        except (AssertionError, RuntimeError) as e:
            logger.debug(e)
        try:
            exp_data = retrieve_vp_pure_data(smiles=smiles)
            if exp_data is not None:
                exp_data[:, 1] *= 1000
                plot_data["TML"] = exp_data.T.tolist()
            else:
                plot_data["TML"] = ([], [])
        except (AssertionError, RuntimeError) as e:
            logger.debug(e)
        all_plots.append(
            (
                json.dumps(plot_data),
                "Temperature (K)",
                "Vapor pressure (Pa)",
                "Vapor pressure",
                "vp_plot",
            )
        )

    if (
        "h_lv_plot" in selected_checkboxes
        and temp_min is not None
        and temp_max is not None
    ):
        plot_data = {}
        try:
            plot_data["GNN"] = pure_h_lv(
                smiles=smiles,
                min_temp=temp_min,
                max_temp=temp_max,
                npoints=npoints,
            )
        except (AssertionError, RuntimeError) as e:
            logger.debug(e)
        plot_data["TML"] = ([], [])
        all_plots.append(
            (
                json.dumps(plot_data),
                "Temperature (K)",
                "Enthalpy of vaporization (kJ/mol)",
                "Enthalpy of vaporization",
                "h_lv_plot",
            )
        )

    if (
        "s_lv_plot" in selected_checkboxes
        and temp_min is not None
        and temp_max is not None
    ):
        plot_data = {}
        try:
            plot_data["GNN"] = pure_s_lv(
                smiles=smiles,
                min_temp=temp_min,
                max_temp=temp_max,
                npoints=npoints,
            )
        except (AssertionError, RuntimeError) as e:
            logger.debug(e)
        plot_data["TML"] = ([], [])
        all_plots.append(
            (
                json.dumps(plot_data),
                "Temperature (K)",
                "Entropy of vaporization (J/mol/K)",
                "Entropy of vaporization",
                "s_lv_plot",
            )
        )

    if "st_plot" in selected_checkboxes and temp_min is not None:
        plot_data = {}
        try:
            plot_data["GNN"] = pure_surface_tension(
                smiles=smiles,
                min_temp=temp_min,
            )
        except (AssertionError, RuntimeError) as e:
            logger.debug(e)
        try:
            exp_data = retrieve_st_pure_data(smiles=smiles)
            if exp_data is not None:
                exp_data[:, 1] *= 1000
                plot_data["TML"] = exp_data.T.tolist()
            else:
                plot_data["TML"] = ([], [])
        except (AssertionError, RuntimeError) as e:
            logger.debug(e)
        all_plots.append(
            (
                json.dumps(plot_data),
                "Temperature (K)",
                "Surface Tension (mN/m)",
                "Surface Tension",
                "st_plot",
            )
        )

    return all_plots


def get_pred(smiles: str) -> List[float]:
    "get prediction"
    pred = predict_pcsaft_parameters(smiles)

    try:
        critical_points = critical_points_feos(pred.copy())
    except RuntimeError:
        critical_points = [0.0, 0.0]

    pred.append(critical_points[0])
    pred.append(critical_points[1] * 0.00001)  # convert from Pa to Bar
    return pred


def get_forms(request):
    "get forms"
    return (
        InChIorSMILESinput(request.POST),
        CustomPlotConfigForm(request.POST),
        CustomPlotCheckForm(request.POST),
        RhoCheckForm(request.POST),
        VPCheckForm(request.POST),
        HlvCheckForm(request.POST),
        SlvCheckForm(request.POST),
        PhaseDiagramCheckForm(request.POST),
        STCheckForm(request.POST),
    )


def get_pure_plots_data(
    smiles: str,
    plot_config: CustomPlotConfigForm,
    checkboxes: tuple[
        RhoCheckForm,
        VPCheckForm,
        HlvCheckForm,
        SlvCheckForm,
        PhaseDiagramCheckForm,
        STCheckForm,
    ],
) -> Tuple[List[List[float]], List]:
    "get custom plots data for pure component"

    (
        rho_checkbox_,
        vp_checkbox_,
        h_lv_checkbox_,
        s_lv_checkbox_,
        phase_diagram_checkbox_,
        st_checkbox_,
    ) = checkboxes
    plot_config.full_clean()
    rho_checkbox_.full_clean()
    vp_checkbox_.full_clean()
    h_lv_checkbox_.full_clean()
    s_lv_checkbox_.full_clean()
    phase_diagram_checkbox_.full_clean()
    st_checkbox_.full_clean()
    custom_plots = []
    phase_diagrams = []

    selected_checkboxes = []
    if rho_checkbox_.cleaned_data["rho_checkbox"]:
        selected_checkboxes.append("den_plot")
    if vp_checkbox_.cleaned_data["vp_checkbox"]:
        selected_checkboxes.append("vp_plot")
    if h_lv_checkbox_.cleaned_data["h_lv_checkbox"]:
        selected_checkboxes.append("h_lv_plot")
    if s_lv_checkbox_.cleaned_data["s_lv_checkbox"]:
        selected_checkboxes.append("s_lv_plot")
    if st_checkbox_.cleaned_data["st_checkbox"]:
        selected_checkboxes.append("st_plot")

    try:
        custom_plots = pure_plots(
            smiles=smiles,
            temp_min=plot_config.cleaned_data["temp_min"],
            temp_max=plot_config.cleaned_data["temp_max"],
            pressure=plot_config.cleaned_data["pressure"],
            selected_checkboxes=selected_checkboxes,
        )
    except RuntimeError as err:
        logger.debug(err)

    if phase_diagram_checkbox_.cleaned_data["phase_diagram_checkbox"]:
        try:
            if plot_config.cleaned_data["temp_min"] is not None:
                phase_diagrams = pure_phase_diagram(
                    smiles=smiles, min_temp=plot_config.cleaned_data["temp_min"]
                )
        except RuntimeError as err:
            logger.debug(err)
    return phase_diagrams, custom_plots


def get_mixture_plots_data(
    smiles_list: List[str],
    mole_fractions_list: List[float],
    config: Tuple[
        CustomPlotConfigForm,
        TernaryLLECheckForm,
        BinaryVLECheckForm,
        BinaryLLECheckForm,
    ],
    kij_matrix: List[List[float]],
) -> Tuple[List[Tuple[str, str, str, str, str]], str, str, str]:
    "get mixture plots data"

    plot_config, ternary_lle_checkform, binary_vle_checkform, binary_lle_checkform = (
        config
    )

    plot_config.full_clean()
    ternary_lle_checkform.full_clean()
    binary_vle_checkform.full_clean()
    binary_lle_checkform.full_clean()
    mixture_plot = []
    ternary_lle_phase_diagram_data = ""
    binary_lle_phase_diagram_data = ""
    vle_phase_diagram_data = ""

    mixture_plot = mixture_plots(
        smiles_list=smiles_list,
        state_list=(
            mole_fractions_list,
            plot_config.cleaned_data["temp_min"],
            plot_config.cleaned_data["temp_max"],
            plot_config.cleaned_data["pressure"],
        ),
        kij_matrix=kij_matrix,
    )

    try:
        if ternary_lle_checkform.cleaned_data["ternary_lle_checkbox"] is False:
            raise ValueError("Ternary LLE checkbox not selected.")
        if len(smiles_list) != 3:
            raise ValueError("LLE/VLE phase diagram only for ternary mixtures.")
        if plot_config.cleaned_data["temp_min"] is None:
            raise ValueError("Missing minimum temperature")
        if plot_config.cleaned_data["pressure"] is None:
            raise ValueError("Missing pressure")

        _ternary_lle_phase_diagram_data = mix_ternary_lle(
            smiles_list=smiles_list,
            kij_matrix=kij_matrix,
            temperature=plot_config.cleaned_data["temp_min"],
            pressure=plot_config.cleaned_data["pressure"],
            npoints=20,
        )
        ternary_lle_phase_diagram_data = json.dumps(_ternary_lle_phase_diagram_data)
    except (ValueError, RuntimeError) as err:
        logger.debug(err)

    try:
        if binary_lle_checkform.cleaned_data["binary_lle_checkbox"] is False:
            raise ValueError("Binary LLE checkbox not selected.")
        if len(smiles_list) != 2:
            raise ValueError("LLE phase diagram only for binary mixtures.")
        if plot_config.cleaned_data["temp_min"] is None:
            raise ValueError("Missing minimum temperature")
        if plot_config.cleaned_data["pressure"] is None:
            raise ValueError("Missing pressure")

        _binary_lle_phase_diagram_data = mix_lle(
            MixLLEParams(
                smiles_list=smiles_list,
                mole_fractions=mole_fractions_list,
                kij_matrix=kij_matrix,
                temperature=plot_config.cleaned_data["temp_min"],
                pressure=plot_config.cleaned_data["pressure"],
                npoints=10,
            )
        )
        binary_lle_phase_diagram_data = json.dumps(_binary_lle_phase_diagram_data)

    except (ValueError, RuntimeError) as err:
        logger.debug(err)

    try:
        if binary_vle_checkform.cleaned_data["binary_vle_checkbox"] is False:
            raise ValueError("Binary VLE checkbox not selected.")
        if len(smiles_list) != 2:
            raise ValueError("VLE phase diagram only for binary mixtures.")
        if plot_config.cleaned_data["pressure"] is None:
            raise ValueError("Missing pressure")

        _vle_phase_diagram_data = mix_vle(
            smiles_list=smiles_list,
            kij_matrix=kij_matrix,
            pressure=plot_config.cleaned_data["pressure"],
            npoints=20,
        )
        vle_phase_diagram_data = json.dumps(_vle_phase_diagram_data)
    except (ValueError, RuntimeError) as err:
        logger.debug(err)

    return (
        mixture_plot,
        binary_lle_phase_diagram_data,
        vle_phase_diagram_data,
        ternary_lle_phase_diagram_data,
    )


def mixture_plots(
    smiles_list: List[str],
    state_list: Tuple[List[float], float, float, float],
    kij_matrix: List[List[float]],
) -> List[Tuple[str, str, str, str, str]]:
    "get mixture plots data"
    mole_fractions_list, temp_min, temp_max, pressure = state_list
    all_plots = []

    if temp_min is not None and temp_max is not None and pressure is not None:
        plot_data = {}
        try:
            plot_data["GNN"] = mix_den(
                MixDenParams(
                    smiles_list=smiles_list,
                    mole_fractions=mole_fractions_list,
                    kij_matrix=kij_matrix,
                    min_temp=temp_min,
                    max_temp=temp_max,
                    pressure=pressure,
                    npoints=20,
                )
            )
        except (AssertionError, RuntimeError) as e:
            logger.debug(e)
        try:
            plot_data["TML"] = ([], [])
            if len(smiles_list) == 2:
                exp_data = retrieve_rho_binary_data(
                    smiles_list=smiles_list,
                    pressure=pressure / 1000,
                    x1=mole_fractions_list[0],
                )
                if exp_data is not None:
                    plot_data["TML"] = exp_data.T.tolist()

            elif len(smiles_list) == 3:
                exp_data = retrieve_rho_ternary_data(
                    smiles_list=smiles_list,
                    pressure=pressure / 1000,
                    x1=mole_fractions_list[0],
                    x2=mole_fractions_list[1],
                )
                if exp_data is not None:
                    plot_data["TML"] = exp_data.T.tolist()

        except (AssertionError, RuntimeError) as e:
            logger.debug(e)

        all_plots.append(
            (
                json.dumps(plot_data),
                "Temperature (K)",
                "Liquid Density (mol / m³)",
                f"Density at {pressure} Pa",
                "mix_den_plot",
            )
        )

    if temp_min is not None and temp_max is not None:
        plot_data = {}
        try:
            plot_data["GNN"] = mix_vp(
                MixVpParams(
                    smiles_list=smiles_list,
                    mole_fractions=mole_fractions_list,
                    kij_matrix=kij_matrix,
                    min_temp=temp_min,
                    max_temp=temp_max,
                    npoints=20,
                )
            )
        except (AssertionError, RuntimeError) as e:
            logger.debug(e)
        try:
            plot_data["TML"] = ([], [])
            if len(smiles_list) == 2:
                exp_data = retrieve_bubble_pressure_data(
                    smiles_list=smiles_list,
                    x1=mole_fractions_list[0],
                )
                if exp_data is not None:
                    exp_data[:, 1] *= 1000
                    plot_data["TML"] = exp_data.T.tolist()
        except (AssertionError, RuntimeError) as e:
            logger.debug(e)

        all_plots.append(
            (
                json.dumps(plot_data),
                "Temperature (K)",
                "Pressure (Pa)",
                "VLE",
                "mix_vp_plot",
            )
        )

    return all_plots


def init_pure_forms(post_data=None):
    "init pure forms"
    if post_data:
        return (
            InChIorSMILESinput(post_data),
            CustomPlotConfigForm(post_data),
            CustomPlotCheckForm(post_data),
            RhoCheckForm(post_data),
            VPCheckForm(post_data),
            HlvCheckForm(post_data),
            SlvCheckForm(post_data),
            PhaseDiagramCheckForm(post_data),
            STCheckForm(post_data),
        )
    return (
        InChIorSMILESinput(),
        CustomPlotConfigForm(),
        CustomPlotCheckForm(),
        RhoCheckForm(),
        VPCheckForm(),
        HlvCheckForm(),
        SlvCheckForm(),
        PhaseDiagramCheckForm(),
        STCheckForm(),
    )


def process_pure_post(
    forms: Tuple[
        InChIorSMILESinput,
        CustomPlotConfigForm,
        CustomPlotCheckForm,
        RhoCheckForm,
        VPCheckForm,
        HlvCheckForm,
        SlvCheckForm,
        PhaseDiagramCheckForm,
        STCheckForm,
    ],
) -> Optional[Dict[str, Any]]:
    "process the post data from the pure page"
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
    ) = forms

    if not form.is_valid():
        return None

    smiles, inchi = form.cleaned_data["query"]
    pred = get_pred(smiles)
    output = True

    plot_checkbox.full_clean()
    phase_diagrams, custom_plots = [], []
    if plot_checkbox.cleaned_data["custom_plot_checkbox"]:
        phase_diagrams, custom_plots = get_pure_plots_data(
            smiles=smiles,
            plot_config=plot_config,
            checkboxes=(
                rho_checkbox,
                vp_checkbox,
                h_lv_checkbox,
                s_lv_checkbox,
                phase_diagram_checkbox,
                st_checkbox,
            ),
        )
    return {
        "smiles": smiles,
        "inchi": inchi,
        "pred": pred,
        "output": output,
        "pure_plots": custom_plots,
        "phase_diagrams": phase_diagrams,
    }


def build_pure_context(forms, post_data=None):
    "build context for pure component page"
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
    ) = forms

    data = {
        "form": form,
        "plot_config": plot_config,
        "plot_checkboxes": [
            plot_checkbox,
            rho_checkbox,
            vp_checkbox,
            h_lv_checkbox,
            s_lv_checkbox,
            st_checkbox if settings.PLATFORM == "desktop" else None,
            phase_diagram_checkbox,
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
            }
        )
    return data


def init_mixture_forms(post_data=None):
    "init mixture forms"
    if post_data:
        return (
            InChIorSMILESareaInputforMixture(post_data),
            CustomPlotConfigForm(post_data),
            BinaryVLECheckForm(post_data),
            BinaryLLECheckForm(post_data),
            TernaryLLECheckForm(post_data),
        )
    return (
        InChIorSMILESareaInputforMixture(),
        CustomPlotConfigForm(),
        BinaryVLECheckForm(),
        BinaryLLECheckForm(),
        TernaryLLECheckForm(),
    )


def process_mixture_post(
    forms: Tuple[
        InChIorSMILESareaInputforMixture,
        CustomPlotConfigForm,
        BinaryVLECheckForm,
        BinaryLLECheckForm,
        TernaryLLECheckForm,
    ],
):
    "process the post data from the mixture page"
    (
        form,
        plot_config,
        binary_vle_checkform,
        binary_lle_checkform,
        ternary_lle_checkform,
    ) = forms
    para_pred_list = []
    mole_fractions_list = []
    mixture_plots_ = (([], []), "", "", "")
    output = False

    if form.is_valid():
        _, smiles_list, mole_fractions_list, kij = form.cleaned_data["text_area"]
        kij_matrix = [
            [0.0 for _ in range(len(smiles_list))] for _ in range(len(smiles_list))
        ]
        # fill kij_matrix
        k_idx = 0
        for i in range(len(smiles_list)):
            for j in range(i + 1, len(smiles_list)):
                kij_matrix[i][j] = kij[k_idx]
                kij_matrix[j][i] = kij[k_idx]
                k_idx += 1
        for smiles in smiles_list:
            para_pred = [round(para, 5) for para in get_pred(smiles)]
            para_pred_list.append(para_pred)
        mixture_plots_ = get_mixture_plots_data(
            smiles_list=smiles_list,
            mole_fractions_list=mole_fractions_list,
            config=(
                plot_config,
                ternary_lle_checkform,
                binary_vle_checkform,
                binary_lle_checkform,
            ),
            kij_matrix=kij_matrix,
        )
        output = True

    return {
        "form": form,
        "binary_vle_checkform": binary_vle_checkform,
        "binary_lle_checkform": binary_lle_checkform,
        "ternary_lle_checkform": ternary_lle_checkform,
        "plot_config": plot_config,
        "para_pred_list": para_pred_list,
        "mole_fractions_list": mole_fractions_list,
        "mixture_plots": mixture_plots_,
        "output": output,
    }


def build_mixture_context(post_data=None):
    "build context for mixture page"
    if post_data:
        return {
            "form": post_data["form"],
            "plot_config": post_data["plot_config"],
            "binary_vle_checkform": post_data["binary_vle_checkform"],
            "binary_lle_checkform": post_data["binary_lle_checkform"],
            "ternary_lle_checkform": post_data["ternary_lle_checkform"],
            "available_params": available_params,
            "parameters_molefractions_list": list(
                zip(post_data["para_pred_list"], post_data["mole_fractions_list"])
            ),
            "mixture_plots": post_data["mixture_plots"][0],
            "binary_lle_phase_diagram_data": post_data["mixture_plots"][1],
            "vle_phase_diagram_data": post_data["mixture_plots"][2],
            "ternary_lle_phase_diagram_data": post_data["mixture_plots"][3],
            "output": post_data["output"],
        }
    return {
        "form": InChIorSMILESareaInputforMixture(),
        "binary_vle_checkform": BinaryVLECheckForm(),
        "binary_lle_checkform": BinaryLLECheckForm(),
        "ternary_lle_checkform": TernaryLLECheckForm(),
        "plot_config": CustomPlotConfigForm(),
        "available_params": available_params,
        "parameters_molefractions_list": [],
        "mixture_plots": [],
        "vp_plots": [],
        "binary_lle_phase_diagram_data": "",
        "vle_phase_diagram_data": "",
        "output": False,
    }
