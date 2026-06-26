"Module for utils like plotting data."

import json
import os.path as osp
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import onnxruntime as ort
import polars as pl
from django.conf import settings
from gnnepcsaft.pcsaft.pcsaft_feos import (
    critical_points_feos,
    mix_den_feos,
    mix_lle_diagram_feos,
    mix_lle_feos,
    mix_vle_diagram_feos,
    mix_vp_feos,
    phase_diagram_feos,
    pure_den_feos,
    pure_h_lv_feos,
    pure_s_lv_feos,
    pure_surface_tension_feos,
    pure_vp_feos,
)
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


# pylint: disable=R0914
def plotdata(
    para: List[float], rho: np.ndarray, vp: np.ndarray
) -> Tuple[Dict[str, List[float]], Dict[str, List[float]]]:
    "Organize data for plotting."
    plotden, plotvp = {}, {}
    pred_rho, pred_vp = rhovp_data(para, rho, vp)

    if rho.shape[0] > 2:

        plotden = {
            "T": rho[:, 0].tolist(),
            "TML": rho[:, -1].tolist(),
            "GNN": pred_rho,
        }
    # plot vp data
    if vp.shape[0] > 2:
        plotvp = {
            "T": vp[:, 0].tolist(),
            "TML": (vp[:, -1]).tolist(),
            "GNN": pred_vp,
        }

    return plotden, plotvp


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


def rhovp_data(
    parameters: List[float], rho: np.ndarray, vp: np.ndarray
) -> Tuple[List[float], List[float]]:
    """Calculates density and vapor pressure with PC-SAFT"""

    all_pred_den = []
    if rho.shape[0] > 0:
        for state in rho:
            try:
                all_pred_den += [pure_den_feos(parameters, state)]
            except (AssertionError, RuntimeError):
                all_pred_den += [float("nan")]

    all_pred_vp = []
    if vp.shape[0] > 0:
        for state in vp:
            try:
                all_pred_vp += [pure_vp_feos(parameters, state)]
            except (AssertionError, RuntimeError):
                all_pred_vp += [float("nan")]

    return all_pred_den, all_pred_vp


def custom_plot(
    parameters: list,
    temp_min: float,
    temp_max: float,
    pressure: float,
    checkboxes: list[bool],
) -> Union[list[tuple[str, int, str, str]], list]:
    """
    Custom plot function for PC-SAFT parameters.

    args
    ---------
    parameters: list
      list with PC-SAFT parameters
    temp_min: float
      minimum temperature in Kelvin
    temp_max: float
      maximum temperature in Kelvin
    pressure: float
      pressure in Pa
    checkboxes: list
      list with checks to plot in the order [density, vapor pressure, enthalpy, entropy]
    """
    temp_range = np.linspace(temp_min, temp_max, 100, dtype=np.float64)
    p_range = np.asarray([pressure] * 100, dtype=np.float64)
    states = np.stack((temp_range, p_range), 1)
    prop_fns = [
        pure_den_feos,
        pure_vp_feos,
        pure_h_lv_feos,
        pure_s_lv_feos,
    ]  # more prop later
    prop_ids = ["den_plot", "vp_plot", "h_lv_plot", "s_lv_plot"]
    prop_names = [
        "Liquid Density (mol / m³)",
        "Vapor pressure (Pa)",
        "Enthalpy of vaporization (kJ/mol)",
        "Entropy of vaporization (J/mol/K)",
    ]
    xlegendpos = [0, 0, 0, 0]
    all_plots = []
    plot_sf = checkboxes.pop()

    for prop_fn, prop_id, prop_name, xpos, checkbox in zip(
        prop_fns,
        prop_ids,
        prop_names,
        xlegendpos,
        checkboxes,
    ):
        if checkbox:
            plot_data = {"T": [], "GNN": [], "TML": []}
            for state in states:
                try:

                    prop_for_state = prop_fn(parameters.copy(), state)
                    plot_data["T"].append(state[0])
                    plot_data["GNN"].append(prop_for_state)
                except (AssertionError, RuntimeError) as e:
                    logger.debug(e)
            all_plots.append((json.dumps(plot_data), xpos, prop_name, prop_id))
    if plot_sf:
        surface_tension, temp_st = pure_surface_tension_feos(parameters, [temp_min])
        plot_data = {
            "T": temp_st.tolist(),
            "GNN": surface_tension.tolist(),
            "TML": [],
        }
        all_plots.append(
            (json.dumps(plot_data), 0, "Surface Tension (mN/m)", "st_plot")
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


def get_main_plots_data(inchi: str, parameters: List) -> Tuple[str, str, str]:
    "get main plot data"

    rho_data = pl.read_parquet(osp.join(dataset_dir, "rho_pure.parquet")).filter(
        pl.col("inchi1") == inchi, pl.col("P_kPa").is_close(101.325)
    )
    vp_data = pl.read_parquet(osp.join(dataset_dir, "vp_pure.parquet")).filter(
        pl.col("inchi1") == inchi
    )
    if vp_data.height > 0:
        vp_data = (
            vp_data.with_columns((pl.col("VP_kPa") * 1000).alias("P_Pa"))
            .select("T_K", "P_Pa")
            .sort("T_K")
            .to_numpy()
        )
        _, plotvp = plotdata(parameters, np.array([]), vp_data)
        plotvp = json.dumps(plotvp)
    else:
        plotvp = ""

    if rho_data.height > 0:
        den_data = (
            rho_data.with_columns(
                (pl.col("P_kPa") * 1000).alias("P_Pa"),
                (pl.col("rho") * 1000 / pl.col("molweight1")).alias("den"),
            )
            .select("T_K", "P_Pa", "den")
            .sort("T_K")
            .to_numpy()
        )
        plotden, _ = plotdata(parameters, den_data, np.array([]))
        plotden = json.dumps(plotden)
    else:
        plotden = ""

    molimg = plotmol(inchi)
    return plotden, plotvp, molimg


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


def get_custom_plots_data(
    pred: list,
    plot_config: CustomPlotConfigForm,
    checkboxes: tuple[
        RhoCheckForm,
        VPCheckForm,
        HlvCheckForm,
        SlvCheckForm,
        PhaseDiagramCheckForm,
        STCheckForm,
    ],
) -> tuple[list, list]:
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
    try:
        custom_plots = custom_plot(
            pred,
            plot_config.cleaned_data["temp_min"],
            plot_config.cleaned_data["temp_max"],
            plot_config.cleaned_data["pressure"],
            [
                rho_checkbox_.cleaned_data["rho_checkbox"],
                vp_checkbox_.cleaned_data["vp_checkbox"],
                h_lv_checkbox_.cleaned_data["h_lv_checkbox"],
                s_lv_checkbox_.cleaned_data["s_lv_checkbox"],
                st_checkbox_.cleaned_data["st_checkbox"],
            ],
        )
    except RuntimeError as err:
        logger.debug(err)
        custom_plots = []
    phase_diagrams = []
    if phase_diagram_checkbox_.cleaned_data["phase_diagram_checkbox"]:
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
        except RuntimeError as err:
            logger.debug(err)
            phase_diagrams = []
    return phase_diagrams, custom_plots


def _get_ternary_lle_data(
    params: List[List[float]], kij_matrix: List[List[float]], state: List[float]
) -> Dict[str, List[float]]:
    t, p = state  # Temperatura (K) e pressão (Pa)

    def _grid(n_pts: int = 25):
        xi = np.linspace(1e-5, 0.999, n_pts, dtype=np.float64)
        x1_m, x2_m = np.meshgrid(xi, xi, indexing="xy")
        x3_m = 1.0 - x1_m - x2_m
        return x1_m, x2_m, x3_m, (x3_m >= 0.0)

    def _collect_tie_lines(x1_m, x2_m, x3_m, mask):
        valid_idx = np.argwhere(mask)
        ternary_data = {"x0": [], "x1": [], "x2": [], "y0": [], "y1": [], "y2": []}
        for i, j in valid_idx:
            try:
                lle = mix_lle_feos(
                    params,
                    [t, p, x1_m[i, j].item(), x2_m[i, j].item(), x3_m[i, j].item()],
                    kij_matrix,
                )
            except (RuntimeError, ValueError):
                continue
            # For LLE, y is one phase and x is the other phase
            ternary_data["x0"].extend(lle["x0"])
            ternary_data["x1"].extend(lle["x1"])
            ternary_data["x2"].extend(lle["x2"])
            ternary_data["y0"].extend(lle["y0"])
            ternary_data["y1"].extend(lle["y1"])
            ternary_data["y2"].extend(lle["y2"])
        return ternary_data

    x1, x2, x3, mask = _grid()
    return _collect_tie_lines(x1, x2, x3, mask)


def get_mixture_plots_data(
    para_pred_list: List[List],
    mole_fractions_list: List[float],
    config: Tuple[
        CustomPlotConfigForm,
        TernaryLLECheckForm,
        BinaryVLECheckForm,
        BinaryLLECheckForm,
    ],
    kij_matrix: List[List[float]],
) -> Tuple[Tuple[List[Tuple[str, int, str]], List[str]], str, str, str]:
    "get mixture plots data"

    plot_config, ternary_lle_checkform, binary_vle_checkform, binary_lle_checkform = (
        config
    )

    plot_config.full_clean()
    ternary_lle_checkform.full_clean()
    binary_vle_checkform.full_clean()
    binary_lle_checkform.full_clean()
    mixture_plot = mixture_plots(
        para_pred_list,
        (
            mole_fractions_list,
            plot_config.cleaned_data["temp_min"],
            plot_config.cleaned_data["temp_max"],
            plot_config.cleaned_data["pressure"],
        ),
        kij_matrix,
    )

    try:
        if ternary_lle_checkform.cleaned_data["ternary_lle_checkbox"] is False:
            raise ValueError("Ternary LLE checkbox not selected.")
        if len(para_pred_list) != 3:
            raise ValueError("LLE phase diagram only for ternary mixtures.")

        ternary_lle_phase_diagram_data = _get_ternary_lle_data(
            para_pred_list,
            kij_matrix,
            [
                plot_config.cleaned_data["temp_min"],
                plot_config.cleaned_data["pressure"],
            ],
        )
        ternary_lle_phase_diagram_data = json.dumps(ternary_lle_phase_diagram_data)
    except (ValueError, RuntimeError) as err:
        logger.debug(err)
        ternary_lle_phase_diagram_data = ""

    try:
        if binary_lle_checkform.cleaned_data["binary_lle_checkbox"] is False:
            raise ValueError("Binary LLE checkbox not selected.")
        if len(para_pred_list) != 2:
            raise ValueError("LLE phase diagram only for binary mixtures.")

        try:
            binary_lle_phase_diagram_data = mix_lle_diagram_feos(
                para_pred_list,
                [
                    plot_config.cleaned_data["temp_min"],
                    plot_config.cleaned_data["pressure"],
                    *mole_fractions_list,
                ],
                kij_matrix,
            )
        except ValueError:
            binary_lle_phase_diagram_data = {"x0": [], "y0": [], "temperature": []}
        binary_lle_phase_diagram_data = json.dumps(binary_lle_phase_diagram_data)
    except (ValueError, RuntimeError) as err:
        logger.debug(err)
        binary_lle_phase_diagram_data = ""

    try:
        if binary_vle_checkform.cleaned_data["binary_vle_checkbox"] is False:
            raise ValueError("Binary VLE checkbox not selected.")
        if len(para_pred_list) != 2:
            raise ValueError("VLE phase diagram only for binary mixtures.")

        try:
            vle_phase_diagram_data = mix_vle_diagram_feos(
                para_pred_list,
                [
                    plot_config.cleaned_data["pressure"],
                ],
                kij_matrix,
            )
        except ValueError:
            vle_phase_diagram_data = {"x0": [], "y0": [], "temperature": []}
        vle_phase_diagram_data = json.dumps(vle_phase_diagram_data)
    except (ValueError, RuntimeError) as err:
        logger.debug(err)
        vle_phase_diagram_data = ""

    return (
        mixture_plot,
        binary_lle_phase_diagram_data,
        vle_phase_diagram_data,
        ternary_lle_phase_diagram_data,
    )


def mixture_plots(
    para_pred_list: List[List[float]],
    state_list: Tuple[List[float], float, float, float],
    kij_matrix: List[List[float]],
) -> Tuple[List[Tuple[str, int, str]], List[str]]:
    "get mixture plots data"
    mole_fractions_list, temp_min, temp_max, pressure = state_list
    temp_range = np.linspace(temp_min, temp_max, 100, dtype=np.float64)
    p_range = np.asarray([pressure] * 100, dtype=np.float64)
    mole_fractions = np.asarray([mole_fractions_list] * 100, dtype=np.float64)
    states = np.stack((temp_range, p_range), 1)
    states = np.hstack((states, mole_fractions))
    all_plots = []
    vp_plots = []

    plot_data_den = {"T": [], "GNN": [], "TML": []}
    plot_data_bubble = {"T": [], "GNN": [], "TML": []}
    plot_data_dew = {"T": [], "GNN": [], "TML": []}
    for state in states:
        try:

            prop_for_state = mix_den_feos(
                para_pred_list.copy(),
                state,
                kij_matrix.copy(),
            )
            bubble_for_state, dew_for_state = mix_vp_feos(
                para_pred_list.copy(),
                state,
                kij_matrix.copy(),
            )
            plot_data_den["T"].append(state[0])
            plot_data_den["GNN"].append(prop_for_state)
            plot_data_bubble["T"].append(state[0])
            plot_data_bubble["GNN"].append(bubble_for_state)
            plot_data_dew["T"].append(state[0])
            plot_data_dew["GNN"].append(dew_for_state)
        except (AssertionError, RuntimeError) as e:
            logger.debug(e)
    all_plots.append(
        (json.dumps(plot_data_den), 0, "Liquid Density (mol / m³)", "den_plot")
    )
    vp_plots.append(json.dumps(plot_data_bubble))
    vp_plots.append(json.dumps(plot_data_dew))

    return all_plots, vp_plots


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
    plotden, plotvp, molimg = get_main_plots_data(inchi, pred)
    output = True

    plot_checkbox.full_clean()
    phase_diagrams, custom_plots = [], []
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
    return {
        "smiles": smiles,
        "inchi": inchi,
        "pred": pred,
        "plotden": plotden,
        "plotvp": plotvp,
        "molimg": molimg,
        "output": output,
        "custom_plots": custom_plots,
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
        "plotden": False,
        "plotvp": False,
        "den_data": "",
        "vp_data": "",
        "mol_data": "",
        "custom_plots": [],
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
                "plotden": post_data["plotden"] != "",
                "plotvp": post_data["plotvp"] != "",
                "den_data": post_data["plotden"],
                "vp_data": post_data["plotvp"],
                "mol_data": post_data["molimg"],
                "custom_plots": post_data["custom_plots"],
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
    para_pred_for_plot = []
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
            para_pred_for_plot.append(para_pred)
        mixture_plots_ = get_mixture_plots_data(
            para_pred_for_plot,
            mole_fractions_list,
            (
                plot_config,
                ternary_lle_checkform,
                binary_vle_checkform,
                binary_lle_checkform,
            ),
            kij_matrix,
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
            "mixture_plots": post_data["mixture_plots"][0][0],
            "vp_plots": post_data["mixture_plots"][0][1],
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
