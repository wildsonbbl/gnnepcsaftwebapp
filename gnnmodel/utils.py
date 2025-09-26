"Module for utils like plotting data."

import csv
import json
import os.path as osp
from json import loads
from typing import Any, Dict, List, Literal, Optional, Tuple, Union
from urllib.parse import quote
from urllib.request import HTTPError, urlopen

import numpy as np
import onnxruntime as ort
from django.conf import settings
from gnnepcsaft.data.rdkit_util import inchitosmiles, mw
from gnnepcsaft.epcsaft.epcsaft_feos import (
    critical_points_feos,
    mix_den_feos,
    mix_lle_diagram_feos,
    mix_vle_diagram_feos,
    mix_vp_feos,
    phase_diagram_feos,
    pure_den_feos,
    pure_h_lv_feos,
    pure_s_lv_feos,
    pure_surface_tension_feos,
    pure_vp_feos,
)
from rdkit.Chem import AllChem as Chem

from gnnepcsaft_mcp_server.utils import predict_epcsaft_parameters

from . import logger
from .forms import (
    CustomPlotCheckForm,
    CustomPlotConfigForm,
    GoogleAPIKeyForm,
    HlvCheckForm,
    InChIorSMILESareaInputforMixture,
    InChIorSMILESinput,
    PhaseDiagramCheckForm,
    RhoCheckForm,
    SlvCheckForm,
    STCheckForm,
    VPCheckForm,
)
from .models import GnnepcsaftPara, ThermoMLDenData, ThermoMLVPData

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
    "Critical temperature (K)",
    "Critical pressure (Bar)",
]

file_dir = osp.dirname(__file__)
dataset_dir = osp.join(file_dir, "data")


def make_dataset() -> dict[str, tuple[List[List[float]], List[List[float]]]]:
    "Make dict dataset for inference."
    import polars as pl  # # pylint: disable = C0415

    data = pl.read_parquet(osp.join(dataset_dir, "thermoml/raw/pure.parquet"))
    inchis = data.unique("inchi1")["inchi1"].to_list()
    _tml_data = {}
    for inchi in inchis:
        try:
            molw = mw(inchi)
        except (TypeError, ValueError) as e:
            logger.debug("Error for InChI: %s \n\n %s", inchi, e)
            continue
        vp = (
            data.filter(pl.col("inchi1") == inchi, pl.col("tp") == 3)
            .select("TK", "PPa", "phase", "tp", "m")
            .sort("TK")
            .to_numpy()
            .tolist()
        )
        rho = np.copy(
            data.filter(pl.col("inchi1") == inchi, pl.col("tp") == 1)
            .select("TK", "PPa", "phase", "tp", "m")
            .filter(pl.col("PPa") == 101325)
            .sort("TK")
            .to_numpy()
        )

        rho[:, -1] *= 1000 / molw  # convert to mol/ m³

        _tml_data[inchi] = (rho.tolist(), vp)
    return _tml_data


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
            "TML": (vp[:, -1] / 1000).tolist(),  # Pressure in kPa
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


def para_update_database(app, schema_editor):  # pylint: disable=W0613
    "fn to update database with pcsaft parameters"
    from tqdm import tqdm  # pylint: disable=C0415

    tml_data = make_dataset()

    data = []
    for inchi in tqdm(tml_data):
        try:
            smiles = inchitosmiles(inchi, False, False)
            para = predict_epcsaft_parameters(smiles)
        except ValueError as e:
            logger.debug("%s: \n\n%s", e, inchi)
            continue
        if para[0] is None:
            continue
        new_comp = GnnepcsaftPara(
            inchi=inchi,
            smiles=smiles,
            m=para[0],
            sigma=para[1],
            e=para[2],
            k_ab=para[3],
            e_ab=para[4],
            mu=para[5],
            na=para[6],
            nb=para[7],
        )
        data.append(
            [
                para[0],
                para[1],
                para[2],
                para[3],
                para[4],
                para[5],
                para[6],
                para[7],
                inchi,
                smiles,
            ]
        )
        new_comp.save()
    with open(
        settings.BASE_DIR / "gnnmodel/static/mydata.csv", "w", encoding="UTF-8"
    ) as f:
        writer = csv.writer(f, delimiter="|", lineterminator="\n")

        writer.writerow(
            ["m", "sigma", "e", "k_ab", "e_ab", "mu", "na", "nb", "inchi", "smiles"]
        )
        writer.writerows(data)
    logger.info("Updated database with pcsaft parameters")


def thermo_update_database(app, schema_editor):  # pylint: disable=W0613
    "fn to update database with plotden, plotvp"
    from tqdm import tqdm  # pylint: disable=C0415

    tml_data = make_dataset()

    for inchi in tqdm(tml_data):
        plotden, plotvp = tml_data[inchi]
        if plotden:
            new_comp = ThermoMLDenData(inchi=inchi, den=json.dumps(plotden))
            new_comp.save()
        if plotvp:
            new_comp = ThermoMLVPData(inchi=inchi, vp=json.dumps(plotvp))
            new_comp.save()
    logger.info("Updated database with plotden, plotvp")


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
                all_pred_vp += [pure_vp_feos(parameters, state) / 1000]  # to kPA
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
        surface_tension, temp_st = pure_surface_tension_feos(
            parameters + ["unk", "unk", 1], [temp_min]
        )
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
    pred = predict_epcsaft_parameters(smiles)

    try:
        critical_points = critical_points_feos(pred.copy())
    except RuntimeError:
        critical_points = [0.0, 0.0]

    pred.append(critical_points[0])
    pred.append(critical_points[1] * 0.00001)  # convert from Pa to Bar
    return pred


def get_main_plots_data(inchi: str, parameters: List) -> Tuple[str, str, str]:
    "get main plot data"

    alldata = ThermoMLVPData.objects.filter(inchi=inchi).all()  # pylint: disable=E1101
    if len(alldata) > 0:
        vp_data = np.asarray(json.loads(alldata[0].vp))
        _, plotvp = plotdata(parameters, np.array([]), vp_data)
        plotvp = json.dumps(plotvp)
    else:
        plotvp = ""

    alldata = ThermoMLDenData.objects.filter(inchi=inchi).all()  # pylint: disable=E1101
    if len(alldata) > 0:
        den_data = np.asarray(json.loads(alldata[0].den))
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
        GoogleAPIKeyForm(request.POST),
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


def get_mixture_plots_data(
    para_pred_list: List[List],
    mole_fractions_list: List[float],
    plot_config: CustomPlotConfigForm,
    kij_matrix: List[List[float]],
) -> Tuple[Tuple[List[Tuple[str, int, str]], List[str]], str, str]:
    "get mixture plots data"

    plot_config.full_clean()
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

    lle_phase_diagram_data = mix_lle_diagram_feos(
        para_pred_list,
        [
            plot_config.cleaned_data["temp_min"],
            plot_config.cleaned_data["pressure"],
            *mole_fractions_list,
        ],
        kij_matrix,
    )

    vle_phase_diagram_data = mix_vle_diagram_feos(
        para_pred_list,
        [
            plot_config.cleaned_data["pressure"],
        ],
        kij_matrix,
    )

    return (
        mixture_plot,
        json.dumps(lle_phase_diagram_data),
        json.dumps(vle_phase_diagram_data),
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


def pure_phase(
    vapor_pressure: float, system_pressure: float
) -> Literal["liquid", "vapor"]:
    """
    Given the vapor pressure and system pressure, return the phase of the molecule.
    Both pressures must be in the same unit.

    Args:
        vapor_pressure (float): The calculated vapor pressure of the pure component.
        system_pressure (float): The actual system pressure.

    """
    assert isinstance(vapor_pressure, (int, float)), "vapor_pressure must be a number"
    assert isinstance(system_pressure, (int, float)), "system_pressure must be a number"
    assert vapor_pressure > 0, "vapor_pressure must be positive"
    assert system_pressure > 0, "system_pressure must be positive"

    return "liquid" if vapor_pressure < system_pressure else "vapor"


def mixture_phase(
    bubble_point: float,
    dew_point: float,
    system_pressure: float,
) -> Literal["liquid", "vapor", "two-phase"]:
    """
    Given the bubble/dew point of the mixture and the system pressure,
    return the phase of the mixture.
    All pressures must be in the same unit.

    Args:
        bubble_point (float): The calculated bubble point of the mixture.
        dew_point (float): The calculated dew point of the mixture.
        system_pressure (float): The actual system pressure.
    """
    assert isinstance(bubble_point, (int, float)), "bubble_point must be a number"
    assert isinstance(dew_point, (int, float)), "dew_point must be a number"
    assert isinstance(system_pressure, (int, float)), "system_pressure must be a number"
    assert bubble_point > 0, "bubble_point must be positive"
    assert dew_point > 0, "dew_point must be positive"
    assert system_pressure > 0, "system_pressure must be positive"
    return (
        "liquid"
        if bubble_point < system_pressure
        else ("two-phase" if dew_point <= system_pressure else "vapor")
    )


def pubchem_description(inchi: str) -> str:
    """
    Look for information on PubChem for the InChI.

    Args:
        inchi (str): The InChI of the molecule.
    """
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/description/json?inchi="
        + quote(inchi, safe="")
    )
    try:
        with urlopen(url) as ans:
            ans = loads(ans.read().decode("utf8").strip())
    except (TypeError, HTTPError, ValueError):
        ans = "no data available on this molecule in PubChem."
    return ans


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
        )
    return (
        InChIorSMILESareaInputforMixture(),
        CustomPlotConfigForm(),
    )


def process_mixture_post(
    forms: Tuple[InChIorSMILESareaInputforMixture, CustomPlotConfigForm],
):
    "process the post data from the mixture page"
    form, plot_config = forms
    para_pred_list = []
    para_pred_for_plot = []
    mole_fractions_list = []
    mixture_plots_ = (([], []), "", "")
    output = False

    if form.is_valid():
        inchi_list, smiles_list, mole_fractions_list, kij = form.cleaned_data[
            "text_area"
        ]
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
        for smiles, inchi in zip(smiles_list, inchi_list):
            para_pred = [round(para, 5) for para in get_pred(smiles)]
            para_pred_list.append(para_pred)
            para_pred_for_plot.append(para_pred + [mw(inchi)])
        mixture_plots_ = get_mixture_plots_data(
            para_pred_for_plot,
            mole_fractions_list,
            plot_config,
            kij_matrix,
        )
        output = True

    return {
        "form": form,
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
            "available_params": available_params,
            "parameters_molefractions_list": list(
                zip(post_data["para_pred_list"], post_data["mole_fractions_list"])
            ),
            "mixture_plots": post_data["mixture_plots"][0][0],
            "vp_plots": post_data["mixture_plots"][0][1],
            "lle_phase_diagram_data": post_data["mixture_plots"][1],
            "vle_phase_diagram_data": post_data["mixture_plots"][2],
            "output": post_data["output"],
        }
    return {
        "form": InChIorSMILESareaInputforMixture(),
        "plot_config": CustomPlotConfigForm(),
        "available_params": available_params,
        "parameters_molefractions_list": [],
        "mixture_plots": [],
        "vp_plots": [],
        "lle_phase_diagram_data": "",
        "vle_phase_diagram_data": "",
        "output": False,
    }
