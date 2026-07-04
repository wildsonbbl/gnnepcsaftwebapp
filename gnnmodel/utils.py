"Module for utils like plotting data."

import json
import os.path as osp
from typing import List, Optional, Tuple, Union

import onnxruntime as ort
from gnnepcsaft.pcsaft.pcsaft_feos import critical_points_feos
from gnnepcsaft_mcp_server.utils import predict_pcsaft_parameters
from rdkit.Chem import AllChem as Chem

from . import logger
from .forms import MixtureForms, PureForms
from .utils_data import (
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
    TernaryVleTxParams,
    mix_den,
    mix_lle,
    mix_ternary_lle,
    mix_ternary_vle_tx_fixed,
    mix_vle,
    mix_vle_pxy,
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


def _get_vp_data(smiles: str, temp_min: float, temp_max: float, npoints: int):
    plot_data = {}
    plot_data["legends"] = [
        "GNN",
        "GNN",
        "ThermoML Archive**",
    ]
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
    return plot_data


def _get_h_lv_data(smiles: str, temp_min: float, temp_max: float, npoints: int):
    plot_data = {}
    plot_data["legends"] = [
        "GNN",
        "GNN",
        "ThermoML Archive**",
    ]
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
    return plot_data


def _get_s_lv_data(smiles: str, temp_min: float, temp_max: float, npoints: int):
    plot_data = {}
    plot_data["legends"] = [
        "GNN",
        "GNN",
        "ThermoML Archive**",
    ]
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
    return plot_data


def _get_st_data(smiles: str, temp_min: float):
    plot_data = {}
    plot_data["legends"] = [
        "GNN",
        "GNN",
        "ThermoML Archive**",
    ]
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
    return plot_data


def _get_rho_data(
    smiles: str, temp_min: float, temp_max: float, pressure: float, npoints: int
):
    plot_data = {}
    plot_data["legends"] = [
        "GNN",
        "GNN",
        "ThermoML Archive**",
    ]
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
    return plot_data


def _get_mix_vle_b_data(
    smiles_list: List[str], kij_matrix: List[List[float]], npoints: int, temp_min: float
):
    plot_data = {}
    plot_data["legends"] = [
        "GNN Bubble P.",
        "GNN Dew P.",
        "Exp. Bubble P. (ThermoML Archive**)",
    ]
    try:
        plot_data["GNN"] = mix_vle_pxy(
            smiles_list=smiles_list,
            kij_matrix=kij_matrix,
            temperature=temp_min,
            npoints=npoints,
        )
    except (AssertionError, RuntimeError, ValueError) as e:
        logger.debug(e)
    try:
        plot_data["TML"] = ([], [])
        exp_data = retrieve_vle_pxy_binary_data(
            smiles_list=smiles_list,
            temperature=temp_min,
        )
        if exp_data is not None:
            exp_data[:, 1] *= 1000
            plot_data["TML"] = exp_data.T.tolist()
    except (AssertionError, RuntimeError, ValueError) as e:
        logger.debug(e)
    return plot_data


def _get_mix_vle_t_data(
    smiles_list: List[str],
    kij_matrix: List[List[float]],
    npoints: int,
    mole_fractions_list: List[float],
    temp_min: float,
):
    plot_data = {}
    plot_data["legends"] = [
        "GNN Bubble P.",
        "GNN Dew P.",
        "Exp. Bubble P. (ThermoML Archive**)",
    ]
    solvent_ratio = mole_fractions_list[1] / (
        mole_fractions_list[1] + mole_fractions_list[2]
    )
    try:
        plot_data["GNN"] = mix_ternary_vle_tx_fixed(
            TernaryVleTxParams(
                smiles_list=smiles_list,
                kij_matrix=kij_matrix,
                temperature=temp_min,
                solvent_ratio=solvent_ratio,
                npoints=npoints,
            )
        )
    except (AssertionError, RuntimeError, ValueError) as e:
        logger.debug(e)
    try:
        plot_data["TML"] = ([], [])
        exp_data = retrieve_vle_ternary_tx_fixed_data(
            smiles_list=smiles_list,
            temperature=temp_min,
            solvent_ratio=solvent_ratio,
        )
        if exp_data is not None:
            exp_data[:, 1] *= 1000
            plot_data["TML"] = exp_data.T.tolist()
    except (AssertionError, RuntimeError, ValueError) as e:
        logger.debug(e)
    return plot_data, solvent_ratio


def _get_mix_vp_data(
    smiles_list: List[str],
    kij_matrix: List[List[float]],
    npoints: int,
    mole_fractions_list: List[float],
    temp_min: float,
    temp_max: float,
):
    plot_data = {}
    plot_data["legends"] = [
        "GNN Bubble P.",
        "GNN Dew P.",
        "Exp. Bubble P. (ThermoML Archive**)",
    ]
    try:
        plot_data["GNN"] = mix_vp(
            MixVpParams(
                smiles_list=smiles_list,
                mole_fractions=mole_fractions_list,
                kij_matrix=kij_matrix,
                min_temp=temp_min,
                max_temp=temp_max,
                npoints=npoints,
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
    except (AssertionError, RuntimeError, ValueError) as e:
        logger.debug(e)
    return plot_data


def _get_mix_rho_data(
    smiles_list: List[str],
    kij_matrix: List[List[float]],
    npoints: int,
    mole_fractions_list: List[float],
    temp_min: float,
    temp_max: float,
    pressure: float,
):
    plot_data = {}
    plot_data["legends"] = [
        "GNN",
        "GNN",
        "ThermoML Archive**",
    ]
    try:
        plot_data["GNN"] = mix_den(
            MixDenParams(
                smiles_list=smiles_list,
                mole_fractions=mole_fractions_list,
                kij_matrix=kij_matrix,
                min_temp=temp_min,
                max_temp=temp_max,
                pressure=pressure,
                npoints=npoints,
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
    return plot_data


def _get_binary_vle_data(
    smiles_list: List[str], forms: MixtureForms, kij_matrix: List[List[float]]
) -> str:
    vle_phase_diagram_data = ""
    try:
        if forms.binary_vle_checkform.cleaned_data["binary_vle_checkbox"] is False:
            raise ValueError("Binary VLE checkbox not selected.")
        if len(smiles_list) != 2:
            raise ValueError("VLE phase diagram only for binary mixtures.")
        if forms.plot_config.cleaned_data["pressure"] is None:
            raise ValueError("Missing pressure")

        _vle_phase_diagram_data = mix_vle(
            smiles_list=smiles_list,
            kij_matrix=kij_matrix,
            pressure=forms.plot_config.cleaned_data["pressure"],
            npoints=forms.plot_config.cleaned_data["npoints"] or 10,
        )
        _vle_phase_diagram_data["exp_T"] = []
        _vle_phase_diagram_data["exp_x0"] = []
        try:
            exp_data = retrieve_vle_binary_data(
                smiles_list=smiles_list,
                pressure=forms.plot_config.cleaned_data["pressure"] / 1000,
            )
            if exp_data is not None:
                _vle_phase_diagram_data["exp_T"].extend(exp_data[:, 0].tolist())
                _vle_phase_diagram_data["exp_x0"].extend(exp_data[:, 1].tolist())

        except (AssertionError, RuntimeError) as e:
            logger.debug(e)
        vle_phase_diagram_data = json.dumps(_vle_phase_diagram_data)
    except (ValueError, RuntimeError) as err:
        logger.debug(err)
    return vle_phase_diagram_data


def _get_binary_lle_data(
    smiles_list: List[str],
    mole_fractions_list: List[float],
    forms: MixtureForms,
    kij_matrix: List[List[float]],
) -> str:
    binary_lle_phase_diagram_data = ""
    try:
        if forms.binary_lle_checkform.cleaned_data["binary_lle_checkbox"] is False:
            raise ValueError("Binary LLE checkbox not selected.")
        if len(smiles_list) != 2:
            raise ValueError("LLE phase diagram only for binary mixtures.")
        if forms.plot_config.cleaned_data["temp_min"] is None:
            raise ValueError("Missing minimum temperature")
        if forms.plot_config.cleaned_data["pressure"] is None:
            raise ValueError("Missing pressure")

        _binary_lle_phase_diagram_data = mix_lle(
            MixLLEParams(
                smiles_list=smiles_list,
                mole_fractions=mole_fractions_list,
                kij_matrix=kij_matrix,
                temperature=forms.plot_config.cleaned_data["temp_min"],
                pressure=forms.plot_config.cleaned_data["pressure"],
                npoints=forms.plot_config.cleaned_data["npoints"] or 10,
            )
        )
        _binary_lle_phase_diagram_data["exp_T"] = []
        _binary_lle_phase_diagram_data["exp_x0"] = []
        try:
            exp_data = retrieve_lle_binary_data(
                smiles_list=smiles_list,
                pressure=forms.plot_config.cleaned_data["pressure"] / 1000,
            )
            if exp_data is not None:
                _binary_lle_phase_diagram_data["exp_T"].extend(exp_data[:, 0].tolist())
                _binary_lle_phase_diagram_data["exp_x0"].extend(exp_data[:, 1].tolist())

        except (AssertionError, RuntimeError) as e:
            logger.debug(e)
        try:
            exp_data = retrieve_vle_binary_data(
                smiles_list=smiles_list,
                pressure=forms.plot_config.cleaned_data["pressure"] / 1000,
            )
            if exp_data is not None:
                _binary_lle_phase_diagram_data["exp_T"].extend(exp_data[:, 0].tolist())
                _binary_lle_phase_diagram_data["exp_x0"].extend(exp_data[:, 1].tolist())

        except (AssertionError, RuntimeError) as e:
            logger.debug(e)
        binary_lle_phase_diagram_data = json.dumps(_binary_lle_phase_diagram_data)
    except (ValueError, RuntimeError) as err:
        logger.debug(err)
    return binary_lle_phase_diagram_data


def _get_ternary_lle_data(
    smiles_list: List[str], forms: MixtureForms, kij_matrix: List[List[float]]
) -> str:
    ternary_lle_phase_diagram_data = ""
    try:
        if forms.ternary_lle_checkform.cleaned_data["ternary_lle_checkbox"] is False:
            raise ValueError("Ternary LLE checkbox not selected.")
        if len(smiles_list) != 3:
            raise ValueError("LLE/VLE phase diagram only for ternary mixtures.")
        if forms.plot_config.cleaned_data["temp_min"] is None:
            raise ValueError("Missing minimum temperature")
        if forms.plot_config.cleaned_data["pressure"] is None:
            raise ValueError("Missing pressure")

        _ternary_lle_phase_diagram_data = mix_ternary_lle(
            smiles_list=smiles_list,
            kij_matrix=kij_matrix,
            temperature=forms.plot_config.cleaned_data["temp_min"],
            pressure=forms.plot_config.cleaned_data["pressure"],
            npoints=forms.plot_config.cleaned_data["npoints"] or 10,
        )
        _ternary_lle_phase_diagram_data["exp_x0"] = []
        _ternary_lle_phase_diagram_data["exp_x1"] = []
        _ternary_lle_phase_diagram_data["exp_x2"] = []
        try:
            exp_data = retrieve_lle_ternary_data(
                smiles_list=smiles_list,
                pressure=forms.plot_config.cleaned_data["pressure"] / 1000,
                temperature=forms.plot_config.cleaned_data["temp_min"],
            )
            if exp_data is not None:
                _ternary_lle_phase_diagram_data["exp_x0"].extend(
                    exp_data[:, 0].tolist()
                )
                _ternary_lle_phase_diagram_data["exp_x1"].extend(
                    exp_data[:, 1].tolist()
                )
                _ternary_lle_phase_diagram_data["exp_x2"].extend(
                    (1 - exp_data[:, 0] - exp_data[:, 1]).tolist()
                )

        except (AssertionError, RuntimeError) as e:
            logger.debug(e)
        try:
            exp_data = retrieve_vle_ternary_data(
                smiles_list=smiles_list,
                pressure=forms.plot_config.cleaned_data["pressure"] / 1000,
                temperature=forms.plot_config.cleaned_data["temp_min"],
            )
            if exp_data is not None:
                _ternary_lle_phase_diagram_data["exp_x0"].extend(
                    exp_data[:, 0].tolist()
                )
                _ternary_lle_phase_diagram_data["exp_x1"].extend(
                    exp_data[:, 1].tolist()
                )
                _ternary_lle_phase_diagram_data["exp_x2"].extend(
                    (1 - exp_data[:, 0] - exp_data[:, 1]).tolist()
                )

        except (AssertionError, RuntimeError) as e:
            logger.debug(e)
        ternary_lle_phase_diagram_data = json.dumps(_ternary_lle_phase_diagram_data)
    except (ValueError, RuntimeError) as err:
        logger.debug(err)
    return ternary_lle_phase_diagram_data


def _full_clean_forms(forms: MixtureForms):
    forms.plot_config.full_clean()
    forms.rho_checkform.full_clean()
    forms.vp_checkform.full_clean()
    forms.ternary_lle_checkform.full_clean()
    forms.binary_vle_checkform.full_clean()
    forms.binary_lle_checkform.full_clean()
    forms.binary_vlepxy_checkform.full_clean()
    forms.ternary_vlepx_checkform.full_clean()


def _get_mix_checkboxes(forms: MixtureForms) -> List[str]:
    checkboxes = []
    if forms.ternary_vlepx_checkform.cleaned_data["ternary_vlepx_checkbox"]:
        checkboxes.append("tvlepx")
    if forms.binary_vlepxy_checkform.cleaned_data["binary_vlepxy_checkbox"]:
        checkboxes.append("bvlepxy")
    if forms.rho_checkform.cleaned_data["rho_checkbox"]:
        checkboxes.append("rho")
    if forms.vp_checkform.cleaned_data["vp_checkbox"]:
        checkboxes.append("vp")
    return checkboxes


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
    selected_checkboxes: list
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
        plot_data = _get_rho_data(smiles, temp_min, temp_max, pressure, npoints)
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
        plot_data = _get_vp_data(smiles, temp_min, temp_max, npoints)
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
        plot_data = _get_h_lv_data(smiles, temp_min, temp_max, npoints)
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
        plot_data = _get_s_lv_data(smiles, temp_min, temp_max, npoints)
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
        plot_data = _get_st_data(smiles, temp_min)
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


def get_pure_plots_data(
    smiles: str, forms: PureForms
) -> Tuple[List[List[float]], List]:
    "get custom plots data for pure component"

    forms.plot_config.full_clean()
    forms.rho_checkbox.full_clean()
    forms.vp_checkbox.full_clean()
    forms.h_lv_checkbox.full_clean()
    forms.s_lv_checkbox.full_clean()
    forms.phase_diagram_checkbox.full_clean()
    forms.st_checkbox.full_clean()
    custom_plots = []
    phase_diagrams = []

    selected_checkboxes = []
    if forms.rho_checkbox.cleaned_data["rho_checkbox"]:
        selected_checkboxes.append("den_plot")
    if forms.vp_checkbox.cleaned_data["vp_checkbox"]:
        selected_checkboxes.append("vp_plot")
    if forms.h_lv_checkbox.cleaned_data["h_lv_checkbox"]:
        selected_checkboxes.append("h_lv_plot")
    if forms.s_lv_checkbox.cleaned_data["s_lv_checkbox"]:
        selected_checkboxes.append("s_lv_plot")
    if forms.st_checkbox.cleaned_data["st_checkbox"]:
        selected_checkboxes.append("st_plot")

    try:
        custom_plots = pure_plots(
            smiles=smiles,
            temp_min=forms.plot_config.cleaned_data["temp_min"],
            temp_max=forms.plot_config.cleaned_data["temp_max"],
            pressure=forms.plot_config.cleaned_data["pressure"],
            selected_checkboxes=selected_checkboxes,
            npoints=forms.plot_config.cleaned_data["npoints"] or 10,
        )
    except RuntimeError as err:
        logger.debug(err)

    if forms.phase_diagram_checkbox.cleaned_data["phase_diagram_checkbox"]:
        try:
            if forms.plot_config.cleaned_data["temp_min"] is not None:
                phase_diagrams = pure_phase_diagram(
                    smiles=smiles, min_temp=forms.plot_config.cleaned_data["temp_min"]
                )
        except RuntimeError as err:
            logger.debug(err)
    return phase_diagrams, custom_plots


def get_mixture_plots_data(
    smiles_list: List[str],
    mole_fractions_list: List[float],
    forms: MixtureForms,
    kij_matrix: List[List[float]],
) -> Tuple[List[Tuple[str, str, str, str, str]], str, str, str]:
    "get mixture plots data"

    _full_clean_forms(forms)
    checkboxes = _get_mix_checkboxes(forms)

    mixture_plot = mixture_plots(
        smiles_list=smiles_list,
        state_list=(
            mole_fractions_list,
            forms.plot_config.cleaned_data["temp_min"],
            forms.plot_config.cleaned_data["temp_max"],
            forms.plot_config.cleaned_data["pressure"],
        ),
        kij_matrix=kij_matrix,
        checkboxes=checkboxes,
        npoints=forms.plot_config.cleaned_data["npoints"] or 10,
    )

    ternary_lle_phase_diagram_data = _get_ternary_lle_data(
        smiles_list, forms, kij_matrix
    )

    binary_lle_phase_diagram_data = _get_binary_lle_data(
        smiles_list, mole_fractions_list, forms, kij_matrix
    )

    vle_phase_diagram_data = _get_binary_vle_data(smiles_list, forms, kij_matrix)

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
    checkboxes: List[str],
    npoints: int = 20,
) -> List[Tuple[str, str, str, str, str]]:
    "get mixture plots data"
    mole_fractions_list, temp_min, temp_max, pressure = state_list
    all_plots = []

    if (
        temp_min is not None
        and temp_max is not None
        and pressure is not None
        and "rho" in checkboxes
    ):
        plot_data = _get_mix_rho_data(
            smiles_list,
            kij_matrix,
            npoints,
            mole_fractions_list,
            temp_min,
            temp_max,
            pressure,
        )

        all_plots.append(
            (
                json.dumps(plot_data),
                "Temperature (K)",
                "Liquid Density (mol / m³)",
                f"Density at {pressure} Pa and mole fractions={mole_fractions_list}",
                "mix_den_plot",
            )
        )

    if temp_min is not None and temp_max is not None and "vp" in checkboxes:
        plot_data = _get_mix_vp_data(
            smiles_list, kij_matrix, npoints, mole_fractions_list, temp_min, temp_max
        )

        all_plots.append(
            (
                json.dumps(plot_data),
                "Temperature (K)",
                "Pressure (Pa)",
                f"VLE at mole fractions={mole_fractions_list}",
                "mix_vp_plot",
            )
        )

    if "tvlepx" in checkboxes and len(smiles_list) == 3:
        plot_data, solvent_ratio = _get_mix_vle_t_data(
            smiles_list, kij_matrix, npoints, mole_fractions_list, temp_min
        )

        all_plots.append(
            (
                json.dumps(plot_data),
                "x1",
                "Pressure (Pa)",
                f"VLE P-x1 at T={temp_min} K and x2/(x2+x3)={solvent_ratio:.3f}",
                "mix_tvlepx_plot",
            )
        )

    if "bvlepxy" in checkboxes and len(smiles_list) == 2 and temp_min is not None:
        plot_data = _get_mix_vle_b_data(smiles_list, kij_matrix, npoints, temp_min)

        all_plots.append(
            (
                json.dumps(plot_data),
                "x1",
                "Pressure (Pa)",
                f"VLE P-x-y at T={temp_min} K",
                "mix_bvlepxy_plot",
            )
        )

    return all_plots
