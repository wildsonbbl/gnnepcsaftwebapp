"Pure screen utilities"

from typing import List, Tuple

import numpy as np
from gnnepcsaft.pcsaft.pcsaft_feos import (
    phase_diagram_feos,
    pure_den_feos,
    pure_h_lv_feos,
    pure_s_lv_feos,
    pure_surface_tension_feos,
    pure_vp_feos,
)
from gnnepcsaft_mcp_server.utils import predict_pcsaft_parameters


def pure_den(
    smiles: str, min_temp: float, max_temp: float, pressure: float, npoints: int
) -> Tuple[List[float], List[float]]:
    "Calculate pure-component density using PC-SAFT EOS"
    parameters = predict_pcsaft_parameters(smiles)
    temperatures = np.linspace(min_temp, max_temp, num=npoints).tolist()

    densities = [pure_den_feos(parameters, [T, pressure]) for T in temperatures]
    return temperatures, densities


def pure_vp(
    smiles: str, min_temp: float, max_temp: float, npoints: int
) -> Tuple[List[float], List[float]]:
    "Calculate pure-component vapor pressure using PC-SAFT EOS"
    parameters = predict_pcsaft_parameters(smiles)
    temperatures = np.linspace(min_temp, max_temp, num=npoints).tolist()

    vapor_pressures = [pure_vp_feos(parameters, [T]) for T in temperatures]
    return temperatures, vapor_pressures


def pure_h_lv(
    smiles: str, min_temp: float, max_temp: float, npoints: int
) -> Tuple[List[float], List[float]]:
    "Calculate pure-component enthalpy of vaporization using PC-SAFT EOS"
    parameters = predict_pcsaft_parameters(smiles)
    temperatures = np.linspace(min_temp, max_temp, num=npoints).tolist()

    h_lvs = [pure_h_lv_feos(parameters, [T]) for T in temperatures]
    return temperatures, h_lvs


def pure_s_lv(
    smiles: str, min_temp: float, max_temp: float, npoints: int
) -> Tuple[List[float], List[float]]:
    "Calculate pure-component entropy of vaporization using PC-SAFT EOS"
    parameters = predict_pcsaft_parameters(smiles)
    temperatures = np.linspace(min_temp, max_temp, num=npoints).tolist()

    s_lvs = [pure_s_lv_feos(parameters, [T]) for T in temperatures]
    return temperatures, s_lvs


def pure_surface_tension(
    smiles: str,
    min_temp: float,
) -> Tuple[List[float], List[float]]:
    "Calculate pure-component surface tension (mN/m) using PC-SAFT EOS"
    parameters = predict_pcsaft_parameters(smiles)

    surface_tensions, temperatures = pure_surface_tension_feos(parameters, [min_temp])

    return temperatures.tolist(), surface_tensions.tolist()


def pure_phase_diagram(
    smiles: str,
    min_temp: float,
) -> List[List[float]]:
    "Calculate pure-component phase diagram using PC-SAFT EOS"
    parameters = predict_pcsaft_parameters(smiles)

    output = phase_diagram_feos(parameters, [min_temp])

    return [
        output["temperature"],
        output.get("pressure", output.get("pressure vapor", [0])),
        output["density liquid"],
        output["density vapor"],
    ]
