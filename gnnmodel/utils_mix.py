"Mixture screen utilities"

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
from gnnepcsaft.pcsaft.pcsaft_feos import (
    critical_points_feos,
    mix_den_feos,
    mix_lle_diagram_feos,
    mix_lle_feos,
    mix_vle_diagram_feos,
    mix_vp_feos,
    pure_vp_feos,
)
from gnnepcsaft_mcp_server.utils import predict_pcsaft_parameters

from . import logger


@dataclass
class BinaryVlePxyContext:
    """Shared context for binary VLE P-x-y calculations."""

    parameters_list: List[List[float]]
    kij_matrix: List[List[float]]
    temperature: float
    min_pc: float
    max_pc: float


@dataclass
class BinaryVleSeries:
    """Series data for binary VLE P-x-y output."""

    xs: List[float]
    bps: List[float]
    dps: List[float]


@dataclass
class TernaryVleContext:
    """Shared context for ternary VLE P-x calculations."""

    parameters_list: List[List[float]]
    kij_matrix: List[List[float]]
    temperature: float
    min_pc: float
    max_pc: float


@dataclass
class TernaryVleSeries:
    """Series data for ternary VLE P-x output."""

    x1_values: List[float]
    bubble_pressures: List[float]
    dew_pressures: List[float]


@dataclass
class MixDenParams:
    """Inputs for mixture density calculations."""

    smiles_list: List[str]
    mole_fractions: List[float]
    kij_matrix: List[List[float]]
    min_temp: float
    max_temp: float
    pressure: float
    npoints: int


@dataclass
class MixLLEParams:
    """Inputs for mixture LLE calculations."""

    smiles_list: List[str]
    mole_fractions: List[float]
    kij_matrix: List[List[float]]
    temperature: float
    pressure: float
    npoints: int


@dataclass
class TernaryVleTxParams:
    """Inputs for ternary VLE P-x calculations."""

    smiles_list: List[str]
    kij_matrix: List[List[float]]
    temperature: float
    solvent_ratio: float
    npoints: int
    mole_fractions: Optional[List[float]] = None


@dataclass
class MixVpParams:
    """Inputs for mixture vapor pressure calculations."""

    smiles_list: List[str]
    mole_fractions: List[float]
    kij_matrix: List[List[float]]
    min_temp: float
    max_temp: float
    npoints: int


def _binary_critical_points(
    parameters_list: List[List[float]],
) -> Tuple[List[float], List[float]]:
    tcs = []
    pcs = []
    for params in parameters_list:
        tc, pc, _ = critical_points_feos(parameters=params)
        tcs.append(tc)
        pcs.append(pc)
    return tcs, pcs


def _binary_pure_vps(
    parameters_list: List[List[float]],
    temperature: float,
    tcs: List[float],
) -> List[Optional[float]]:
    vps = []
    for params, tc in zip(parameters_list, tcs):
        try:
            if temperature < tc:
                vps.append(pure_vp_feos(parameters=params, state=[temperature]))
            else:
                vps.append(None)
        except RuntimeError:
            vps.append(None)
    return vps


def _validate_binary_vle_pxy(
    temperature: float,
    tcs: List[float],
    vps: List[Optional[float]],
):
    if temperature >= tcs[0] and temperature >= tcs[1]:
        raise ValueError(
            f"Temperature {temperature} K is above the critical temperature "
            f"of both components, which is {tcs[0]:.2f} K and {tcs[1]:.2f} K. "
            "VLE calculation is not meaningful.",
        )
    if vps[0] is None and vps[1] is None:
        raise ValueError(
            f"Unable to compute pure vapor pressures at {temperature} K for either component; "
            "VLE calculation is not meaningful.",
        )


def _binary_more_volatile_is_first(
    vps: List[Optional[float]],
    tcs: List[float],
    temperature: float,
) -> bool:
    vp0, vp1 = vps
    if vp0 is not None and vp1 is not None:
        return vp0 >= vp1
    if vp0 is None and vp1 is not None:
        return True
    if vp1 is None and vp0 is not None:
        return False
    if temperature >= tcs[0] and temperature >= tcs[1]:
        return True
    return True


def _append_binary_pure_point(
    temperature: float,
    tcs: List[float],
    vps: List[Optional[float]],
    series: BinaryVleSeries,
):
    if temperature < tcs[1] and vps[1] is not None:
        series.xs.append(0.0)
        series.bps.append(vps[1])
        series.dps.append(vps[1])


def _binary_vle_pxy_point(
    context: BinaryVlePxyContext,
    x0: float,
):
    try:
        bp, dp = mix_vp_feos(
            parameters=context.parameters_list,
            state=[context.temperature, np.nan, x0, 1 - x0],
            kij_matrix=context.kij_matrix,
        )
    except RuntimeError as exc:
        logger.debug("mix_vle_pxy: Runtime Error at x0=%.4f: %s", x0, exc)
        return None
    except BaseException as exc:  # pylint: disable=W0718
        exception_type = type(exc).__name__
        if exception_type == "PanicException":
            logger.warning("mix_vle_pxy: PanicException at x0=%.4f: %s", x0, exc)
            return None
        logger.exception(
            "mix_vle_pxy: unexpected %s at x0=%.4f",
            exception_type,
            x0,
        )
        raise

    if bp > context.max_pc or dp > context.max_pc:
        logger.warning(
            "mix_vle_pxy: breaking from point above both Pc at x0=%.4f: bp=%.2f, dp=%.2f",
            x0,
            bp,
            dp,
        )
        return "break"

    warn_above_min = bp > context.min_pc or dp > context.min_pc
    return bp, dp, warn_above_min


def _build_binary_vle_pxy_series(
    x0s: List[float],
    temperature: float,
    tcs: List[float],
    vps: List[Optional[float]],
    context: BinaryVlePxyContext,
) -> BinaryVleSeries:
    series = BinaryVleSeries(xs=[], bps=[], dps=[])
    _append_binary_pure_point(temperature, tcs, vps, series)

    for x0 in x0s:
        point = _binary_vle_pxy_point(context, x0)
        if point is None:
            continue
        if point == "break":
            break
        bp, dp, warn_above_min = point
        if warn_above_min:
            logger.warning(
                "mix_vle_pxy: point above one component Pc at x0=%.4f: bp=%.2f, dp=%.2f",
                x0,
                bp,
                dp,
            )
        series.bps.append(max(bp, dp))
        series.dps.append(min(bp, dp))
        series.xs.append(x0)

    return series


def _ternary_critical_points(
    parameters_list: List[List[float]],
) -> Tuple[List[float], List[float]]:
    tcs = []
    pcs = []
    for params in parameters_list:
        tc, pc, _ = critical_points_feos(parameters=params)
        tcs.append(tc)
        pcs.append(pc)
    return tcs, pcs


def _ternary_pure_vps(
    parameters_list: List[List[float]],
    temperature: float,
    tcs: List[float],
) -> List[Optional[float]]:
    vps = []
    for params, tc in zip(parameters_list, tcs):
        try:
            if temperature < tc:
                vps.append(pure_vp_feos(parameters=params, state=[temperature]))
            else:
                vps.append(None)
        except RuntimeError:
            vps.append(None)
    return vps


def _validate_ternary_vle_tx_fixed(
    temperature: float,
    tcs: List[float],
    vps: List[Optional[float]],
):
    if all(temperature >= tc for tc in tcs):
        raise ValueError(
            f"Temperature {temperature} K is above the critical temperature of all components; "
            "VLE calculation is not meaningful.",
        )
    if all(vp is None for vp in vps):
        raise ValueError(
            f"Unable to compute pure vapor pressures at {temperature} K for any component; "
            "VLE calculation is not meaningful.",
        )


def _volatility_rank(vp: Optional[float], tc: float, temperature: float) -> float:
    if vp is not None:
        return vp
    if temperature >= tc:
        return float("inf")
    return float("-inf")


def _validate_ternary_volatility(
    smiles_list: List[str],
    vps: List[Optional[float]],
    tcs: List[float],
    temperature: float,
):
    ranks = [_volatility_rank(vps[i], tcs[i], temperature) for i in range(3)]
    most_volatile_index = int(np.argmax(ranks))
    if most_volatile_index != 0:
        raise ValueError(
            f"More volatile component ({smiles_list[most_volatile_index]}) must be listed first "
            "in ternary P-x calculation",
        )


def _validate_solvent_ratio(solvent_ratio: float):
    if not 0.0 < solvent_ratio < 1.0:
        raise ValueError(
            f"For ternary P-x, solvent ratio must be between 0 and 1, got ratio = {solvent_ratio}"
        )


def _ternary_solvent_split(x1: float, solvent_ratio: float) -> Tuple[float, float]:
    remaining = 1.0 - x1
    x2 = remaining * solvent_ratio
    x3 = remaining * (1.0 - solvent_ratio)
    return x2, x3


def _ternary_vle_tx_point(
    context: TernaryVleContext,
    x1: float,
    x2: float,
    x3: float,
) -> Optional[Tuple[float, float]]:
    try:
        return mix_vp_feos(
            parameters=context.parameters_list,
            state=[context.temperature, 0.0, float(x1), float(x2), float(x3)],
            kij_matrix=context.kij_matrix,
        )
    except RuntimeError as exc:
        logger.debug(
            "mix_ternary_vle_tx_fixed: Runtime Error at x1=%.4f, x2=%.4f, x3=%.4f: %s",
            x1,
            x2,
            x3,
            exc,
        )
        return None
    except BaseException as exc:  # pylint: disable=W0718
        exception_type = type(exc).__name__
        if exception_type == "PanicException":
            logger.warning(
                "mix_ternary_vle_tx_fixed: PanicException at x1=%.4f, x2=%.4f, x3=%.4f: %s",
                x1,
                x2,
                x3,
                exc,
            )
            return None
        logger.exception(
            "mix_ternary_vle_tx_fixed: unexpected %s at x1=%.4f, x2=%.4f, x3=%.4f",
            exception_type,
            x1,
            x2,
            x3,
        )
        raise


def _append_ternary_pressures(
    x1: float,
    bubble_p: float,
    dew_p: float,
    context: TernaryVleContext,
    series: TernaryVleSeries,
):
    if bubble_p > context.max_pc or dew_p > context.max_pc:
        logger.warning(
            "mix_ternary_vle_tx_fixed: breaking from "
            "point above all Pc at x1=%.4f: bp=%.2f, dp=%.2f",
            x1,
            bubble_p,
            dew_p,
        )
        return "break"
    if bubble_p > context.min_pc or dew_p > context.min_pc:
        logger.warning(
            "mix_ternary_vle_tx_fixed: point above at least "
            "one Pc at x1=%.4f: bp=%.2f, dp=%.2f",
            x1,
            bubble_p,
            dew_p,
        )

    if np.isfinite(bubble_p) and np.isfinite(dew_p) and bubble_p > 0.0 and dew_p > 0.0:
        series.x1_values.append(float(x1))
        series.bubble_pressures.append(float(max(bubble_p, dew_p)))
        series.dew_pressures.append(float(min(bubble_p, dew_p)))
    return None


def _build_ternary_vle_series(
    x1_grid: List[float],
    solvent_ratio: float,
    context: TernaryVleContext,
) -> TernaryVleSeries:
    series = TernaryVleSeries(x1_values=[], bubble_pressures=[], dew_pressures=[])

    for x1 in x1_grid:
        x2, x3 = _ternary_solvent_split(x1, solvent_ratio)
        if x2 <= 0.0 or x3 <= 0.0:
            continue

        point = _ternary_vle_tx_point(context, x1, x2, x3)
        if point is None:
            continue

        bubble_p, dew_p = point
        result = _append_ternary_pressures(x1, bubble_p, dew_p, context, series)
        if result == "break":
            break

    return series


def _build_fraction_grid(
    mole_fractions: Optional[List[float]],
    n_points: int,
) -> List[float]:
    x_grid = np.linspace(0.0, 1.0, num=n_points, dtype=np.float64).tolist()
    if mole_fractions:
        x_grid.extend(mole_fractions)
        x_grid = sorted(x_grid)
    return x_grid


def mix_den(params: MixDenParams) -> Tuple[List[float], List[float]]:
    """Calculate mixture density using PC-SAFT EOS."""
    parameters_list = [
        predict_pcsaft_parameters(smiles) for smiles in params.smiles_list
    ]
    temperatures = np.linspace(
        params.min_temp, params.max_temp, num=params.npoints
    ).tolist()

    densities = [
        mix_den_feos(
            parameters=parameters_list,
            state=[T, params.pressure] + params.mole_fractions,
            kij_matrix=params.kij_matrix,
        )
        for T in temperatures
    ]
    return temperatures, densities


def mix_vp(params: MixVpParams) -> Tuple[List[float], List[float], List[float]]:
    """Calculate mixture vapor pressure using PC-SAFT EOS."""
    parameters_list = [
        predict_pcsaft_parameters(smiles) for smiles in params.smiles_list
    ]
    temperatures = np.linspace(
        params.min_temp, params.max_temp, num=params.npoints
    ).tolist()

    buble_points = []
    dew_point = []
    valid_temperatures = []
    for temp in temperatures:
        try:
            bp, dp = mix_vp_feos(
                parameters=parameters_list,
                state=[temp, 0] + params.mole_fractions,
                kij_matrix=params.kij_matrix,
            )
            if bp > dp:
                buble_points.append(bp)
                dew_point.append(dp)
            else:
                buble_points.append(dp)
                dew_point.append(bp)
            valid_temperatures.append(temp)
        except RuntimeError as exc:
            logger.debug("mix_vp: Runtime Error at temperature=%.4f: %s", temp, exc)
        except BaseException as exc:  # pylint: disable=W0718
            exception_type = type(exc).__name__
            if exception_type == "PanicException":
                logger.warning(
                    "mix_vp: PanicException at temperature=%.4f: %s", temp, exc
                )
            else:
                logger.exception(
                    "mix_vp: unexpected %s at temperature=%.4f",
                    exception_type,
                    temp,
                )
                raise

    return valid_temperatures, buble_points, dew_point


def mix_vle(
    smiles_list: List[str], kij_matrix: List[List[float]], pressure: float, npoints: int
) -> Dict[str, List[float]]:
    "Calculate mixture VLE (T-x-y) using PC-SAFT EOS"
    parameters_list = [predict_pcsaft_parameters(smiles) for smiles in smiles_list]

    return mix_vle_diagram_feos(
        parameters=parameters_list,
        state=[pressure],
        kij_matrix=kij_matrix,
        npoints=npoints,
    )


def mix_vle_pxy(
    smiles_list: List[str],
    kij_matrix: List[List[float]],
    temperature: float,
    npoints: int,
    mole_fractions: Optional[List[float]] = None,
) -> Tuple[List[float], List[float], List[float]]:
    "Calculate mixture VLE (P-x-y) using PC-SAFT EOS"
    parameters_list = [predict_pcsaft_parameters(smiles) for smiles in smiles_list]
    x0s = _build_fraction_grid(mole_fractions, npoints)
    tcs, pcs = _binary_critical_points(parameters_list)
    vps = _binary_pure_vps(parameters_list, temperature, tcs)

    _validate_binary_vle_pxy(temperature, tcs, vps)
    if not _binary_more_volatile_is_first(vps, tcs, temperature):
        raise ValueError(
            f"More volatile component ({smiles_list[1]}) must be listed first in P-x-y calculation",
        )

    min_pc = min(pcs)
    max_pc = max(pcs)

    context = BinaryVlePxyContext(
        parameters_list=parameters_list,
        kij_matrix=kij_matrix,
        temperature=temperature,
        min_pc=min_pc,
        max_pc=max_pc,
    )
    series = _build_binary_vle_pxy_series(x0s, temperature, tcs, vps, context)

    return series.xs, series.bps, series.dps


def mix_lle(
    params: MixLLEParams,
) -> Dict[str, List[float]]:
    "Calculate mixture LLE using PC-SAFT EOS"
    parameters_list = [
        predict_pcsaft_parameters(smiles) for smiles in params.smiles_list
    ]

    return mix_lle_diagram_feos(
        parameters=parameters_list,
        state=[params.temperature, params.pressure, *params.mole_fractions],
        kij_matrix=params.kij_matrix,
        npoints=params.npoints,
    )


def _get_ternary_lle_data(
    params: List[List[float]],
    state: List[float],
    kij_matrix: List[List[float]],
    npoints: int,
) -> Dict[str, List[float]]:
    t, p = state  # Temperatura (K) e pressão (Pa)

    def _grid(n_pts: int = npoints):
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
            if lle["density liquid"][0] > lle["density vapor"][0]:
                ternary_data["x0"].extend(lle["x0"])
                ternary_data["x1"].extend(lle["x1"])
                ternary_data["x2"].extend(lle["x2"])
                ternary_data["y0"].extend(lle["y0"])
                ternary_data["y1"].extend(lle["y1"])
                ternary_data["y2"].extend(lle["y2"])
            else:
                ternary_data["x0"].extend(lle["y0"])
                ternary_data["x1"].extend(lle["y1"])
                ternary_data["x2"].extend(lle["y2"])
                ternary_data["y0"].extend(lle["x0"])
                ternary_data["y1"].extend(lle["x1"])
                ternary_data["y2"].extend(lle["x2"])
        return ternary_data

    x1, x2, x3, mask = _grid()
    return _collect_tie_lines(x1, x2, x3, mask)


def mix_ternary_lle(
    smiles_list: List[str],
    kij_matrix: List[List[float]],
    temperature: float,
    pressure: float,
    npoints: int,
) -> Dict[str, List[float]]:
    "Calculate ternary LLE/VLE using PC-SAFT EOS"
    parameters_list = [predict_pcsaft_parameters(smiles) for smiles in smiles_list]

    return _get_ternary_lle_data(
        params=parameters_list,
        state=[temperature, pressure],
        kij_matrix=kij_matrix,
        npoints=npoints,
    )


def mix_ternary_vle_tx_fixed(
    params: TernaryVleTxParams,
) -> Tuple[List[float], List[float], List[float]]:
    """
    Calculate ternary isothermal VLE curve (P-x) at fixed solvent ratio.

    solvent_ratio = x2 / (x2 + x3). The first component is scanned in composition.
    Uses values from the provided params dataclass.
    """
    parameters_list = [
        predict_pcsaft_parameters(smiles) for smiles in params.smiles_list
    ]
    tcs, pcs = _ternary_critical_points(parameters_list)
    vps = _ternary_pure_vps(parameters_list, params.temperature, tcs)

    _validate_ternary_vle_tx_fixed(params.temperature, tcs, vps)
    _validate_ternary_volatility(params.smiles_list, vps, tcs, params.temperature)
    _validate_solvent_ratio(params.solvent_ratio)

    context = TernaryVleContext(
        parameters_list=parameters_list,
        kij_matrix=params.kij_matrix,
        temperature=params.temperature,
        min_pc=min(pcs),
        max_pc=max(pcs),
    )
    x1_grid = _build_fraction_grid(params.mole_fractions, n_points=params.npoints)
    series = _build_ternary_vle_series(x1_grid, params.solvent_ratio, context)

    return series.x1_values, series.bubble_pressures, series.dew_pressures
