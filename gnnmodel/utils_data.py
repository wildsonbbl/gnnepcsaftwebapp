"Experimental data utilities"

import os.path as osp
from typing import List, Optional, Tuple

import polars as pl
from gnnepcsaft_mcp_server.utils import smilestoinchi
from numpy import float64
from numpy.typing import NDArray

application_path = osp.dirname(osp.abspath(__file__))

# Tolerance constants used across data retrieval functions
# Fraction tolerance matches rounding to 4 decimal places used in grouping
TOL_FRACTION = 1e-4
# Solvent ratio tolerance matches rounding to 2 decimal places
TOL_SOLVENT_RATIO = 0.01
# Temperature tolerance for coarse matching (K)
TOL_TEMP = 0.1
# Small tolerance for pressure/temperature fine matching
TOL_PRESSURE_TEMP = 0.01


def _read_parquet_if_exists(path_list: List[str]) -> Optional[pl.DataFrame]:
    if any(not osp.exists(path) for path in path_list):
        return None
    return pl.scan_parquet(
        path_list, extra_columns="ignore", missing_columns="insert"
    ).collect()


def _filter_binary_pair(
    dframe: pl.DataFrame,
    inchi1: str,
    inchi2: str,
    col_x1: str,
    col_x2: str,
) -> pl.DataFrame:
    return (
        dframe.filter(
            ((pl.col("inchi1") == inchi1) & (pl.col("inchi2") == inchi2))
            | ((pl.col("inchi1") == inchi2) & (pl.col("inchi2") == inchi1))
        )
        .with_columns(
            pl.when(pl.col("inchi1") == inchi1)
            .then(pl.col(col_x1))
            .otherwise(pl.col(col_x2))
            .alias("x_c1"),
        )
        .with_columns(
            pl.col("x_c1").round(4).alias("x_approx"),
            pl.col("T_K").round(1).alias("T_approx"),
        )
        .unique()
    )


def _build_binary_rho_data(
    inchi1: str,
    inchi2: str,
) -> Optional[NDArray[float64]]:
    rho_bin = pl.read_parquet(osp.join(application_path, "_data", "rho_binary.parquet"))
    filtered = _filter_binary_pair(
        rho_bin, inchi1, inchi2, "mole_fraction_c1", "mole_fraction_c2"
    )
    if filtered.height == 0:
        return None
    return (
        filtered.group_by(["P_kPa", "x_approx"])
        .agg(
            pl.col("T_K").min().alias("T_min"),
            pl.col("T_K").max().alias("T_max"),
            pl.len().alias("count"),
        )
        .sort(["P_kPa", "x_approx"])
        .to_numpy()
    )


def _build_binary_bubble_data(
    inchi1: str,
    inchi2: str,
) -> Optional[NDArray[float64]]:
    path_vle = osp.join(application_path, "_data", "vle_binary.parquet")
    path_vp = osp.join(application_path, "_data", "vp_binary.parquet")
    df = _read_parquet_if_exists([path_vle, path_vp])

    if df is not None:
        filtered = _filter_binary_pair(
            df, inchi1, inchi2, "mole_fraction_c1p2", "mole_fraction_c2p2"
        )
        if filtered.height > 0:
            return (
                filtered.group_by("x_approx")
                .agg(
                    pl.col("T_K").min().alias("T_min"),
                    pl.col("T_K").max().alias("T_max"),
                    pl.len().alias("count"),
                )
                .sort("x_approx")
                .to_numpy()
            )
    return None


def _build_binary_lle_data(
    inchi1: str,
    inchi2: str,
) -> Optional[NDArray[float64]]:
    path_lle = osp.join(application_path, "_data", "lle_binary.parquet")
    path_lle_temp = osp.join(application_path, "_data", "lle_binary_temp.parquet")
    df = _read_parquet_if_exists([path_lle, path_lle_temp])

    if df is not None:
        filtered = _filter_binary_pair(
            df, inchi1, inchi2, "mole_fraction_c1", "mole_fraction_c2"
        )
        if filtered.height > 0:
            return (
                filtered.group_by("P_kPa")
                .agg(
                    pl.col("T_K").min().alias("T_min"),
                    pl.col("T_K").max().alias("T_max"),
                    pl.len().alias("count"),
                )
                .sort("P_kPa")
                .to_numpy()
            )
    return None


def _build_binary_vle_data(
    inchi1: str,
    inchi2: str,
) -> Tuple[Optional[NDArray[float64]], Optional[NDArray[float64]]]:
    path_vle = osp.join(application_path, "_data", "vle_binary.parquet")
    path_vp = osp.join(application_path, "_data", "vp_binary.parquet")
    df = _read_parquet_if_exists([path_vle, path_vp])

    if df is not None:
        filtered = _filter_binary_pair(
            df, inchi1, inchi2, "mole_fraction_c1p2", "mole_fraction_c2p2"
        )
        if filtered.height > 0:
            return (
                filtered.group_by("P_kPa")
                .agg(
                    pl.col("T_K").min().alias("T_min"),
                    pl.col("T_K").max().alias("T_max"),
                    pl.len().alias("count"),
                )
                .sort("P_kPa")
                .to_numpy(),
                filtered.group_by("T_approx")
                .agg(
                    pl.col("P_kPa").min().alias("P_min"),
                    pl.col("P_kPa").max().alias("P_max"),
                    pl.len().alias("count"),
                )
                .sort("T_approx")
                .to_numpy(),
            )

    return None, None


def _filter_ternary_set(
    dframe: pl.DataFrame, target_set: list, col_x1: str, col_x2: str, col_x3: str
) -> pl.DataFrame:
    return (
        dframe.filter(
            pl.col("inchi1").is_in(target_set)
            & pl.col("inchi2").is_in(target_set)
            & pl.col("inchi3").is_in(target_set)
        )
        .with_columns(
            pl.when(pl.col("inchi1") == target_set[0])
            .then(pl.col(col_x1))
            .otherwise(
                pl.when(pl.col("inchi2") == target_set[0])
                .then(pl.col(col_x2))
                .otherwise(pl.col(col_x3))
            )
            .alias("x_mapped_1"),
            pl.when(pl.col("inchi1") == target_set[1])
            .then(pl.col(col_x1))
            .otherwise(
                pl.when(pl.col("inchi2") == target_set[1])
                .then(pl.col(col_x2))
                .otherwise(pl.col(col_x3))
            )
            .alias("x_mapped_2"),
            pl.when(pl.col("inchi1") == target_set[2])
            .then(pl.col(col_x1))
            .otherwise(
                pl.when(pl.col("inchi2") == target_set[2])
                .then(pl.col(col_x2))
                .otherwise(pl.col(col_x3))
            )
            .alias("x_mapped_3"),
        )
        .with_columns(
            pl.col("x_mapped_1").round(4).alias("x_approx_1"),
            pl.col("x_mapped_2").round(4).alias("x_approx_2"),
            (pl.col("x_mapped_2") / (pl.col("x_mapped_2") + pl.col("x_mapped_3")))
            .round(2)
            .alias("solvent_ratio"),
        )
        .unique()
    )


def _build_ternary_rho_data(
    target_set: list,
) -> Optional[NDArray[float64]]:
    path_rho = osp.join(application_path, "_data", "rho_ternary.parquet")
    df = _read_parquet_if_exists([path_rho])
    if df is None:
        return None

    filtered = _filter_ternary_set(
        df, target_set, "mole_fraction_c1", "mole_fraction_c2", "mole_fraction_c3"
    )
    if filtered.height > 0:
        return (
            filtered.group_by(["P_kPa", "x_approx_1", "x_approx_2"])
            .agg(
                pl.col("T_K").min().alias("T_min"),
                pl.col("T_K").max().alias("T_max"),
                pl.len().alias("count"),
            )
            .sort(["P_kPa", "x_approx_1", "x_approx_2"])
            .to_numpy()
        )
    return None


def _build_ternary_lle_data(target_set: list) -> Optional[NDArray[float64]]:
    path_lle = osp.join(application_path, "_data", "lle_ternary.parquet")
    path_lle_mass = osp.join(application_path, "_data", "lle_mass_ternary.parquet")
    df = _read_parquet_if_exists([path_lle, path_lle_mass])
    if df is None:
        return None

    filtered = _filter_ternary_set(
        df, target_set, "mole_fraction_c1", "mole_fraction_c2", "mole_fraction_c3"
    )
    if filtered.height == 0:
        return None

    return (
        filtered.group_by(["P_kPa", "T_K"])
        .agg(pl.len().alias("count"))
        .sort(["P_kPa", "T_K"])
        .to_numpy()
    )


def _build_ternary_vle_data(
    target_set: list,
) -> Tuple[Optional[NDArray[float64]], Optional[NDArray[float64]]]:
    path_vle = osp.join(application_path, "_data", "vle_ternary.parquet")
    path_vp = osp.join(application_path, "_data", "vp_ternary.parquet")
    df = _read_parquet_if_exists([path_vle, path_vp])
    if df is None:
        return None, None

    filtered = _filter_ternary_set(
        df, target_set, "mole_fraction_c1p2", "mole_fraction_c2p2", "mole_fraction_c3p2"
    )
    if filtered.height == 0:
        return None, None

    vle_data = (
        filtered.group_by(["P_kPa", "T_K"])
        .agg(pl.len().alias("count"))
        .sort(["P_kPa", "T_K"])
        .to_numpy()
    )

    vle_tx_data = (
        filtered.group_by(["T_K", "solvent_ratio"])
        .agg(
            pl.col("P_kPa").min().alias("P_min"),
            pl.col("P_kPa").max().alias("P_max"),
            pl.len().alias("count"),
        )
        .sort(["T_K", "solvent_ratio"])
        .to_numpy()
    )

    return vle_data, vle_tx_data


def default_mixture_output_args():
    """Return the default output_args dict for mixture plots."""
    return {
        "rho_data": None,
        "bubble_data": None,
        "lle_data": None,
        "vle_data": None,
        "vle_pxy_data": None,
        "rho_data_t": None,
        "lle_data_t": None,
        "vle_data_t": None,
        "vle_tx_data_t": None,
        "preds": [],
    }


def retrieve_rho_pure_data(smiles: str, pressure: float) -> Optional[NDArray[float64]]:
    "retrieve density data for plots"

    df = pl.read_parquet(osp.join(application_path, "_data", "rho_pure.parquet"))

    return (
        df.filter(
            pl.col("inchi1") == smilestoinchi(smiles), pl.col("P_kPa") == pressure
        )
        .select(
            pl.col("T_K"),
            (pl.col("rho") * 1000 / pl.col("molweight1")),
        )
        .to_numpy()
    )


def retrieve_vp_pure_data(smiles: str) -> Optional[NDArray[float64]]:
    "retrieve vapor pressure data for plots"

    df = pl.read_parquet(osp.join(application_path, "_data", "vp_pure.parquet"))

    return (
        df.filter(pl.col("inchi1") == smilestoinchi(smiles))
        .select("T_K", "VP_kPa")
        .to_numpy()
    )


def retrieve_st_pure_data(smiles: str) -> Optional[NDArray[float64]]:
    "retrieve surface tension (N/m) data for plots"

    df = pl.read_parquet(osp.join(application_path, "_data", "st_pure.parquet"))

    return (
        df.filter(pl.col("inchi1") == smilestoinchi(smiles))
        .select("T_K", "st")
        .to_numpy()
    )


def retrieve_available_data_pure(
    smiles: str,
) -> Tuple[
    Optional[NDArray[float64]],
    Tuple[
        Optional[NDArray[float64]],
        Optional[NDArray[float64]],
        int,
    ],
    Tuple[
        Optional[NDArray[float64]],
        Optional[NDArray[float64]],
        int,
    ],
]:
    "retrieve available pure data for smiles"

    rho_pure = pl.read_parquet(osp.join(application_path, "_data", "rho_pure.parquet"))
    vp_pure = pl.read_parquet(osp.join(application_path, "_data", "vp_pure.parquet"))
    st_pure = pl.read_parquet(osp.join(application_path, "_data", "st_pure.parquet"))

    inchi = smilestoinchi(smiles)

    rho_filtered = rho_pure.filter(pl.col("inchi1") == inchi)
    if rho_filtered.height > 0:
        pure_data = (
            rho_filtered.select("T_K", "P_kPa")
            .group_by(pl.col("P_kPa"))
            .agg(
                pl.col("T_K").min().alias("T_min"),
                pl.col("T_K").max().alias("T_max"),
                pl.len().alias("count"),
            )
            .sort(pl.col("P_kPa"))
            .to_numpy()
        )
    else:
        pure_data = None

    vp_filtered = vp_pure.filter(pl.col("inchi1") == inchi)
    if vp_filtered.height > 0:
        vp_data = vp_filtered.select("T_K").to_numpy()
        vp_range = (vp_data.min(), vp_data.max(), vp_filtered.height)
    else:
        vp_range = (None, None, 0)

    st_filtered = st_pure.filter(pl.col("inchi1") == inchi)
    if st_filtered.height > 0:
        st_data = st_filtered.select("T_K").to_numpy()
        st_range = (st_data.min(), st_data.max(), st_filtered.height)
    else:
        st_range = (None, None, 0)

    return pure_data, vp_range, st_range


def retrieve_rho_binary_data(
    smiles_list: list, pressure: float, x1: float
) -> Optional[NDArray[float64]]:
    "retrieve binary density data"
    if len(smiles_list) != 2:
        return None

    df = pl.read_parquet(osp.join(application_path, "_data", "rho_binary.parquet"))
    i1, i2 = smilestoinchi(smiles_list[0]), smilestoinchi(smiles_list[1])

    # Tolerance
    tol_x = TOL_FRACTION

    # Normalize x1 to strictly match input order
    # If file has (i1, i2) -> use mole_fraction_c1
    # If file has (i2, i1) -> use mole_fraction_c2
    filtered = _filter_binary_pair(
        df, i1, i2, "mole_fraction_c1", "mole_fraction_c2"
    ).filter(
        pl.col("P_kPa").is_close(pressure)
        & pl.col("x_c1").is_between(x1 - tol_x, x1 + tol_x)
    )

    if filtered.height == 0:
        return None

    return filtered.select(
        pl.col("T_K"),
        pl.col("rho")
        * 1000
        / (
            pl.col("molweight1") * pl.col("mole_fraction_c1")
            + pl.col("molweight2") * (1 - pl.col("mole_fraction_c1"))
        ),
    ).to_numpy()


def retrieve_bubble_pressure_data(
    smiles_list: list, x1: float
) -> Optional[NDArray[float64]]:
    "retrieve binary bubble point pressure data (P-T at constant x)"
    if len(smiles_list) != 2:
        return None

    path_vp = osp.join(application_path, "_data", "vp_binary.parquet")
    path_vle = osp.join(application_path, "_data", "vle_binary.parquet")
    df = _read_parquet_if_exists([path_vp, path_vle])
    if df is None:
        return None
    i1, i2 = smilestoinchi(smiles_list[0]), smilestoinchi(smiles_list[1])

    tol_x = TOL_FRACTION

    filtered = _filter_binary_pair(
        df, i1, i2, "mole_fraction_c1p2", "mole_fraction_c2p2"
    ).filter(pl.col("x_c1").is_close(x1, abs_tol=tol_x))

    if filtered.height == 0:
        return None

    # Return T and P_kPa
    return filtered.select("T_K", "P_kPa").to_numpy()


def retrieve_vle_binary_data(
    smiles_list: list, pressure: float
) -> Optional[NDArray[float64]]:
    """
    retrieve binary VLE data (T-x-y)
    """
    if len(smiles_list) != 2:
        return None

    i1, i2 = smilestoinchi(smiles_list[0]), smilestoinchi(smiles_list[1])

    path_vle = osp.join(application_path, "_data", "vle_binary.parquet")
    path_vp = osp.join(application_path, "_data", "vp_binary.parquet")
    df = _read_parquet_if_exists([path_vle, path_vp])
    if df is not None:
        filtered = _filter_binary_pair(
            df, i1, i2, "mole_fraction_c1p2", "mole_fraction_c2p2"
        ).filter(pl.col("P_kPa").is_close(pressure))

        if filtered.height > 0:
            return filtered.select("T_K", "x_c1").to_numpy()

    return None


def retrieve_vle_pxy_binary_data(
    smiles_list: list, temperature: float
) -> Optional[NDArray[float64]]:
    """
    retrieve binary VLE data (P-x-y) at constant T.
    """
    if len(smiles_list) != 2:
        return None

    i1, i2 = smilestoinchi(smiles_list[0]), smilestoinchi(smiles_list[1])
    tol_t = TOL_TEMP  # Tolerance for temperature

    path_vle = osp.join(application_path, "_data", "vle_binary.parquet")
    path_vp = osp.join(application_path, "_data", "vp_binary.parquet")
    df = _read_parquet_if_exists([path_vle, path_vp])
    if df is not None:
        filtered = _filter_binary_pair(
            df, i1, i2, "mole_fraction_c1p2", "mole_fraction_c2p2"
        ).filter(pl.col("T_K").is_close(temperature, abs_tol=tol_t))
        if filtered.height > 0:
            return filtered.select("P_kPa", "x_c1").to_numpy()

    return None


def retrieve_lle_binary_data(
    smiles_list: list, pressure: float
) -> Optional[NDArray[float64]]:
    "retrieve binary LLE data (T-x-x)"
    if len(smiles_list) != 2:
        return None

    path = osp.join(application_path, "_data", "lle_binary.parquet")
    path_temp = osp.join(application_path, "_data", "lle_binary_temp.parquet")

    df = _read_parquet_if_exists([path, path_temp])
    if df is None:
        return None
    i1, i2 = smilestoinchi(smiles_list[0]), smilestoinchi(smiles_list[1])

    filtered = _filter_binary_pair(
        df, i1, i2, "mole_fraction_c1", "mole_fraction_c2"
    ).filter(pl.col("P_kPa").is_close(pressure))

    if filtered.height == 0:
        return None

    return filtered.select("T_K", "x_c1").to_numpy()


def retrieve_available_data_binary(smiles_list: list) -> Tuple[
    Optional[NDArray[float64]],
    Optional[NDArray[float64]],
    Optional[NDArray[float64]],
    Optional[NDArray[float64]],
    Optional[NDArray[float64]],
]:
    "retrieve available binary data"
    if len(smiles_list) != 2:
        return None, None, None, None, None

    i1, i2 = smilestoinchi(smiles_list[0]), smilestoinchi(smiles_list[1])

    rho_data = _build_binary_rho_data(i1, i2)
    bubble_data = _build_binary_bubble_data(i1, i2)
    lle_data = _build_binary_lle_data(i1, i2)
    vle_data, vle_pxy_data = _build_binary_vle_data(i1, i2)

    return rho_data, bubble_data, lle_data, vle_data, vle_pxy_data


def retrieve_available_data_ternary(
    smiles_list: list,
) -> Tuple[
    Optional[NDArray[float64]],
    Optional[NDArray[float64]],
    Optional[NDArray[float64]],
    Optional[NDArray[float64]],
]:
    "retrieve available ternary data"
    if len(smiles_list) != 3:
        return None, None, None, None

    target_set = [
        smilestoinchi(smiles_list[0]),
        smilestoinchi(smiles_list[1]),
        smilestoinchi(smiles_list[2]),
    ]

    rho_data = _build_ternary_rho_data(target_set)
    lle_data = _build_ternary_lle_data(target_set)
    vle_data, vle_tx_data = _build_ternary_vle_data(target_set)

    return rho_data, lle_data, vle_data, vle_tx_data


def retrieve_rho_ternary_data(
    smiles_list: list, pressure: float, x1: float, x2: float
) -> Optional[NDArray[float64]]:
    "retrieve ternary density data"
    if len(smiles_list) != 3:
        return None

    path_rho = osp.join(application_path, "_data", "rho_ternary.parquet")

    target_set = [
        smilestoinchi(smiles_list[0]),
        smilestoinchi(smiles_list[1]),
        smilestoinchi(smiles_list[2]),
    ]

    df = _read_parquet_if_exists([path_rho])
    if df is None:
        return None

    # Tolerance
    tol_x = TOL_FRACTION

    filtered = _filter_ternary_set(
        df, target_set, "mole_fraction_c1", "mole_fraction_c2", "mole_fraction_c3"
    ).filter(
        pl.col("P_kPa").is_close(pressure)
        & (pl.col("x_mapped_1").is_close(x1, abs_tol=tol_x))
        & (pl.col("x_mapped_2").is_close(x2, abs_tol=tol_x))
    )

    if filtered.height == 0:
        return None

    # molar density = mass_rho * 1000 / avg_mw
    return filtered.select(
        pl.col("T_K"),
        pl.col("rho")
        * 1000.0
        / (
            pl.col("mole_fraction_c1") * pl.col("molweight1")
            + pl.col("mole_fraction_c2") * pl.col("molweight2")
            + pl.col("mole_fraction_c3") * pl.col("molweight3")
        ),
    ).to_numpy()


def retrieve_lle_ternary_data(
    smiles_list: list, pressure: float, temperature: float
) -> Optional[NDArray[float64]]:
    "retrieve ternary lle data (tie lines/binodal points)"
    if len(smiles_list) != 3:
        return None

    path_lle = osp.join(application_path, "_data", "lle_ternary.parquet")
    path_lle_mass = osp.join(application_path, "_data", "lle_mass_ternary.parquet")
    df = _read_parquet_if_exists([path_lle, path_lle_mass])
    if df is None:
        return None

    target_set = [
        smilestoinchi(smiles_list[0]),
        smilestoinchi(smiles_list[1]),
        smilestoinchi(smiles_list[2]),
    ]

    tol = TOL_PRESSURE_TEMP
    return (
        _filter_ternary_set(
            df, target_set, "mole_fraction_c1", "mole_fraction_c2", "mole_fraction_c3"
        )
        .filter(
            (pl.col("P_kPa").is_close(pressure, abs_tol=tol))
            & (pl.col("T_K").is_close(temperature, abs_tol=tol))
        )
        .select("x_mapped_1", "x_mapped_2")
        .to_numpy()
    )


def retrieve_vle_ternary_data(
    smiles_list: list, pressure: float, temperature: float
) -> Optional[NDArray[float64]]:
    "retrieve ternary vle data (liquid phase composition points)"
    if len(smiles_list) != 3:
        return None

    path_vle = osp.join(application_path, "_data", "vle_ternary.parquet")
    path_vp = osp.join(application_path, "_data", "vp_ternary.parquet")
    df = _read_parquet_if_exists([path_vle, path_vp])
    if df is None:
        return None

    target_set = [
        smilestoinchi(smiles_list[0]),
        smilestoinchi(smiles_list[1]),
        smilestoinchi(smiles_list[2]),
    ]

    tol = TOL_PRESSURE_TEMP
    return (
        _filter_ternary_set(
            df,
            target_set,
            "mole_fraction_c1p2",
            "mole_fraction_c2p2",
            "mole_fraction_c3p2",
        )
        .filter(
            pl.col("P_kPa").is_close(pressure, abs_tol=tol)
            & pl.col("T_K").is_close(temperature, abs_tol=tol)
        )
        .select("x_mapped_1", "x_mapped_2")
        .to_numpy()
    )


def retrieve_vle_ternary_tx_fixed_data(
    smiles_list: list, temperature: float, solvent_ratio: float
) -> Optional[NDArray[float64]]:
    """
    retrieve ternary VLE data for isothermal P-x analysis with fixed solvent ratio.

    solvent_ratio = x2 / (x2 + x3) in liquid phase composition (p2 columns).
    """
    if len(smiles_list) != 3:
        return None

    path_vle = osp.join(application_path, "_data", "vle_ternary.parquet")
    path_vp = osp.join(application_path, "_data", "vp_ternary.parquet")
    df = _read_parquet_if_exists([path_vle, path_vp])
    if df is None:
        return None

    target_set = [
        smilestoinchi(smiles_list[0]),
        smilestoinchi(smiles_list[1]),
        smilestoinchi(smiles_list[2]),
    ]

    tol_t = TOL_TEMP
    tol_ratio = TOL_SOLVENT_RATIO
    return (
        _filter_ternary_set(
            df,
            target_set,
            "mole_fraction_c1p2",
            "mole_fraction_c2p2",
            "mole_fraction_c3p2",
        )
        .filter(
            pl.col("T_K").is_close(temperature, abs_tol=tol_t)
            & (pl.col("x_mapped_2") + pl.col("x_mapped_3") > 1e-12)
            & pl.col("solvent_ratio").is_close(solvent_ratio, abs_tol=tol_ratio)
        )
        .select("x_mapped_1", "P_kPa")
        .to_numpy()
    )
