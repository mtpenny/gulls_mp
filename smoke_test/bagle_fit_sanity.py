"""Joint BAGLE photometric+astrometric fit sanity checks."""
from __future__ import annotations

from dataclasses import dataclass, field
import json
import math
import os
import re
from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence, Tuple

import numpy as np
import pandas as pd

from .errors import SmokeTestError

try:
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    _HAS_ASTROPY = True
except Exception:  # pragma: no cover - optional dependency guard
    _HAS_ASTROPY = False


MAS_PER_ARCSEC = 1000.0
MAS_PER_DEG = 3600.0 * 1000.0
LOG10_FACTOR = 2.5 / math.log(10.0)


@dataclass
class BagleJointFitSummary:
    """Summary of the BAGLE joint-fit sanity check."""

    run_dir: Path
    lc_file: Path
    event_id: int
    subrun: int
    field: int
    out_chi2_single_lens: float
    n_phot_points: int
    n_ast_points: int
    fit_reduced_chi2: float
    mu_ref_e: float
    mu_ref_n: float
    mu_fit_e: float
    mu_fit_n: float
    mu_amp_ref: float
    mu_amp_fit: float
    mu_dir_diff_deg: float
    piE_ref_e: float
    piE_ref_n: float
    piE_fit_e: float
    piE_fit_n: float
    piE_amp_ref: float
    piE_amp_fit: float
    piE_dir_diff_deg: float
    true_ast_rms_mas: float
    true_ast_sigma_equiv: float
    plot_path: Path
    result_json_path: Path
    warnings: List[str] = field(default_factory=list)


def _angle_diff_deg(lhs: float, rhs: float) -> float:
    return (lhs - rhs + 180.0) % 360.0 - 180.0


def _uniform_sample(indices: np.ndarray, max_points: int) -> np.ndarray:
    if max_points is None or max_points <= 0:
        return indices
    if len(indices) <= max_points:
        return indices
    pick = np.linspace(0, len(indices) - 1, max_points, dtype=int)
    return indices[np.unique(pick)]


def _derive_event_key(lc_file: Path) -> Tuple[int, int, int | None] | None:
    stem = lc_file.stem.split(".", 1)[0]

    match = re.search(r"_(\d+)_(\d+)_(\d+)$", stem)
    if match:
        subrun = int(match.group(1))
        field = int(match.group(2))
        event = int(match.group(3))
        return event, subrun, field

    match = re.search(r"_(\d+)_(\d+)$", stem)
    if match:
        subrun = int(match.group(1))
        event = int(match.group(2))
        return event, subrun, None

    return None


def _find_lc_for_event(run_dir: Path, event_id: int, subrun: int, field: int) -> Path:
    pattern = f"*_{subrun}_{field}_{event_id}.*.lc"
    matches = sorted(run_dir.glob(pattern))
    if not matches:
        matches = sorted(run_dir.rglob(pattern))
    if not matches:
        raise SmokeTestError(
            f"BAGLE sanity: cannot find .lc for EventID={event_id}, SubRun={subrun}, Field={field} under {run_dir}"
        )
    det = [path for path in matches if path.suffixes[-2:] == [".det", ".lc"]]
    if det:
        return det[0]
    return matches[0]


def _safe_get(row: Mapping[str, Any], key: str) -> float:
    value = row.get(key, math.nan)
    try:
        return float(value)
    except Exception:
        return math.nan


def _parse_astrometry_frame(lc_file: Path) -> Tuple[float | None, float | None]:
    with lc_file.open(encoding="utf-8") as handle:
        for raw in handle:
            if not raw.startswith("#"):
                break
            if not raw.startswith("#Astrometry_Frame:"):
                continue
            ra_deg = None
            dec_deg = None
            for token in raw.strip().split():
                if token.startswith("RA_deg="):
                    try:
                        ra_deg = float(token.split("=", 1)[1])
                    except ValueError:
                        ra_deg = None
                elif token.startswith("Dec_deg="):
                    try:
                        dec_deg = float(token.split("=", 1)[1])
                    except ValueError:
                        dec_deg = None
            if ra_deg is not None and dec_deg is not None:
                return ra_deg, dec_deg
    return None, None


def _parse_astrometry_transform(
    lc_file: Path,
) -> Tuple[float, float, float, float] | None:
    with lc_file.open(encoding="utf-8") as handle:
        for raw in handle:
            if not raw.startswith("#"):
                break
            if not raw.startswith("#Astrometry_Transform:"):
                continue
            coeff_tokens = re.findall(r"\(([-+0-9.eE]+)\)", raw)
            if len(coeff_tokens) < 4:
                continue
            try:
                a11 = float(coeff_tokens[0])
                a12 = float(coeff_tokens[1])
                a21 = float(coeff_tokens[2])
                a22 = float(coeff_tokens[3])
            except ValueError:
                continue
            return a11, a12, a21, a22
    return None


def _parse_astrometry_bagle_model_frame(lc_file: Path) -> str | None:
    with lc_file.open(encoding="utf-8") as handle:
        for raw in handle:
            if not raw.startswith("#"):
                break
            if not raw.startswith("#Astrometry_BAGLE:"):
                continue
            for token in raw.strip().split():
                if token.startswith("model_frame="):
                    return token.split("=", 1)[1].strip().lower()
            return None
    return None


def _pm_gal_to_icrs(l_deg: float, b_deg: float, mu_l: float, mu_b: float) -> Tuple[float, float]:
    if not _HAS_ASTROPY:
        raise RuntimeError("astropy not available")
    coord = SkyCoord(
        l=l_deg * u.deg,
        b=b_deg * u.deg,
        pm_l_cosb=mu_l * u.mas / u.yr,
        pm_b=mu_b * u.mas / u.yr,
        frame="galactic",
    )
    icrs = coord.icrs
    return (
        float(icrs.pm_ra_cosdec.to_value(u.mas / u.yr)),
        float(icrs.pm_dec.to_value(u.mas / u.yr)),
    )


def _set_uniform_prior(
    fitter: Any,
    model_fitter: Any,
    name: str,
    lo: float,
    hi: float,
) -> None:
    if not math.isfinite(lo) or not math.isfinite(hi):
        return
    if hi <= lo:
        return
    fitter.priors[name] = model_fitter.make_gen(lo, hi)


def _vec_from_model_attr(
    model_obj: Any,
    attr_name: str,
) -> np.ndarray | None:
    if not hasattr(model_obj, attr_name):
        return None
    raw = getattr(model_obj, attr_name)
    arr = np.asarray(raw, dtype=float).reshape(-1)
    if arr.size < 2 or not np.all(np.isfinite(arr[:2])):
        return None
    return arr[:2].copy()


def _get_model_astrometry(
    model_obj: Any,
    t_eval: np.ndarray,
    *,
    lens_relative_astrometry: bool,
) -> np.ndarray:
    ast = np.asarray(model_obj.get_astrometry(t_eval), dtype=float)
    if ast.ndim != 2 or ast.shape[1] < 2:
        raise RuntimeError(f"unexpected BAGLE astrometry shape {ast.shape}")
    ast_xy = ast[:, :2]
    if not lens_relative_astrometry:
        return ast_xy
    lens_ast = np.asarray(model_obj.get_lens_astrometry(t_eval), dtype=float)
    if lens_ast.ndim != 2 or lens_ast.shape[1] < 2:
        raise RuntimeError(f"unexpected BAGLE lens astrometry shape {lens_ast.shape}")
    if lens_ast.shape[0] != ast_xy.shape[0]:
        raise RuntimeError(
            f"unexpected BAGLE astrometry length mismatch (source={ast_xy.shape[0]}, lens={lens_ast.shape[0]})"
        )
    return ast_xy - lens_ast[:, :2]


def _should_try_scipy_fallback(exc: Exception) -> bool:
    text = str(exc).lower()
    needles = (
        "multinest",
        "dlsym",
        "pymultinest",
        "symbol not found",
        "libmultinest",
    )
    return any(tok in text for tok in needles)


def _solve_with_scipy_least_squares(
    model_module: Any,
    *,
    ra_deg: float,
    dec_deg: float,
    param_names: Sequence[str],
    init_params: Mapping[str, float],
    bounds: Mapping[str, Tuple[float, float]],
    t_phot: np.ndarray,
    mag_obs: np.ndarray,
    mag_err: np.ndarray,
    t_ast: np.ndarray,
    x_ast_arcsec: np.ndarray,
    y_ast_arcsec: np.ndarray,
    x_ast_err_arcsec: np.ndarray,
    y_ast_err_arcsec: np.ndarray,
    lens_relative_astrometry: bool = False,
    max_nfev: int = 2500,
) -> Dict[str, float]:
    try:
        from scipy.optimize import least_squares
    except Exception as exc:  # pragma: no cover - optional dependency guard
        raise SmokeTestError(
            "BAGLE sanity fallback requires scipy.optimize.least_squares."
        ) from exc

    x0 = np.array([float(init_params[name]) for name in param_names], dtype=float)
    lb = np.empty_like(x0)
    ub = np.empty_like(x0)

    for ii, name in enumerate(param_names):
        if name in bounds:
            lo, hi = bounds[name]
        else:
            center = float(init_params[name])
            half = max(1.0, abs(center) * 10.0)
            lo, hi = center - half, center + half
        lb[ii] = lo
        ub[ii] = hi
        if ub[ii] <= lb[ii]:
            ub[ii] = lb[ii] + 1.0e-6

    # Ensure initial guess sits inside the bounded domain.
    eps = 1.0e-9
    x0 = np.minimum(np.maximum(x0, lb + eps), ub - eps)

    def _residual(vec: np.ndarray) -> np.ndarray:
        p = {name: float(vec[ii]) for ii, name in enumerate(param_names)}
        pspl = model_module.PSPL_PhotAstrom_Par_Param1(
            p["mL"],
            p["t0"],
            p["beta"],
            p["dL"],
            p["dL_dS"],
            p["xS0_E"],
            p["xS0_N"],
            p["muL_E"],
            p["muL_N"],
            p["muS_E"],
            p["muS_N"],
            [p["b_sff1"]],
            [p["mag_src1"]],
            raL=ra_deg,
            decL=dec_deg,
        )

        mag_model = np.asarray(pspl.get_photometry(t_phot), dtype=float)
        ast_model = _get_model_astrometry(
            pspl,
            t_ast,
            lens_relative_astrometry=lens_relative_astrometry,
        )

        res_phot = (mag_obs - mag_model) / mag_err
        res_ast_e = (x_ast_arcsec - ast_model[:, 0]) / x_ast_err_arcsec
        res_ast_n = (y_ast_arcsec - ast_model[:, 1]) / y_ast_err_arcsec
        return np.concatenate([res_phot, res_ast_e, res_ast_n])

    result = least_squares(
        _residual,
        x0,
        bounds=(lb, ub),
        method="trf",
        max_nfev=int(max_nfev),
    )
    if not result.success:
        raise SmokeTestError(
            f"BAGLE sanity fallback least-squares failed: {result.message}"
        )

    return {name: float(result.x[ii]) for ii, name in enumerate(param_names)}


def _plot_joint_fit_diagnostics(
    output_path: Path,
    t_phot: np.ndarray,
    mag_obs: np.ndarray,
    mag_err: np.ndarray,
    t_ast: np.ndarray,
    x_ast_arcsec: np.ndarray,
    y_ast_arcsec: np.ndarray,
    x_ast_err_arcsec: np.ndarray,
    y_ast_err_arcsec: np.ndarray,
    t_ast_true: np.ndarray | None,
    x_ast_true_arcsec: np.ndarray | None,
    y_ast_true_arcsec: np.ndarray | None,
    model_obj: Any,
    t0_ref: float,
    title: str,
    lens_relative_astrometry: bool = False,
    show_noisy_astrometry: bool = True,
    noisy_ast_alpha: float = 0.18,
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - optional dependency guard
        raise SmokeTestError(
            "BAGLE sanity: matplotlib is required to generate fit diagnostics plot"
        ) from exc

    t_plot = np.linspace(
        min(float(np.min(t_phot)), float(np.min(t_ast))),
        max(float(np.max(t_phot)), float(np.max(t_ast))),
        800,
    )
    mag_model = np.asarray(model_obj.get_photometry(t_plot), dtype=float)
    ast_model = _get_model_astrometry(
        model_obj,
        t_plot,
        lens_relative_astrometry=lens_relative_astrometry,
    )

    fig, axes = plt.subplots(3, 1, figsize=(11, 12), constrained_layout=True)

    dt_phot = t_phot - t0_ref
    dt_ast = t_ast - t0_ref
    dt_plot = t_plot - t0_ref

    axes[0].errorbar(
        dt_phot,
        mag_obs,
        yerr=mag_err,
        fmt=".",
        ms=3.0,
        alpha=0.45,
        color="tab:blue",
        label="Data",
        zorder=1,
    )
    axes[0].plot(dt_plot, mag_model, "-", lw=1.8, color="tab:red", label="BAGLE best fit", zorder=5)
    axes[0].invert_yaxis()
    axes[0].set_ylabel("Magnitude")
    axes[0].set_title(title)
    axes[0].legend(loc="best", fontsize=9)

    if show_noisy_astrometry:
        noisy_ast_alpha = max(0.0, min(1.0, float(noisy_ast_alpha)))
        axes[1].errorbar(
            dt_ast,
            x_ast_arcsec * MAS_PER_ARCSEC,
            yerr=x_ast_err_arcsec * MAS_PER_ARCSEC,
            fmt=".",
            ms=2.5,
            alpha=noisy_ast_alpha,
            color="tab:green",
            label="E data",
            zorder=1,
        )
        axes[1].errorbar(
            dt_ast,
            y_ast_arcsec * MAS_PER_ARCSEC,
            yerr=y_ast_err_arcsec * MAS_PER_ARCSEC,
            fmt=".",
            ms=2.5,
            alpha=noisy_ast_alpha,
            color="tab:purple",
            label="N data",
            zorder=1,
        )
    axes[1].plot(
        dt_plot,
        ast_model[:, 0] * MAS_PER_ARCSEC,
        "-",
        lw=1.8,
        color="tab:green",
        label="E model",
        zorder=5,
    )
    axes[1].plot(
        dt_plot,
        ast_model[:, 1] * MAS_PER_ARCSEC,
        "-",
        lw=1.8,
        color="tab:purple",
        label="N model",
        zorder=5,
    )
    if (
        t_ast_true is not None
        and x_ast_true_arcsec is not None
        and y_ast_true_arcsec is not None
        and len(t_ast_true) > 0
    ):
        dt_ast_true = t_ast_true - t0_ref
        axes[1].plot(
            dt_ast_true,
            x_ast_true_arcsec * MAS_PER_ARCSEC,
            ".",
            ms=2.0,
            alpha=0.35,
            color="tab:olive",
            label="E noiseless",
            zorder=6,
        )
        axes[1].plot(
            dt_ast_true,
            y_ast_true_arcsec * MAS_PER_ARCSEC,
            ".",
            ms=2.0,
            alpha=0.35,
            color="tab:brown",
            label="N noiseless",
            zorder=6,
        )
    axes[1].set_ylabel("Centroid (mas)")
    axes[1].legend(loc="best", fontsize=9)

    axes[2].plot(
        ast_model[:, 0] * MAS_PER_ARCSEC,
        ast_model[:, 1] * MAS_PER_ARCSEC,
        "-",
        lw=2.0,
        color="tab:red",
        label="BAGLE best fit",
        zorder=5,
    )
    if show_noisy_astrometry:
        noisy_ast_alpha = max(0.0, min(1.0, float(noisy_ast_alpha)))
        axes[2].errorbar(
            x_ast_arcsec * MAS_PER_ARCSEC,
            y_ast_arcsec * MAS_PER_ARCSEC,
            xerr=x_ast_err_arcsec * MAS_PER_ARCSEC,
            yerr=y_ast_err_arcsec * MAS_PER_ARCSEC,
            fmt=".",
            ms=2.2,
            alpha=noisy_ast_alpha,
            color="tab:blue",
            label="Data",
            zorder=1,
        )
    if x_ast_true_arcsec is not None and y_ast_true_arcsec is not None and len(x_ast_true_arcsec) > 0:
        axes[2].plot(
            x_ast_true_arcsec * MAS_PER_ARCSEC,
            y_ast_true_arcsec * MAS_PER_ARCSEC,
            ".",
            ms=2.0,
            alpha=0.35,
            color="tab:orange",
            label="Noiseless centroid",
            zorder=6,
        )
    axes[2].set_xlabel("E (mas)")
    axes[2].set_ylabel("N (mas)")
    axes[2].set_aspect("equal", adjustable="box")
    axes[2].legend(loc="best", fontsize=9)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def run_bagle_joint_fit_sanity(
    run_dir: Path,
    out_files: Sequence[Path],
    params: Mapping[str, str],
    *,
    chi2_max: float = 20.0,
    event_id: int | None = None,
    max_phot_points: int = 500,
    max_ast_points: int = 500,
    n_live_points: int = 200,
    fit_true_astrometry: bool = False,
    true_ast_err_mas: float = 0.01,
    fit_reduced_chi2_max: float = 3.0,
    mu_amp_frac_tol: float = 0.20,
    mu_dir_tol_deg: float = 10.0,
    piE_amp_frac_tol: float = 0.35,
    piE_dir_tol_deg: float = 15.0,
    true_ast_rms_mas_max: float = 0.50,
    true_ast_sigma_max: float = 1.5,
) -> BagleJointFitSummary:
    """Fit one single-source, single-lens-like event with BAGLE and validate vectors."""
    run_dir = run_dir.resolve()
    fit_dir = run_dir / "bagle_fit_sanity"
    fit_dir.mkdir(parents=True, exist_ok=True)
    if true_ast_err_mas <= 0.0:
        raise SmokeTestError(
            f"BAGLE sanity: true_ast_err_mas must be positive, got {true_ast_err_mas}."
        )

    # BAGLE uses on-disk caches during import/runtime. In sandboxed environments,
    # default cache paths may be read-only, so provide writable defaults.
    runtime_cache = fit_dir / "_runtime_cache"
    runtime_cache.mkdir(parents=True, exist_ok=True)
    if "PARALLAX_CACHE_DIR" not in os.environ:
        parallax_cache = runtime_cache / "parallax_cache"
        parallax_cache.mkdir(parents=True, exist_ok=True)
        os.environ["PARALLAX_CACHE_DIR"] = str(parallax_cache)
    if "MPLCONFIGDIR" not in os.environ:
        mpl_cache = runtime_cache / "mplconfig"
        mpl_cache.mkdir(parents=True, exist_ok=True)
        os.environ["MPLCONFIGDIR"] = str(mpl_cache)
    if "XDG_CACHE_HOME" not in os.environ:
        xdg_cache = runtime_cache / "xdg_cache"
        xdg_cache.mkdir(parents=True, exist_ok=True)
        os.environ["XDG_CACHE_HOME"] = str(xdg_cache)

    try:
        from bagle import model
        from bagle import model_fitter
    except Exception as exc:  # pragma: no cover - external dependency
        raise SmokeTestError(
            "BAGLE sanity requires BAGLE and its fitting dependencies. "
            "Install BAGLE (and typically pymultinest, dynesty, ultranest) before running."
        ) from exc

    out_tables = [pd.read_csv(path, sep=r"\s+") for path in out_files]
    if not out_tables:
        raise SmokeTestError(f"BAGLE sanity: no .out files found under {run_dir}")
    out_df = pd.concat(out_tables, ignore_index=True)

    required_out_cols = [
        "EventID",
        "SubRun",
        "Field",
        "NSource",
        "ObsGroup_0_chi2",
        "u0lens1",
        "t0lens1",
        "tE_ref",
        "thetaE",
        "ra_deg",
        "dec_deg",
        "piEE",
        "piEN",
        "murel_helio_alpha",
        "murel_helio_delta",
    ]
    missing_out = [col for col in required_out_cols if col not in out_df.columns]
    if missing_out:
        raise SmokeTestError(
            f"BAGLE sanity: .out is missing required columns: {', '.join(missing_out)}"
        )

    sim_zero_raw = params.get("SIMULATION_ZERO_TIME")
    if sim_zero_raw is None:
        raise SmokeTestError(
            "BAGLE sanity requires SIMULATION_ZERO_TIME in the parameter file for absolute timing."
        )
    try:
        sim_zero_time = float(sim_zero_raw)
    except ValueError as exc:
        raise SmokeTestError(
            f"BAGLE sanity: invalid SIMULATION_ZERO_TIME value: {sim_zero_raw}"
        ) from exc

    finite_mask = np.isfinite(pd.to_numeric(out_df["ObsGroup_0_chi2"], errors="coerce"))
    select_mask = finite_mask & (pd.to_numeric(out_df["NSource"], errors="coerce") == 1.0)
    select_mask &= pd.to_numeric(out_df["ObsGroup_0_chi2"], errors="coerce") < chi2_max
    if "ObsGroup_0_FiniteSourceflag" in out_df.columns:
        select_mask &= pd.to_numeric(out_df["ObsGroup_0_FiniteSourceflag"], errors="coerce") == 0.0
    if event_id is not None:
        select_mask &= pd.to_numeric(out_df["EventID"], errors="coerce") == float(event_id)

    candidates = out_df.loc[select_mask].copy()
    if candidates.empty:
        shortlist = out_df.loc[finite_mask, ["EventID", "NSource", "ObsGroup_0_chi2"]].sort_values(
            "ObsGroup_0_chi2",
            ascending=True,
        )
        head = shortlist.head(8).to_string(index=False)
        raise SmokeTestError(
            "BAGLE sanity: no event satisfies the selection "
            f"(NSource==1, ObsGroup_0_chi2<{chi2_max}, finite-source flag==0).\n"
            f"Closest events:\n{head}"
        )

    row = candidates.sort_values(
        ["ObsGroup_0_chi2", "EventID"],
        ascending=[True, True],
    ).iloc[0]
    evt = int(float(row["EventID"]))
    subrun = int(float(row["SubRun"]))
    field = int(float(row["Field"]))
    out_chi2 = float(row["ObsGroup_0_chi2"])
    lensing_context_parts: List[str] = []
    if "NLens" in row.index:
        nlens_val = _safe_get(row, "NLens")
        if math.isfinite(nlens_val) and nlens_val > 1.0:
            lensing_context_parts.append(f"underlying event has NLens={int(round(nlens_val))}")
    if "NPlanets" in row.index:
        nplan_val = _safe_get(row, "NPlanets")
        if math.isfinite(nplan_val) and nplan_val > 0.0:
            lensing_context_parts.append(f"NPlanets={int(round(nplan_val))}")
    if "Planet_0_q" in row.index:
        q_val = _safe_get(row, "Planet_0_q")
        if math.isfinite(q_val) and q_val > 0.0:
            lensing_context_parts.append(f"Planet_0_q={q_val:.4g}")

    lc_file = _find_lc_for_event(run_dir, evt, subrun, field)
    df = pd.read_csv(lc_file, sep=r"\s+", comment="#")
    astrometry_model_frame_raw = _parse_astrometry_bagle_model_frame(lc_file)
    astrometry_model_frame = astrometry_model_frame_raw
    if astrometry_model_frame is None:
        astrometry_model_frame = "lens_relative"
    if astrometry_model_frame not in ("lens_relative", "absolute"):
        raise SmokeTestError(
            f"BAGLE sanity: unsupported #Astrometry_BAGLE model_frame={astrometry_model_frame!r} in {lc_file.name}."
        )
    lens_relative_astrometry = astrometry_model_frame == "lens_relative"
    if "lens0_x_thE" in df.columns and "lens0_y_thE" in df.columns and "lens1_x_thE" in df.columns and "lens1_y_thE" in df.columns:
        lens_sep_thE = np.hypot(
            df["lens1_x_thE"].to_numpy(dtype=float, copy=False) - df["lens0_x_thE"].to_numpy(dtype=float, copy=False),
            df["lens1_y_thE"].to_numpy(dtype=float, copy=False) - df["lens0_y_thE"].to_numpy(dtype=float, copy=False),
        )
        finite_sep = lens_sep_thE[np.isfinite(lens_sep_thE)]
        if finite_sep.size > 0:
            lensing_context_parts.append(
                "lens-separation_thE[min/med/max]={:.3f}/{:.3f}/{:.3f}".format(
                    float(np.min(finite_sep)),
                    float(np.median(finite_sep)),
                    float(np.max(finite_sep)),
                )
            )
    lensing_context = ""
    if lensing_context_parts:
        lensing_context = " Likely cause: " + "; ".join(lensing_context_parts) + "."

    required_lc_cols = [
        "Simulation_time",
        "measured_relative_flux",
        "measured_relative_flux_error",
        "RA_centroid_deg",
        "Dec_centroid_deg",
        "x_centroid_error_mas",
        "y_centroid_error_mas",
    ]
    missing_lc = [col for col in required_lc_cols if col not in df.columns]
    if missing_lc:
        raise SmokeTestError(
            f"BAGLE sanity: {lc_file.name} missing required columns: {', '.join(missing_lc)}"
        )

    sim_time = df["Simulation_time"].to_numpy(dtype=float, copy=False)
    warnings: List[str] = []
    if astrometry_model_frame_raw is None:
        warnings.append(
            f"{lc_file.name}: #Astrometry_BAGLE model_frame missing; assuming lens_relative per general output convention."
        )
    warnings.append(
        f"{lc_file.name}: BAGLE astrometry model_frame={astrometry_model_frame} (from #Astrometry_BAGLE contract)."
    )

    # Explicitly convert JD-like timing to MJD for BAGLE.
    # Prefer BJD when it appears sufficiently precise; otherwise use
    # SIMULATION_ZERO_TIME + Simulation_time (JD) and convert to MJD.
    t_jd_sim = sim_zero_time + sim_time
    t_jd = t_jd_sim.copy()
    if "BJD" in df.columns:
        bjd = df["BJD"].to_numpy(dtype=float, copy=False)
        finite_bjd = np.isfinite(bjd)
        n_bjd = int(np.sum(finite_bjd))
        if n_bjd > 0:
            unique_bjd = int(np.unique(np.round(bjd[finite_bjd], 8)).size)
            min_unique = max(50, n_bjd // 20)
            if unique_bjd >= min_unique:
                t_jd[finite_bjd] = bjd[finite_bjd]
            else:
                warnings.append(
                    f"{lc_file.name}: BJD appears quantized ({unique_bjd} unique values over {n_bjd} epochs); "
                    "using SIMULATION_ZERO_TIME + Simulation_time for high-precision timing."
                )
        else:
            warnings.append(
                f"{lc_file.name}: BJD column has no finite values; using SIMULATION_ZERO_TIME + Simulation_time."
            )
    else:
        warnings.append(
            f"{lc_file.name}: BJD column missing; using SIMULATION_ZERO_TIME + Simulation_time."
        )
    t_mjd = t_jd - 2400000.5

    flux = df["measured_relative_flux"].to_numpy(dtype=float, copy=False)
    flux_err = df["measured_relative_flux_error"].to_numpy(dtype=float, copy=False)
    phot_mask = (
        np.isfinite(t_mjd)
        & np.isfinite(flux)
        & np.isfinite(flux_err)
        & (flux > 0.0)
        & (flux_err > 0.0)
    )
    phot_idx = np.where(phot_mask)[0]
    if phot_idx.size < 60:
        raise SmokeTestError(
            f"BAGLE sanity: insufficient valid photometric epochs ({phot_idx.size}) in {lc_file.name}"
        )
    phot_idx = _uniform_sample(phot_idx, max_phot_points)
    t_phot = t_mjd[phot_idx]
    mag_zero_point = 20.0
    mag_obs = mag_zero_point - 2.5 * np.log10(flux[phot_idx])
    mag_err = LOG10_FACTOR * (flux_err[phot_idx] / flux[phot_idx])
    mag_err = np.clip(mag_err, 1.0e-4, None)

    ra_deg = _safe_get(row, "ra_deg")
    dec_deg = _safe_get(row, "dec_deg")
    if math.isnan(ra_deg) or math.isnan(dec_deg):
        ra_frame, dec_frame = _parse_astrometry_frame(lc_file)
        if ra_frame is None or dec_frame is None:
            raise SmokeTestError(
                f"BAGLE sanity: {lc_file.name} missing both .out and #Astrometry_Frame RA/Dec"
            )
        ra_deg = ra_frame
        dec_deg = dec_frame

    ra_obs_deg = df["RA_centroid_deg"].to_numpy(dtype=float, copy=False)
    dec_obs_deg = df["Dec_centroid_deg"].to_numpy(dtype=float, copy=False)
    cos_dec = math.cos(math.radians(dec_deg))
    if abs(cos_dec) < 1.0e-8:
        raise SmokeTestError(
            f"BAGLE sanity: cos(dec) too small at event pointing (dec={dec_deg:.8f} deg)."
        )
    dra_deg = (ra_obs_deg - ra_deg + 180.0) % 360.0 - 180.0
    x_ast_all_arcsec = dra_deg * cos_dec * 3600.0
    y_ast_all_arcsec = (dec_obs_deg - dec_deg) * 3600.0

    x_ast_true_all_arcsec: np.ndarray | None = None
    y_ast_true_all_arcsec: np.ndarray | None = None
    if "RA_centroid_true_deg" in df.columns and "Dec_centroid_true_deg" in df.columns:
        ra_true_deg = df["RA_centroid_true_deg"].to_numpy(dtype=float, copy=False)
        dec_true_deg = df["Dec_centroid_true_deg"].to_numpy(dtype=float, copy=False)
        dra_true_deg = (ra_true_deg - ra_deg + 180.0) % 360.0 - 180.0
        x_ast_true_all_arcsec = dra_true_deg * cos_dec * 3600.0
        y_ast_true_all_arcsec = (dec_true_deg - dec_deg) * 3600.0
    elif "true_x_centroid_mas" in df.columns and "true_y_centroid_mas" in df.columns:
        true_e_mas = df["true_x_centroid_mas"].to_numpy(dtype=float, copy=False)
        true_n_mas = df["true_y_centroid_mas"].to_numpy(dtype=float, copy=False)
        transform = _parse_astrometry_transform(lc_file)
        if transform is not None:
            a11, a12, a21, a22 = transform
            x_ast_true_all_arcsec = (a11 * true_e_mas + a12 * true_n_mas) / MAS_PER_ARCSEC
            y_ast_true_all_arcsec = (a21 * true_e_mas + a22 * true_n_mas) / MAS_PER_ARCSEC
            warnings.append(
                f"{lc_file.name}: derived noiseless BAGLE astrometry from true_x/y centroid columns using #Astrometry_Transform."
            )
        else:
            raise SmokeTestError(
                f"BAGLE sanity: {lc_file.name} has true_x/y centroid columns but no #Astrometry_Transform; "
                "cannot compare/plot noiseless astrometry in BAGLE frame."
            )
    else:
        raise SmokeTestError(
            f"BAGLE sanity: {lc_file.name} missing noiseless astrometry columns "
            "(need RA_centroid_true_deg/Dec_centroid_true_deg or true_x/y with #Astrometry_Transform)."
        )

    x_err_mas = df["x_centroid_error_mas"].to_numpy(dtype=float, copy=False)
    y_err_mas = df["y_centroid_error_mas"].to_numpy(dtype=float, copy=False)
    ast_mask = (
        np.isfinite(t_mjd)
        & np.isfinite(x_ast_all_arcsec)
        & np.isfinite(y_ast_all_arcsec)
        & np.isfinite(x_err_mas)
        & np.isfinite(y_err_mas)
        & (x_err_mas > 0.0)
        & (y_err_mas > 0.0)
    )
    ast_idx = np.where(ast_mask)[0]
    if ast_idx.size < 60:
        raise SmokeTestError(
            f"BAGLE sanity: insufficient valid astrometric epochs ({ast_idx.size}) in {lc_file.name}"
        )
    ast_idx = _uniform_sample(ast_idx, max_ast_points)
    t_ast = t_mjd[ast_idx]
    x_ast_arcsec = x_ast_all_arcsec[ast_idx]
    y_ast_arcsec = y_ast_all_arcsec[ast_idx]
    x_ast_err_arcsec = np.clip(x_err_mas[ast_idx] / MAS_PER_ARCSEC, 1.0e-6, None)
    y_ast_err_arcsec = np.clip(y_err_mas[ast_idx] / MAS_PER_ARCSEC, 1.0e-6, None)
    t_ast_true: np.ndarray | None = None
    x_ast_true_arcsec: np.ndarray | None = None
    y_ast_true_arcsec: np.ndarray | None = None
    true_idx_in_ast: np.ndarray | None = None
    if x_ast_true_all_arcsec is not None and y_ast_true_all_arcsec is not None:
        x_true_sel = x_ast_true_all_arcsec[ast_idx]
        y_true_sel = y_ast_true_all_arcsec[ast_idx]
        finite_true = np.isfinite(t_ast) & np.isfinite(x_true_sel) & np.isfinite(y_true_sel)
        if np.any(finite_true):
            true_idx_in_ast = np.where(finite_true)[0]
            t_ast_true = t_ast[finite_true]
            x_ast_true_arcsec = x_true_sel[finite_true]
            y_ast_true_arcsec = y_true_sel[finite_true]
        else:
            raise SmokeTestError(
                f"BAGLE sanity: {lc_file.name} noiseless astrometry columns are present but contain no finite values "
                "on selected astrometric epochs."
            )

    fit_t_ast = t_ast
    fit_x_ast_arcsec = x_ast_arcsec
    fit_y_ast_arcsec = y_ast_arcsec
    fit_x_ast_err_arcsec = x_ast_err_arcsec
    fit_y_ast_err_arcsec = y_ast_err_arcsec
    if fit_true_astrometry:
        if t_ast_true is None or x_ast_true_arcsec is None or y_ast_true_arcsec is None or len(t_ast_true) == 0:
            raise SmokeTestError(
                f"BAGLE sanity: fit_true_astrometry requested but no finite noiseless astrometry available for {lc_file.name}."
            )
        tiny_err_arcsec = max(true_ast_err_mas / MAS_PER_ARCSEC, 1.0e-9)
        fit_t_ast = t_ast_true
        fit_x_ast_arcsec = x_ast_true_arcsec
        fit_y_ast_arcsec = y_ast_true_arcsec
        fit_x_ast_err_arcsec = np.full_like(fit_x_ast_arcsec, tiny_err_arcsec, dtype=float)
        fit_y_ast_err_arcsec = np.full_like(fit_y_ast_arcsec, tiny_err_arcsec, dtype=float)
        warnings.append(
            f"{lc_file.name}: fitting BAGLE to noiseless astrometry with fixed uncertainty {true_ast_err_mas:.4g} mas."
        )

    t0_guess_jd = sim_zero_time + _safe_get(row, "t0lens1")
    t0_guess = t0_guess_jd - 2400000.5
    tE_guess = abs(_safe_get(row, "tE_ref"))
    thetaE_guess = abs(_safe_get(row, "thetaE"))
    if math.isnan(t0_guess) or math.isnan(tE_guess) or tE_guess <= 0.0:
        raise SmokeTestError(
            f"BAGLE sanity: event {evt} has invalid t0lens1/tE_ref in .out"
        )
    if math.isnan(thetaE_guess) or thetaE_guess <= 0.0:
        raise SmokeTestError(
            f"BAGLE sanity: event {evt} has invalid thetaE in .out"
        )
    beta_guess = _safe_get(row, "u0lens1") * thetaE_guess

    near_t0_idx = int(np.argmin(np.abs(t_mjd - t0_guess)))
    xS0_guess = float(x_ast_all_arcsec[near_t0_idx])
    yS0_guess = float(y_ast_all_arcsec[near_t0_idx])

    fs_guess = _safe_get(row, "Obs_0_fs")
    if not math.isfinite(fs_guess):
        fs_guess = 0.9
    fs_guess = min(max(fs_guess, 0.02), 0.999)

    far_mask = np.abs(sim_time - _safe_get(row, "t0lens1")) > 3.0 * max(tE_guess, 1.0)
    base_flux_pool = flux[(far_mask) & np.isfinite(flux) & (flux > 0.0)]
    if base_flux_pool.size == 0:
        base_flux_pool = flux[np.isfinite(flux) & (flux > 0.0)]
    base_flux = float(np.median(base_flux_pool))
    blend_mag_guess = mag_zero_point - 2.5 * math.log10(base_flux)
    mag_src_guess = blend_mag_guess - 2.5 * math.log10(fs_guess)

    mL_guess = _safe_get(row, "Lens_Mass")
    if not math.isfinite(mL_guess) or mL_guess <= 0.0:
        mL_guess = 0.3

    dL_kpc = _safe_get(row, "Lens_Dist")
    dS_kpc = _safe_get(row, "Source_Dist")
    if not math.isfinite(dL_kpc) or dL_kpc <= 0.0:
        raise SmokeTestError(f"BAGLE sanity: invalid Lens_Dist for event {evt}")
    if not math.isfinite(dS_kpc) or dS_kpc <= 0.0:
        raise SmokeTestError(f"BAGLE sanity: invalid Source_Dist for event {evt}")
    dL_pc = dL_kpc * 1000.0
    dS_pc = dS_kpc * 1000.0
    if dS_pc <= dL_pc:
        dS_pc = dL_pc * 1.05
        warnings.append(
            f"{lc_file.name}: Source_Dist<=Lens_Dist in .out; clamped dS to 1.05*dL for BAGLE prior setup."
        )
    dL_dS_guess = dL_pc / dS_pc

    mu_rel_e_ref = _safe_get(row, "murel_helio_alpha")
    mu_rel_n_ref = _safe_get(row, "murel_helio_delta")
    if math.isnan(mu_rel_e_ref) or math.isnan(mu_rel_n_ref):
        raise SmokeTestError(
            f"BAGLE sanity: event {evt} missing murel_helio_alpha/delta in .out"
        )

    muS_e = 0.0
    muS_n = 0.0
    muL_e = mu_rel_e_ref
    muL_n = mu_rel_n_ref
    if _HAS_ASTROPY:
        try:
            src_l = _safe_get(row, "Source_l")
            src_b = _safe_get(row, "Source_b")
            src_mul = _safe_get(row, "Source_mul")
            src_mub = _safe_get(row, "Source_mub")
            lens_l = _safe_get(row, "Lens_l")
            lens_b = _safe_get(row, "Lens_b")
            lens_mul = _safe_get(row, "Lens_mul")
            lens_mub = _safe_get(row, "Lens_mub")
            if all(
                math.isfinite(v)
                for v in [src_l, src_b, src_mul, src_mub, lens_l, lens_b, lens_mul, lens_mub]
            ):
                muS_e, muS_n = _pm_gal_to_icrs(src_l, src_b, src_mul, src_mub)
                muL_e, muL_n = _pm_gal_to_icrs(lens_l, lens_b, lens_mul, lens_mub)
        except Exception:
            warnings.append(
                f"{lc_file.name}: failed galactic->ICRS PM conversion; using murel_helio for BAGLE initialization."
            )
    else:
        warnings.append(
            f"{lc_file.name}: astropy not installed; using murel_helio for BAGLE initialization."
        )

    data = {
        "target": f"gulls_event_{evt}",
        "phot_data": ["sim1"],
        "ast_data": ["sim1"],
        "phot_files": [str(lc_file)],
        "ast_files": [str(lc_file)],
        "t_phot1": t_phot,
        "mag1": mag_obs,
        "mag_err1": mag_err,
        "t_ast1": fit_t_ast,
        "xpos1": fit_x_ast_arcsec,
        "ypos1": fit_y_ast_arcsec,
        "xpos_err1": fit_x_ast_err_arcsec,
        "ypos_err1": fit_y_ast_err_arcsec,
    }

    basename = fit_dir / f"event_{evt:06d}_"

    fitter = model_fitter.MicrolensSolver(
        data,
        model.PSPL_PhotAstrom_Par_Param1,
        n_live_points=int(n_live_points),
        outputfiles_basename=str(basename),
        resume=False,
    )

    fit_bounds: Dict[str, Tuple[float, float]] = {}

    def _set_prior_and_bounds(name: str, lo: float, hi: float) -> None:
        if not math.isfinite(lo) or not math.isfinite(hi) or hi <= lo:
            return
        fit_bounds[name] = (lo, hi)
        _set_uniform_prior(fitter, model_fitter, name, lo, hi)

    _set_prior_and_bounds("mL", max(0.01, 0.25 * mL_guess), max(0.08, 4.0 * mL_guess))
    half_t0 = max(10.0, 2.5 * tE_guess)
    _set_prior_and_bounds("t0", t0_guess - half_t0, t0_guess + half_t0)
    _set_prior_and_bounds("xS0_E", xS0_guess - 0.02, xS0_guess + 0.02)
    _set_prior_and_bounds("xS0_N", yS0_guess - 0.02, yS0_guess + 0.02)
    half_beta = max(0.5, 3.0 * abs(beta_guess) + 0.2)
    _set_prior_and_bounds("beta", beta_guess - half_beta, beta_guess + half_beta)
    mu_half_window = 15.0 if lens_relative_astrometry else 4.0
    _set_prior_and_bounds("muL_E", muL_e - mu_half_window, muL_e + mu_half_window)
    _set_prior_and_bounds("muL_N", muL_n - mu_half_window, muL_n + mu_half_window)
    _set_prior_and_bounds("muS_E", muS_e - mu_half_window, muS_e + mu_half_window)
    _set_prior_and_bounds("muS_N", muS_n - mu_half_window, muS_n + mu_half_window)
    half_dL = max(200.0, 0.25 * dL_pc)
    _set_prior_and_bounds("dL", dL_pc - half_dL, dL_pc + half_dL)
    _set_prior_and_bounds("dL_dS", max(0.02, dL_dS_guess - 0.12), min(0.98, dL_dS_guess + 0.12))
    _set_prior_and_bounds("b_sff1", max(0.02, fs_guess - 0.18), min(0.999, fs_guess + 0.18))
    _set_prior_and_bounds("mag_src1", mag_src_guess - 2.5, mag_src_guess + 2.5)

    needed = [
        "mL",
        "t0",
        "beta",
        "dL",
        "dL_dS",
        "xS0_E",
        "xS0_N",
        "muL_E",
        "muL_N",
        "muS_E",
        "muS_N",
        "b_sff1",
        "mag_src1",
    ]
    init_params: Dict[str, float] = {
        "mL": mL_guess,
        "t0": t0_guess,
        "beta": beta_guess,
        "dL": dL_pc,
        "dL_dS": dL_dS_guess,
        "xS0_E": xS0_guess,
        "xS0_N": yS0_guess,
        "muL_E": muL_e,
        "muL_N": muL_n,
        "muS_E": muS_e,
        "muS_N": muS_n,
        "b_sff1": fs_guess,
        "mag_src1": mag_src_guess,
    }

    use_scipy_solver = lens_relative_astrometry
    if use_scipy_solver:
        warnings.append(
            f"{lc_file.name}: using scipy BAGLE fallback because lens-relative astrometry is required "
            "for this event and MicrolensSolver assumes absolute astrometry."
        )
    else:
        try:
            fitter.solve()
            best = fitter.get_best_fit()
        except Exception as exc:
            if not _should_try_scipy_fallback(exc):
                raise SmokeTestError(
                    f"BAGLE sanity: joint fit failed for {lc_file.name}: {exc}"
                ) from exc
            warnings.append(
                f"{lc_file.name}: PyMultiNest solve unavailable ({exc}); used scipy least-squares BAGLE fallback."
            )
            use_scipy_solver = True

    if use_scipy_solver:
        best = _solve_with_scipy_least_squares(
            model,
            ra_deg=ra_deg,
            dec_deg=dec_deg,
            param_names=needed,
            init_params=init_params,
            bounds=fit_bounds,
            t_phot=t_phot,
            mag_obs=mag_obs,
            mag_err=mag_err,
            t_ast=fit_t_ast,
            x_ast_arcsec=fit_x_ast_arcsec,
            y_ast_arcsec=fit_y_ast_arcsec,
            x_ast_err_arcsec=fit_x_ast_err_arcsec,
            y_ast_err_arcsec=fit_y_ast_err_arcsec,
            lens_relative_astrometry=lens_relative_astrometry,
        )

    missing_best = [name for name in needed if name not in best]
    if missing_best:
        raise SmokeTestError(
            "BAGLE sanity: best-fit dictionary missing keys: "
            f"{', '.join(missing_best)}"
        )

    best_model = model.PSPL_PhotAstrom_Par_Param1(
        best["mL"],
        best["t0"],
        best["beta"],
        best["dL"],
        best["dL_dS"],
        best["xS0_E"],
        best["xS0_N"],
        best["muL_E"],
        best["muL_N"],
        best["muS_E"],
        best["muS_N"],
        [best["b_sff1"]],
        [best["mag_src1"]],
        raL=ra_deg,
        decL=dec_deg,
    )

    mag_model = np.asarray(best_model.get_photometry(t_phot), dtype=float)
    ast_model = _get_model_astrometry(
        best_model,
        fit_t_ast,
        lens_relative_astrometry=lens_relative_astrometry,
    )

    chi2_phot = float(np.sum(((mag_obs - mag_model) / mag_err) ** 2))
    chi2_ast = float(
        np.sum(((fit_x_ast_arcsec - ast_model[:, 0]) / fit_x_ast_err_arcsec) ** 2)
        + np.sum(((fit_y_ast_arcsec - ast_model[:, 1]) / fit_y_ast_err_arcsec) ** 2)
    )
    chi2_total = chi2_phot + chi2_ast
    n_param_eff = 13
    dof = max(1, int(len(mag_obs) + 2 * len(fit_t_ast) - n_param_eff))
    red_chi2 = chi2_total / dof
    if red_chi2 > fit_reduced_chi2_max:
        raise SmokeTestError(
            f"BAGLE sanity: poor joint-fit quality for {lc_file.name} "
            f"(reduced chi2={red_chi2:.3f}, threshold={fit_reduced_chi2_max:.3f}, "
            f"chi2_phot={chi2_phot:.2f}, chi2_ast={chi2_ast:.2f}, dof={dof})."
            f"{lensing_context}"
        )

    if (
        true_idx_in_ast is None
        or x_ast_true_arcsec is None
        or y_ast_true_arcsec is None
        or len(true_idx_in_ast) == 0
    ):
        raise SmokeTestError(
            "BAGLE sanity: noiseless astrometry was not carried through fitting; cannot perform strict model-vs-true checks."
        )
    if len(true_idx_in_ast) < 20:
        raise SmokeTestError(
            f"BAGLE sanity: only {len(true_idx_in_ast)} noiseless astrometric epochs available for strict validation."
        )

    ast_model_true = _get_model_astrometry(
        best_model,
        t_ast_true,
        lens_relative_astrometry=lens_relative_astrometry,
    )
    if len(ast_model_true) != len(t_ast_true):
        raise SmokeTestError(
            f"BAGLE sanity: unexpected astrometry length for noiseless epochs: {ast_model_true.shape}"
        )
    resid_true_e_mas = (ast_model_true[:, 0] - x_ast_true_arcsec) * MAS_PER_ARCSEC
    resid_true_n_mas = (ast_model_true[:, 1] - y_ast_true_arcsec) * MAS_PER_ARCSEC
    true_ast_rms_mas = float(np.sqrt(np.mean(resid_true_e_mas**2 + resid_true_n_mas**2)))

    pair_err_mas = np.hypot(
        x_ast_err_arcsec[true_idx_in_ast] * MAS_PER_ARCSEC,
        y_ast_err_arcsec[true_idx_in_ast] * MAS_PER_ARCSEC,
    )
    finite_pair_err = pair_err_mas[np.isfinite(pair_err_mas) & (pair_err_mas > 0.0)]
    if finite_pair_err.size == 0:
        raise SmokeTestError(
            "BAGLE sanity: invalid astrometric uncertainties when evaluating model-vs-noiseless residuals."
        )
    true_ast_sigma_equiv = true_ast_rms_mas / float(np.median(finite_pair_err))
    if true_ast_rms_mas > true_ast_rms_mas_max or true_ast_sigma_equiv > true_ast_sigma_max:
        raise SmokeTestError(
            f"BAGLE sanity: model does not track noiseless astrometry for {lc_file.name} "
            f"(RMS_true={true_ast_rms_mas:.3f} mas, limit={true_ast_rms_mas_max:.3f} mas; "
            f"sigma_equiv={true_ast_sigma_equiv:.2f}, limit={true_ast_sigma_max:.2f}; "
            f"epochs={len(true_idx_in_ast)})."
            f"{lensing_context}"
        )

    mu_fit = _vec_from_model_attr(best_model, "muRel")
    if mu_fit is None:
        raise SmokeTestError(
            "BAGLE sanity: best-fit model has no finite muRel vector; cannot compare proper motion."
        )
    # Convention handling:
    # When astrometry is lens-relative, compare against raw .out murel_helio (lens-source).
    # Otherwise compare with source-lens convention used by BAGLE muRel.
    mu_ref_raw = np.array([mu_rel_e_ref, mu_rel_n_ref], dtype=float)
    if lens_relative_astrometry:
        mu_ref = mu_ref_raw.copy()
        warnings.append(
            "PM comparison uses raw .out murel_helio_* (lens-source) because astrometry is lens-relative."
        )
    else:
        mu_ref = -mu_ref_raw
        warnings.append(
            "Applied sign conversion for PM comparison: .out murel_helio_* (lens-source) -> BAGLE muRel convention (source-lens)."
        )
    mu_amp_fit = float(np.hypot(mu_fit[0], mu_fit[1]))
    mu_amp_ref = float(np.hypot(mu_ref[0], mu_ref[1]))
    if mu_amp_ref <= 1.0e-6:
        raise SmokeTestError("BAGLE sanity: .out murel_helio vector is ~0, cannot compare direction.")
    mu_amp_frac = abs(mu_amp_fit - mu_amp_ref) / mu_amp_ref
    mu_dir_fit = math.degrees(math.atan2(mu_fit[1], mu_fit[0]))
    mu_dir_ref = math.degrees(math.atan2(mu_ref[1], mu_ref[0]))
    mu_dir_diff = abs(_angle_diff_deg(mu_dir_fit, mu_dir_ref))
    if mu_amp_frac > mu_amp_frac_tol or mu_dir_diff > mu_dir_tol_deg:
        raise SmokeTestError(
            f"BAGLE sanity: proper-motion mismatch for {lc_file.name} "
            f"(fit_mu=({mu_fit[0]:.4f},{mu_fit[1]:.4f}) mas/yr, "
            f"out_mu_bagle=({mu_ref[0]:.4f},{mu_ref[1]:.4f}) mas/yr, "
            f"out_mu_raw=({mu_ref_raw[0]:.4f},{mu_ref_raw[1]:.4f}) mas/yr, "
            f"|Δamp|/amp={mu_amp_frac:.3f} (tol={mu_amp_frac_tol:.3f}), "
            f"Δdir={mu_dir_diff:.2f} deg (tol={mu_dir_tol_deg:.2f} deg))."
            f"{lensing_context}"
        )

    piE_fit = _vec_from_model_attr(best_model, "piE")
    if piE_fit is None:
        raise SmokeTestError(
            "BAGLE sanity: best-fit model has no finite piE vector; cannot compare parallax."
        )
    # Same convention handling for piE (direction tied to mu_rel definition).
    piE_ref_raw = np.array([_safe_get(row, "piEE"), _safe_get(row, "piEN")], dtype=float)
    if lens_relative_astrometry:
        piE_ref = piE_ref_raw.copy()
        warnings.append(
            "Parallax comparison uses raw .out piEE/piEN (lens-source) because astrometry is lens-relative."
        )
    else:
        piE_ref = -piE_ref_raw
        warnings.append(
            "Applied sign conversion for parallax comparison: .out piEE/piEN (lens-source convention) -> BAGLE piE convention."
        )
    if not np.all(np.isfinite(piE_ref_raw)):
        raise SmokeTestError(
            "BAGLE sanity: .out missing finite piEE/piEN; cannot compare parallax."
        )
    piE_amp_fit = float(np.hypot(piE_fit[0], piE_fit[1]))
    piE_amp_ref = float(np.hypot(piE_ref[0], piE_ref[1]))
    if piE_amp_ref <= 1.0e-6:
        raise SmokeTestError("BAGLE sanity: .out parallax vector is ~0, cannot compare direction.")
    piE_amp_frac = abs(piE_amp_fit - piE_amp_ref) / piE_amp_ref
    piE_dir_fit = math.degrees(math.atan2(piE_fit[1], piE_fit[0]))
    piE_dir_ref = math.degrees(math.atan2(piE_ref[1], piE_ref[0]))
    piE_dir_diff = abs(_angle_diff_deg(piE_dir_fit, piE_dir_ref))
    if piE_amp_frac > piE_amp_frac_tol or piE_dir_diff > piE_dir_tol_deg:
        raise SmokeTestError(
            f"BAGLE sanity: parallax mismatch for {lc_file.name} "
            f"(fit_piE=({piE_fit[0]:.4f},{piE_fit[1]:.4f}), "
            f"out_piE_bagle=({piE_ref[0]:.4f},{piE_ref[1]:.4f}), "
            f"out_piE_raw=({piE_ref_raw[0]:.4f},{piE_ref_raw[1]:.4f}), "
            f"|Δamp|/amp={piE_amp_frac:.3f} (tol={piE_amp_frac_tol:.3f}), "
            f"Δdir={piE_dir_diff:.2f} deg (tol={piE_dir_tol_deg:.2f} deg))."
            f"{lensing_context}"
        )

    plot_path = fit_dir / f"event_{evt:06d}_best_model.png"
    _plot_joint_fit_diagnostics(
        plot_path,
        t_phot,
        mag_obs,
        mag_err,
        t_ast,
        x_ast_arcsec,
        y_ast_arcsec,
        x_ast_err_arcsec,
        y_ast_err_arcsec,
        t_ast_true,
        x_ast_true_arcsec,
        y_ast_true_arcsec,
        best_model,
        t0_guess,
        title=(
            f"BAGLE joint fit: event {evt} "
            f"(single-lens chi2={out_chi2:.3f}, reduced-fit-chi2={red_chi2:.3f})"
        ),
        lens_relative_astrometry=lens_relative_astrometry,
        show_noisy_astrometry=True,
        noisy_ast_alpha=(0.24 if fit_true_astrometry else 0.16),
    )

    result_json_path = fit_dir / f"event_{evt:06d}_summary.json"
    payload: Dict[str, Any] = {
        "event_id": evt,
        "subrun": subrun,
        "field": field,
        "lc_file": str(lc_file),
        "single_lens_chi2_out": out_chi2,
        "lensing_context": lensing_context_parts,
        "astrometry_fit_mode": ("noiseless" if fit_true_astrometry else "noisy"),
        "astrometry_model_frame": ("lens_relative" if lens_relative_astrometry else "absolute"),
        "astrometry_fit_error_mas": (
            float(true_ast_err_mas) if fit_true_astrometry else None
        ),
        "n_phot_points": int(len(t_phot)),
        "n_ast_points": int(len(fit_t_ast)),
        "n_ast_points_noisy": int(len(t_ast)),
        "n_ast_points_true": int(len(t_ast_true) if t_ast_true is not None else 0),
        "fit_reduced_chi2": red_chi2,
        "astrometry_true_check": {
            "rms_mas": true_ast_rms_mas,
            "sigma_equiv": true_ast_sigma_equiv,
            "rms_limit_mas": float(true_ast_rms_mas_max),
            "sigma_limit": float(true_ast_sigma_max),
            "n_epochs": int(len(true_idx_in_ast)),
        },
        "proper_motion": {
            "fit": {"E": float(mu_fit[0]), "N": float(mu_fit[1]), "amp": mu_amp_fit},
            "reference_bagle_convention": {"E": float(mu_ref[0]), "N": float(mu_ref[1]), "amp": mu_amp_ref},
            "reference_out_raw": {"E": float(mu_ref_raw[0]), "N": float(mu_ref_raw[1]), "amp": float(np.hypot(mu_ref_raw[0], mu_ref_raw[1]))},
            "direction_difference_deg": mu_dir_diff,
        },
        "parallax": {
            "fit": {"E": float(piE_fit[0]), "N": float(piE_fit[1]), "amp": piE_amp_fit},
            "reference_bagle_convention": {"E": float(piE_ref[0]), "N": float(piE_ref[1]), "amp": piE_amp_ref},
            "reference_out_raw": {"E": float(piE_ref_raw[0]), "N": float(piE_ref_raw[1]), "amp": float(np.hypot(piE_ref_raw[0], piE_ref_raw[1]))},
            "direction_difference_deg": piE_dir_diff,
        },
        "warnings": warnings,
    }
    result_json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    return BagleJointFitSummary(
        run_dir=run_dir,
        lc_file=lc_file,
        event_id=evt,
        subrun=subrun,
        field=field,
        out_chi2_single_lens=out_chi2,
        n_phot_points=int(len(t_phot)),
        n_ast_points=int(len(fit_t_ast)),
        fit_reduced_chi2=red_chi2,
        mu_ref_e=float(mu_ref[0]),
        mu_ref_n=float(mu_ref[1]),
        mu_fit_e=float(mu_fit[0]),
        mu_fit_n=float(mu_fit[1]),
        mu_amp_ref=mu_amp_ref,
        mu_amp_fit=mu_amp_fit,
        mu_dir_diff_deg=mu_dir_diff,
        piE_ref_e=float(piE_ref[0]),
        piE_ref_n=float(piE_ref[1]),
        piE_fit_e=float(piE_fit[0]),
        piE_fit_n=float(piE_fit[1]),
        piE_amp_ref=piE_amp_ref,
        piE_amp_fit=piE_amp_fit,
        piE_dir_diff_deg=piE_dir_diff,
        true_ast_rms_mas=true_ast_rms_mas,
        true_ast_sigma_equiv=true_ast_sigma_equiv,
        plot_path=plot_path,
        result_json_path=result_json_path,
        warnings=warnings,
    )


__all__ = ["BagleJointFitSummary", "run_bagle_joint_fit_sanity"]
