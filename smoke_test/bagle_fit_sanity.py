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
    obs_location_used: str
    lens_ast_rms_raw_mas: float | None
    lens_ast_rms_demean_mas: float | None
    plot_path: Path
    lens_plot_path: Path | None
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


def _parse_header_keyvals(lc_file: Path, prefix: str) -> Dict[str, str]:
    with lc_file.open(encoding="utf-8") as handle:
        for raw in handle:
            if not raw.startswith("#"):
                break
            if not raw.startswith(prefix):
                continue
            parsed: Dict[str, str] = {}
            for token in raw.strip().split()[1:]:
                if "=" not in token:
                    continue
                key, value = token.split("=", 1)
                parsed[key.strip()] = value.strip()
            return parsed
    return {}


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


def _build_pspl_model(
    model_module: Any,
    params: Mapping[str, float],
    *,
    ra_deg: float,
    dec_deg: float,
    obs_location: str,
) -> Any:
    return model_module.PSPL_PhotAstrom_Par_Param1(
        params["mL"],
        params["t0"],
        params["beta"],
        params["dL"],
        params["dL_dS"],
        params["xS0_E"],
        params["xS0_N"],
        params["muL_E"],
        params["muL_N"],
        params["muS_E"],
        params["muS_N"],
        [params["b_sff1"]],
        [params["mag_src1"]],
        raL=ra_deg,
        decL=dec_deg,
        obsLocation=obs_location,
    )


def _select_obs_location(
    model_module: Any,
    *,
    requested: str,
    init_params: Mapping[str, float],
    t_probe_mjd: float,
    ra_deg: float,
    dec_deg: float,
) -> str:
    requested_clean = requested.strip() if requested else ""
    candidate = requested_clean or "earth"

    probe_time = np.asarray([t_probe_mjd], dtype=float)
    try:
        probe_model = _build_pspl_model(
            model_module,
            init_params,
            ra_deg=ra_deg,
            dec_deg=dec_deg,
            obs_location=candidate,
        )
        _ = np.asarray(probe_model.get_astrometry(probe_time), dtype=float)
        return candidate
    except Exception as exc:
        raise SmokeTestError(
            f"BAGLE sanity: requested obsLocation={candidate!r} failed to initialize: {exc}"
        ) from exc


def _solve_with_scipy_least_squares(
    model_module: Any,
    *,
    ra_deg: float,
    dec_deg: float,
    obs_location: str,
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
        pspl = _build_pspl_model(
            model_module,
            p,
            ra_deg=ra_deg,
            dec_deg=dec_deg,
            obs_location=obs_location,
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
    axes[0].plot(
        dt_plot,
        mag_model,
        "--",
        lw=2.0,
        color="tab:red",
        label="BAGLE best fit",
        zorder=20,
    )
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
        "--",
        lw=2.0,
        color="tab:green",
        label="E model",
        zorder=20,
    )
    axes[1].plot(
        dt_plot,
        ast_model[:, 1] * MAS_PER_ARCSEC,
        "--",
        lw=2.0,
        color="tab:purple",
        label="N model",
        zorder=20,
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
        "--",
        lw=2.0,
        color="tab:red",
        label="BAGLE best fit",
        zorder=20,
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


def _plot_lens_track_diagnostics(
    output_path: Path,
    t_lens: np.ndarray,
    lens_obs_x_arcsec: np.ndarray,
    lens_obs_y_arcsec: np.ndarray,
    lens_formula_x_arcsec: np.ndarray | None,
    lens_formula_y_arcsec: np.ndarray | None,
    model_obj: Any,
    t0_ref: float,
    title: str,
    subtitle: str | None = None,
    lens_rms_raw_mas: float | None = None,
    lens_rms_demean_mas: float | None = None,
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - optional dependency guard
        raise SmokeTestError(
            "BAGLE sanity: matplotlib is required to generate lens-track diagnostics plot"
        ) from exc

    t_lens = np.asarray(t_lens, dtype=float)
    lens_obs_x_arcsec = np.asarray(lens_obs_x_arcsec, dtype=float)
    lens_obs_y_arcsec = np.asarray(lens_obs_y_arcsec, dtype=float)
    finite = (
        np.isfinite(t_lens)
        & np.isfinite(lens_obs_x_arcsec)
        & np.isfinite(lens_obs_y_arcsec)
    )
    if int(np.sum(finite)) < 2:
        raise SmokeTestError(
            "BAGLE sanity: insufficient finite points to plot lens-track diagnostics."
        )
    t_lens = t_lens[finite]
    lens_obs_x_arcsec = lens_obs_x_arcsec[finite]
    lens_obs_y_arcsec = lens_obs_y_arcsec[finite]

    t_dense = np.linspace(float(np.min(t_lens)), float(np.max(t_lens)), 1200)
    lens_model_dense = np.asarray(model_obj.get_lens_astrometry(t_dense), dtype=float)
    lens_model_obs = np.asarray(model_obj.get_lens_astrometry(t_lens), dtype=float)
    if (
        lens_model_dense.ndim != 2
        or lens_model_obs.ndim != 2
        or lens_model_dense.shape[1] < 2
        or lens_model_obs.shape[1] < 2
    ):
        raise SmokeTestError(
            "BAGLE sanity: unexpected get_lens_astrometry output shape while plotting lens tracks."
        )

    lens_model_dense_xy = lens_model_dense[:, :2]
    lens_model_obs_xy = lens_model_obs[:, :2]
    if lens_model_obs_xy.shape[0] != len(t_lens):
        raise SmokeTestError(
            "BAGLE sanity: BAGLE lens astrometry length mismatch while plotting lens tracks."
        )

    lens_formula_ok = (
        lens_formula_x_arcsec is not None
        and lens_formula_y_arcsec is not None
        and len(lens_formula_x_arcsec) == len(t_lens)
        and len(lens_formula_y_arcsec) == len(t_lens)
    )

    resid_e_mas = (lens_obs_x_arcsec - lens_model_obs_xy[:, 0]) * MAS_PER_ARCSEC
    resid_n_mas = (lens_obs_y_arcsec - lens_model_obs_xy[:, 1]) * MAS_PER_ARCSEC
    if lens_rms_raw_mas is None:
        lens_rms_raw_mas = float(np.sqrt(np.mean(resid_e_mas**2 + resid_n_mas**2)))
    if lens_rms_demean_mas is None:
        resid_e_dm = resid_e_mas - float(np.median(resid_e_mas))
        resid_n_dm = resid_n_mas - float(np.median(resid_n_mas))
        lens_rms_demean_mas = float(np.sqrt(np.mean(resid_e_dm**2 + resid_n_dm**2)))

    dt_lens = t_lens - t0_ref
    dt_dense = t_dense - t0_ref

    fig, axes = plt.subplots(2, 1, figsize=(11, 9), constrained_layout=True)

    axes[0].plot(
        dt_dense,
        lens_model_dense_xy[:, 0] * MAS_PER_ARCSEC,
        "-",
        lw=1.9,
        color="tab:red",
        label="E model",
        zorder=5,
    )
    axes[0].plot(
        dt_dense,
        lens_model_dense_xy[:, 1] * MAS_PER_ARCSEC,
        "-",
        lw=1.9,
        color="tab:orange",
        label="N model",
        zorder=5,
    )
    axes[0].plot(
        dt_lens,
        lens_obs_x_arcsec * MAS_PER_ARCSEC,
        ".",
        ms=2.2,
        alpha=0.50,
        color="tab:blue",
        label="E Gulls",
        zorder=1,
    )
    axes[0].plot(
        dt_lens,
        lens_obs_y_arcsec * MAS_PER_ARCSEC,
        ".",
        ms=2.2,
        alpha=0.50,
        color="tab:green",
        label="N Gulls",
        zorder=1,
    )
    if lens_formula_ok:
        axes[0].plot(
            dt_lens,
            np.asarray(lens_formula_x_arcsec) * MAS_PER_ARCSEC,
            "--",
            lw=1.4,
            color="0.35",
            label="E pm+pllx",
            zorder=4,
        )
        axes[0].plot(
            dt_lens,
            np.asarray(lens_formula_y_arcsec) * MAS_PER_ARCSEC,
            "--",
            lw=1.4,
            color="0.55",
            label="N pm+pllx",
            zorder=4,
        )
    axes[0].set_ylabel("Lens position (mas)")
    if subtitle:
        axes[0].set_title(f"{title}\n{subtitle}")
    else:
        axes[0].set_title(title)
    axes[0].legend(loc="best", fontsize=9)

    axes[1].plot(
        lens_model_dense_xy[:, 0] * MAS_PER_ARCSEC,
        lens_model_dense_xy[:, 1] * MAS_PER_ARCSEC,
        "-",
        lw=2.0,
        color="tab:red",
        label="BAGLE lens track",
        zorder=5,
    )
    axes[1].plot(
        lens_obs_x_arcsec * MAS_PER_ARCSEC,
        lens_obs_y_arcsec * MAS_PER_ARCSEC,
        ".",
        ms=2.0,
        alpha=0.45,
        color="tab:blue",
        label="Gulls lens columns",
        zorder=1,
    )
    if lens_formula_ok:
        axes[1].plot(
            np.asarray(lens_formula_x_arcsec) * MAS_PER_ARCSEC,
            np.asarray(lens_formula_y_arcsec) * MAS_PER_ARCSEC,
            "--",
            lw=1.6,
            color="0.35",
            label="Gulls lens pm+pllx",
            zorder=4,
        )
    axes[1].set_xlabel("E (mas)")
    axes[1].set_ylabel("N (mas)")
    axes[1].set_aspect("equal", adjustable="box")
    axes[1].legend(loc="best", fontsize=9)
    axes[1].text(
        0.02,
        0.98,
        (
            f"RMS raw={lens_rms_raw_mas:.3f} mas\n"
            f"RMS after XY offset={lens_rms_demean_mas:.3f} mas"
        ),
        transform=axes[1].transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85, "edgecolor": "0.6"},
    )

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
    max_phot_points: int = 0,
    max_ast_points: int = 0,
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
    obs_location: str = "earth",
    lens_ast_rms_demean_mas_max: float = 0.05,
    microlensing_mask_te: float = 0.0,
) -> BagleJointFitSummary:
    """Fit one single-source, single-lens-like event with BAGLE and validate vectors."""
    run_dir = run_dir.resolve()
    fit_dir = run_dir / "bagle_fit_sanity"
    fit_dir.mkdir(parents=True, exist_ok=True)
    if true_ast_err_mas <= 0.0:
        raise SmokeTestError(
            f"BAGLE sanity: true_ast_err_mas must be positive, got {true_ast_err_mas}."
        )
    if lens_ast_rms_demean_mas_max <= 0.0:
        raise SmokeTestError(
            f"BAGLE sanity: lens_ast_rms_demean_mas_max must be positive, got {lens_ast_rms_demean_mas_max}."
        )
    if microlensing_mask_te < 0.0:
        raise SmokeTestError(
            f"BAGLE sanity: microlensing_mask_te must be non-negative, got {microlensing_mask_te}."
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

    def _read_out_table_strict(path: Path) -> pd.DataFrame:
        lines = path.read_text(encoding="utf-8").splitlines()
        if not lines:
            raise SmokeTestError(f"BAGLE sanity: .out file is empty: {path}")

        header_cols = lines[0].split()
        expected_ncols = len(header_cols)
        if expected_ncols == 0:
            raise SmokeTestError(f"BAGLE sanity: .out header is empty: {path}")

        rows: List[List[str]] = []
        for line_number, raw in enumerate(lines[1:], start=2):
            stripped = raw.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if "\t" in raw:
                raise SmokeTestError(
                    f"BAGLE sanity: malformed .out table {path}:{line_number} contains tab delimiters"
                )
            row = raw.split()
            if len(row) != expected_ncols:
                raise SmokeTestError(
                    f"BAGLE sanity: malformed .out table {path}:{line_number} has {len(row)} columns; "
                    f"header has {expected_ncols}"
                )
            rows.append(row)

        if not rows:
            raise SmokeTestError(f"BAGLE sanity: .out file has no data rows: {path}")

        return pd.DataFrame(rows, columns=header_cols)

    out_tables = [_read_out_table_strict(path) for path in out_files]
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

    candidate_rows = candidates.sort_values(
        ["ObsGroup_0_chi2", "EventID"],
        ascending=[True, True],
    )

    row = None
    lc_file = None
    df = None
    mask_skips: List[str] = []
    min_epochs_after_mask = 60
    for _, candidate_row in candidate_rows.iterrows():
        evt_try = int(float(candidate_row["EventID"]))
        subrun_try = int(float(candidate_row["SubRun"]))
        field_try = int(float(candidate_row["Field"]))
        lc_try = _find_lc_for_event(run_dir, evt_try, subrun_try, field_try)
        df_try = pd.read_csv(lc_try, sep=r"\s+", comment="#")
        if "Simulation_time" not in df_try.columns:
            raise SmokeTestError(
                f"BAGLE sanity: {lc_try.name} missing required Simulation_time column"
            )
        sim_time_try = df_try["Simulation_time"].to_numpy(dtype=float, copy=False)

        t0_try = _safe_get(candidate_row, "t0lens1")
        tE_try = abs(_safe_get(candidate_row, "tE_ref"))
        if math.isnan(t0_try) or math.isnan(tE_try) or tE_try <= 0.0:
            continue

        if microlensing_mask_te > 0.0:
            mask_halfwidth_days_try = microlensing_mask_te * max(tE_try, 1.0e-12)
            outside_event_window_try = np.abs(sim_time_try - t0_try) > mask_halfwidth_days_try
        else:
            outside_event_window_try = np.ones_like(sim_time_try, dtype=bool)

        n_outside = int(np.sum(outside_event_window_try & np.isfinite(sim_time_try)))
        if n_outside < min_epochs_after_mask:
            mask_skips.append(
                f"{lc_try.name} (EventID={evt_try}): only {n_outside} epochs remain after mask"
            )
            continue

        row = candidate_row
        lc_file = lc_try
        df = df_try
        break

    if row is None or lc_file is None or df is None:
        detail = "\n".join(mask_skips[:8]) if mask_skips else "No candidate passed mask-based epoch availability checks."
        raise SmokeTestError(
            "BAGLE sanity: all selected events were skipped because the microlensing mask left insufficient data.\n"
            + detail
        )

    evt = int(float(row["EventID"]))
    subrun = int(float(row["SubRun"]))
    field = int(float(row["Field"]))
    out_chi2 = float(row["ObsGroup_0_chi2"])
    lensing_context_parts: List[str] = []

    bagle_contract = _parse_header_keyvals(lc_file, "#Astrometry_BAGLE:")
    contract = _parse_header_keyvals(lc_file, "#Astrometry_Contract:")
    contract_cols = _parse_header_keyvals(lc_file, "#Astrometry_Columns:")
    blendless_cols_raw = bagle_contract.get("blendless_columns", "").strip()
    warnings: List[str] = []
    deferred_failures: List[str] = []
    for skip_msg in mask_skips[:5]:
        warnings.append(f"BAGLE candidate skipped: {skip_msg}.")

    def _resolve_col(contract_keys: tuple[str, ...], *fallbacks: str, allow_none: bool = False) -> str | None:
        candidates: List[str] = []
        for contract_key in contract_keys:
            if contract_key not in contract_cols:
                continue
            mapped = contract_cols[contract_key]
            if mapped.strip().lower() == "none":
                return None if allow_none else None
            candidates.append(mapped)
        candidates.extend(fallbacks)
        for col in candidates:
            if col in df.columns:
                return col
        return None

    astrometry_model_frame_raw = bagle_contract.get("model_frame")
    if astrometry_model_frame_raw is not None:
        astrometry_model_frame_raw = astrometry_model_frame_raw.strip().lower()
    else:
        astrometry_model_frame_raw = _parse_astrometry_bagle_model_frame(lc_file)
    if astrometry_model_frame_raw is None:
        astrometry_model_frame_raw = contract.get("model_frame")
        if astrometry_model_frame_raw is None:
            astrometry_model_frame_raw = "absolute"
            warnings.append(
                f"{lc_file.name}: missing #Astrometry_BAGLE model_frame; defaulting BAGLE comparison to absolute frame."
            )
        else:
            warnings.append(
                f"{lc_file.name}: BAGLE model_frame inferred from #Astrometry_Contract ({astrometry_model_frame_raw})."
            )
    astrometry_model_frame = astrometry_model_frame_raw
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

    col_ra_obs = _resolve_col(
        ("sky_ra_measured_deg", "sky_ra_obs_deg"),
        "RA_measured_deg",
        "RA_obs_deg",
        "RA_centroid_deg",
    )
    col_dec_obs = _resolve_col(
        ("sky_dec_measured_deg", "sky_dec_obs_deg"),
        "Dec_measured_deg",
        "Dec_obs_deg",
        "Dec_centroid_deg",
    )
    col_ra_det = _resolve_col(
        ("sky_ra_noiseless_deg", "sky_ra_det_deg"),
        "RA_noiseless_deg",
        "RA_det_deg",
        "RA_centroid_true_deg",
    )
    col_dec_det = _resolve_col(
        ("sky_dec_noiseless_deg", "sky_dec_det_deg"),
        "Dec_noiseless_deg",
        "Dec_det_deg",
        "Dec_centroid_true_deg",
    )
    col_sigma = _resolve_col(("sky_sigma_mas",), "sigma_astrometric_mas")

    required_lc_cols = [
        "Simulation_time",
        "measured_relative_flux",
        "measured_relative_flux_error",
        col_ra_obs,
        col_dec_obs,
        col_sigma,
    ]
    missing_lc = [col for col in required_lc_cols if col is None or col not in df.columns]
    if missing_lc:
        raise SmokeTestError(
            f"BAGLE sanity: {lc_file.name} missing required columns: {', '.join(missing_lc)}"
        )

    sim_time = df["Simulation_time"].to_numpy(dtype=float, copy=False)
    warnings.append(
        f"{lc_file.name}: BAGLE astrometry model_frame={astrometry_model_frame} (from #Astrometry_BAGLE contract)."
    )

    # Explicitly convert the published simulation clock into JD/MJD for BAGLE.
    # BJD is treated as a cross-check diagnostic only; the BAGLE fit uses the
    # deterministic simulation timing.
    t_jd_sim = sim_zero_time + sim_time
    t_jd = t_jd_sim.copy()
    timing_bjd_diff_stats: Dict[str, float | int] | None = None
    timing_bjd_diff_sec_max = 1.0e-4
    if "BJD" in df.columns:
        bjd = df["BJD"].to_numpy(dtype=float, copy=False)
        finite_bjd = np.isfinite(bjd) & np.isfinite(t_jd_sim)
        n_bjd = int(np.sum(finite_bjd))
        if n_bjd > 0:
            diff_sec = (bjd[finite_bjd] - t_jd_sim[finite_bjd]) * 86400.0
            max_abs_diff_sec = float(np.max(np.abs(diff_sec)))
            timing_bjd_diff_stats = {
                "n_finite": n_bjd,
                "min_sec": float(np.min(diff_sec)),
                "max_sec": float(np.max(diff_sec)),
                "mean_sec": float(np.mean(diff_sec)),
                "std_sec": float(np.std(diff_sec)),
                "max_abs_sec": max_abs_diff_sec,
                "threshold_sec": float(timing_bjd_diff_sec_max),
            }
            print(
                f"{lc_file.name}: BJD - (SIMULATION_ZERO_TIME + Simulation_time) [sec] "
                f"min={timing_bjd_diff_stats['min_sec']:.9g} "
                f"max={timing_bjd_diff_stats['max_sec']:.9g} "
                f"mean={timing_bjd_diff_stats['mean_sec']:.9g} "
                f"std={timing_bjd_diff_stats['std_sec']:.9g}"
            )
            finite_idx = np.flatnonzero(finite_bjd)
            for sample_idx in finite_idx[:5]:
                diff_i_sec = (bjd[sample_idx] - t_jd_sim[sample_idx]) * 86400.0
                print(
                    f"  ep={sample_idx} sim_jd={t_jd_sim[sample_idx]:.16f} "
                    f"BJD={bjd[sample_idx]:.16f} diff_sec={diff_i_sec:.9g}"
                )
            if max_abs_diff_sec > timing_bjd_diff_sec_max:
                deferred_failures.append(
                    f"BAGLE sanity: BJD differs from (SIMULATION_ZERO_TIME + Simulation_time) "
                    f"by up to {max_abs_diff_sec:.9g} s in {lc_file.name}, "
                    f"exceeding the {timing_bjd_diff_sec_max:.9g} s limit."
                )
        else:
            warnings.append(f"{lc_file.name}: BJD column is present but contains no finite values.")
    else:
        warnings.append(f"{lc_file.name}: BJD column missing; skipped timing cross-check.")
    t_mjd = t_jd - 2400000.5

    t0_event_days = _safe_get(row, "t0lens1")
    tE_event_days = abs(_safe_get(row, "tE_ref"))
    if math.isnan(t0_event_days) or math.isnan(tE_event_days) or tE_event_days <= 0.0:
        raise SmokeTestError(
            f"BAGLE sanity: event {evt} has invalid t0lens1/tE_ref in .out"
        )
    if microlensing_mask_te > 0.0:
        mask_halfwidth_days = microlensing_mask_te * max(tE_event_days, 1.0e-12)
        outside_event_window = np.abs(sim_time - t0_event_days) > mask_halfwidth_days
        warnings.append(
            f"{lc_file.name}: BAGLE fit/validation mask applied: excluded epochs with "
            f"|t-t0|<={microlensing_mask_te:.2f}*tE (|t-t0|<={mask_halfwidth_days:.3f} day)."
        )
    else:
        outside_event_window = np.ones_like(sim_time, dtype=bool)
        warnings.append(
            f"{lc_file.name}: BAGLE fit/validation mask disabled (microlensing_mask_te=0)."
        )

    flux = df["measured_relative_flux"].to_numpy(dtype=float, copy=False)
    flux_err = df["measured_relative_flux_error"].to_numpy(dtype=float, copy=False)
    phot_mask = (
        outside_event_window
        &
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

    ra_obs_deg = df[col_ra_obs].to_numpy(dtype=float, copy=False)
    dec_obs_deg = df[col_dec_obs].to_numpy(dtype=float, copy=False)
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
    true_ast_source: str | None = None
    if blendless_cols_raw and blendless_cols_raw.lower() != "none":
        split_blendless = [part.strip() for part in blendless_cols_raw.split(",") if part.strip()]
        if len(split_blendless) != 2:
            raise SmokeTestError(
                f"BAGLE sanity: invalid #Astrometry_BAGLE blendless_columns={blendless_cols_raw!r} "
                f"in {lc_file.name}; expected two comma-separated columns."
            )
        blendless_x_col, blendless_y_col = split_blendless
        if blendless_x_col not in df.columns or blendless_y_col not in df.columns:
            raise SmokeTestError(
                f"BAGLE sanity: declared blendless astrometry columns "
                f"({blendless_x_col}/{blendless_y_col}) are missing in {lc_file.name}."
            )

        if blendless_x_col == "RA_centroid_src_only_deg" and blendless_y_col == "Dec_centroid_src_only_deg":
            ra_src_only_deg = df[blendless_x_col].to_numpy(dtype=float, copy=False)
            dec_src_only_deg = df[blendless_y_col].to_numpy(dtype=float, copy=False)
            dra_src_only_deg = (ra_src_only_deg - ra_deg + 180.0) % 360.0 - 180.0
            x_ast_true_all_arcsec = dra_src_only_deg * cos_dec * 3600.0
            y_ast_true_all_arcsec = (dec_src_only_deg - dec_deg) * 3600.0
            true_ast_source = "contract blendless RA/Dec columns"
        elif blendless_x_col == "centroid_src_x_mas" and blendless_y_col == "centroid_src_y_mas":
            src_only_e_mas = df[blendless_x_col].to_numpy(dtype=float, copy=False)
            src_only_n_mas = df[blendless_y_col].to_numpy(dtype=float, copy=False)
            transform = _parse_astrometry_transform(lc_file)
            if transform is None:
                raise SmokeTestError(
                    f"BAGLE sanity: {lc_file.name} declares blendless centroid columns "
                    f"({blendless_x_col}/{blendless_y_col}) but no #Astrometry_Transform is available."
                )
            a11, a12, a21, a22 = transform
            x_ast_true_all_arcsec = (a11 * src_only_e_mas + a12 * src_only_n_mas) / MAS_PER_ARCSEC
            y_ast_true_all_arcsec = (a21 * src_only_e_mas + a22 * src_only_n_mas) / MAS_PER_ARCSEC
            true_ast_source = "contract blendless centroid columns"
        else:
            raise SmokeTestError(
                f"BAGLE sanity: unsupported #Astrometry_BAGLE blendless_columns "
                f"({blendless_x_col}/{blendless_y_col}) in {lc_file.name}."
            )
        warnings.append(
            f"{lc_file.name}: using contract-declared blendless astrometry columns "
            f"({blendless_x_col}/{blendless_y_col})."
        )
    elif col_ra_det is not None and col_dec_det is not None:
        ra_true_deg = df[col_ra_det].to_numpy(dtype=float, copy=False)
        dec_true_deg = df[col_dec_det].to_numpy(dtype=float, copy=False)
        dra_true_deg = (ra_true_deg - ra_deg + 180.0) % 360.0 - 180.0
        x_ast_true_all_arcsec = dra_true_deg * cos_dec * 3600.0
        y_ast_true_all_arcsec = (dec_true_deg - dec_deg) * 3600.0
        true_ast_source = "deterministic total-centroid RA/Dec columns"
    else:
        raise SmokeTestError(
            f"BAGLE sanity: {lc_file.name} missing usable noiseless astrometry columns "
            "(need contract-declared blendless columns, or deterministic total-centroid RA/Dec columns)."
        )
    if true_ast_source is not None:
        warnings.append(f"{lc_file.name}: noiseless astrometry source = {true_ast_source}.")

    sigma_mas = df[col_sigma].to_numpy(dtype=float, copy=False)
    x_err_mas = sigma_mas.copy()
    y_err_mas = sigma_mas.copy()
    ast_mask = (
        outside_event_window
        &
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

    t0_guess_jd = sim_zero_time + t0_event_days
    t0_guess = t0_guess_jd - 2400000.5
    tE_guess = tE_event_days
    thetaE_guess = abs(_safe_get(row, "thetaE"))
    if math.isnan(thetaE_guess) or thetaE_guess <= 0.0:
        raise SmokeTestError(
            f"BAGLE sanity: event {evt} has invalid thetaE in .out"
        )
    beta_guess = -_safe_get(row, "u0lens1") * thetaE_guess

    near_t0_idx = int(np.argmin(np.abs(t_mjd - t0_guess)))
    xS0_guess = float(x_ast_all_arcsec[near_t0_idx])
    yS0_guess = float(y_ast_all_arcsec[near_t0_idx])

    fs_guess = _safe_get(row, "Obs_0_fs")
    if not math.isfinite(fs_guess):
        fs_guess = 0.9
    fs_guess = min(max(fs_guess, 0.02), 0.999)

    far_mask = np.abs(sim_time - t0_event_days) > 3.0 * max(tE_guess, 1.0)
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

    obs_probe_params: Dict[str, float] = {
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
    obs_location_used = _select_obs_location(
        model,
        requested=obs_location,
        init_params=obs_probe_params,
        t_probe_mjd=t0_guess,
        ra_deg=ra_deg,
        dec_deg=dec_deg,
    )
    warnings.append(
        f"{lc_file.name}: BAGLE observer location set to {obs_location_used!r}."
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
        "obsLocation": obs_location_used,
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
            obs_location=obs_location_used,
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

    best_model = _build_pspl_model(
        model,
        {name: float(best[name]) for name in needed},
        ra_deg=ra_deg,
        dec_deg=dec_deg,
        obs_location=obs_location_used,
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
        deferred_failures.append(
            f"BAGLE sanity: poor joint-fit quality for {lc_file.name} "
            f"(reduced chi2={red_chi2:.3f}, threshold={fit_reduced_chi2_max:.3f}, "
            f"chi2_phot={chi2_phot:.2f}, chi2_ast={chi2_ast:.2f}, dof={dof})."
            f"{lensing_context}"
        )

    true_ast_rms_mas = float("nan")
    true_ast_sigma_equiv = float("nan")
    true_ast_epochs = 0
    if (
        true_idx_in_ast is None
        or x_ast_true_arcsec is None
        or y_ast_true_arcsec is None
        or len(true_idx_in_ast) == 0
    ):
        deferred_failures.append(
            "BAGLE sanity: noiseless astrometry was not carried through fitting; cannot perform strict model-vs-true checks."
        )
    else:
        true_ast_epochs = int(len(true_idx_in_ast))
        if len(true_idx_in_ast) < 20:
            deferred_failures.append(
                f"BAGLE sanity: only {len(true_idx_in_ast)} noiseless astrometric epochs available for strict validation."
            )
        ast_model_true = _get_model_astrometry(
            best_model,
            t_ast_true,
            lens_relative_astrometry=lens_relative_astrometry,
        )
        if len(ast_model_true) != len(t_ast_true):
            deferred_failures.append(
                f"BAGLE sanity: unexpected astrometry length for noiseless epochs: {ast_model_true.shape}"
            )
        else:
            resid_true_e_mas = (ast_model_true[:, 0] - x_ast_true_arcsec) * MAS_PER_ARCSEC
            resid_true_n_mas = (ast_model_true[:, 1] - y_ast_true_arcsec) * MAS_PER_ARCSEC
            true_ast_rms_mas = float(np.sqrt(np.mean(resid_true_e_mas**2 + resid_true_n_mas**2)))

            pair_err_mas = np.hypot(
                x_ast_err_arcsec[true_idx_in_ast] * MAS_PER_ARCSEC,
                y_ast_err_arcsec[true_idx_in_ast] * MAS_PER_ARCSEC,
            )
            finite_pair_err = pair_err_mas[np.isfinite(pair_err_mas) & (pair_err_mas > 0.0)]
            if finite_pair_err.size == 0:
                deferred_failures.append(
                    "BAGLE sanity: invalid astrometric uncertainties when evaluating model-vs-noiseless residuals."
                )
            else:
                true_ast_sigma_equiv = true_ast_rms_mas / float(np.median(finite_pair_err))
                if true_ast_rms_mas > true_ast_rms_mas_max or true_ast_sigma_equiv > true_ast_sigma_max:
                    deferred_failures.append(
                        f"BAGLE sanity: model does not track noiseless astrometry for {lc_file.name} "
                        f"(RMS_true={true_ast_rms_mas:.3f} mas, limit={true_ast_rms_mas_max:.3f} mas; "
                        f"sigma_equiv={true_ast_sigma_equiv:.2f}, limit={true_ast_sigma_max:.2f}; "
                        f"epochs={len(true_idx_in_ast)})."
                        f"{lensing_context}"
                    )

    lens_ast_rms_raw_mas: float | None = None
    lens_ast_rms_demean_mas: float | None = None
    lens_t_plot: np.ndarray | None = None
    lens_obs_x_plot_arcsec: np.ndarray | None = None
    lens_obs_y_plot_arcsec: np.ndarray | None = None
    lens_formula_x_plot_arcsec: np.ndarray | None = None
    lens_formula_y_plot_arcsec: np.ndarray | None = None
    lens_track_comparable = True
    lens_cols_raw = bagle_contract.get("lens_columns", "").strip()
    lens_frame = bagle_contract.get("lens_frame", "xy_thetaE").strip().lower()
    lens_col_x = "lens0_x"
    lens_col_y = "lens0_y"
    if lens_cols_raw:
        split_cols = [part.strip() for part in lens_cols_raw.split(",") if part.strip()]
        if len(split_cols) == 2:
            lens_col_x, lens_col_y = split_cols
        else:
            deferred_failures.append(
                f"BAGLE sanity: invalid #Astrometry_BAGLE lens_columns={lens_cols_raw!r} in {lc_file.name}; expected two comma-separated columns."
            )

    if lens_col_x in df.columns and lens_col_y in df.columns:
        if lens_frame in ("ra_dec_deg", "radec_deg", "sky_ra_dec_deg"):
            lens_ra_all = df[lens_col_x].to_numpy(dtype=float, copy=False)
            lens_dec_all = df[lens_col_y].to_numpy(dtype=float, copy=False)
            lens_dra_deg = (lens_ra_all - ra_deg + 180.0) % 360.0 - 180.0
            lens_x_all_arcsec = lens_dra_deg * cos_dec * 3600.0
            lens_y_all_arcsec = (lens_dec_all - dec_deg) * 3600.0
        elif lens_frame in ("xy_arcsec", "bagle_xy_arcsec"):
            lens_x_all_arcsec = df[lens_col_x].to_numpy(dtype=float, copy=False)
            lens_y_all_arcsec = df[lens_col_y].to_numpy(dtype=float, copy=False)
        elif lens_frame in ("xy_thetae", "bagle_xy_thetae", "thetae"):
            thetae_mas = abs(_safe_get(row, "thetaE"))
            if (not math.isfinite(thetae_mas)) or thetae_mas <= 0.0:
                deferred_failures.append(
                    f"BAGLE sanity: lens_frame={lens_frame} requested in {lc_file.name} but thetaE is invalid in .out."
                )
                lens_x_all_arcsec = np.full_like(t_ast, np.nan, dtype=float)
                lens_y_all_arcsec = np.full_like(t_ast, np.nan, dtype=float)
            else:
                conv = thetae_mas / MAS_PER_ARCSEC
                lens_x_all_arcsec = df[lens_col_x].to_numpy(dtype=float, copy=False) * conv
                lens_y_all_arcsec = df[lens_col_y].to_numpy(dtype=float, copy=False) * conv
        elif lens_frame in ("event_xy_thetae", "event_thetae", "event_frame_thetae"):
            lens_track_comparable = False
            warnings.append(
                f"{lc_file.name}: BAGLE lens_columns are event-frame thetaE diagnostics "
                f"({lens_col_x}/{lens_col_y}); skipped absolute BAGLE lens-track comparison."
            )
            lens_x_all_arcsec = np.full_like(t_ast, np.nan, dtype=float)
            lens_y_all_arcsec = np.full_like(t_ast, np.nan, dtype=float)
        else:
            deferred_failures.append(
                f"BAGLE sanity: unsupported #Astrometry_BAGLE lens_frame={lens_frame!r} in {lc_file.name}."
            )
            lens_x_all_arcsec = np.full_like(t_ast, np.nan, dtype=float)
            lens_y_all_arcsec = np.full_like(t_ast, np.nan, dtype=float)

        lens_x_sel_arcsec = lens_x_all_arcsec[ast_idx]
        lens_y_sel_arcsec = lens_y_all_arcsec[ast_idx]
        lens_finite = np.isfinite(t_ast) & np.isfinite(lens_x_sel_arcsec) & np.isfinite(lens_y_sel_arcsec)
        n_lens_finite = int(np.sum(lens_finite))
        if not lens_track_comparable:
            pass
        elif n_lens_finite == 0:
            warnings.append(
                f"{lc_file.name}: primary-lens astrometry comparison skipped because mask/quality cuts left zero usable epochs."
            )
        elif n_lens_finite < 20:
            warnings.append(
                f"{lc_file.name}: primary-lens astrometry comparison skipped because only {n_lens_finite} "
                "usable epochs remain after mask/quality cuts (<20)."
            )
        else:
            t_lens = t_ast[lens_finite]
            lens_obs_x_arcsec = lens_x_sel_arcsec[lens_finite]
            lens_obs_y_arcsec = lens_y_sel_arcsec[lens_finite]
            lens_t_plot = t_lens.copy()
            lens_obs_x_plot_arcsec = lens_obs_x_arcsec.copy()
            lens_obs_y_plot_arcsec = lens_obs_y_arcsec.copy()
            if "lens_pm_parallax_dRAcosDec_mas" in df.columns and "lens_pm_parallax_dDec_mas" in df.columns:
                lens_pm_dra_mas = df["lens_pm_parallax_dRAcosDec_mas"].to_numpy(dtype=float, copy=False)
                lens_pm_ddec_mas = df["lens_pm_parallax_dDec_mas"].to_numpy(dtype=float, copy=False)
                lens_pm_x_all_arcsec = lens_pm_dra_mas / MAS_PER_ARCSEC
                lens_pm_y_all_arcsec = lens_pm_ddec_mas / MAS_PER_ARCSEC
                lens_pm_x_sel = lens_pm_x_all_arcsec[ast_idx]
                lens_pm_y_sel = lens_pm_y_all_arcsec[ast_idx]
                lens_formula_x_plot_arcsec = lens_pm_x_sel[lens_finite]
                lens_formula_y_plot_arcsec = lens_pm_y_sel[lens_finite]
            lens_model = np.asarray(best_model.get_lens_astrometry(t_lens), dtype=float)
            if lens_model.ndim != 2 or lens_model.shape[1] < 2:
                raise SmokeTestError(
                    f"BAGLE sanity: unexpected BAGLE lens astrometry shape {lens_model.shape}."
                )
            lens_model_xy = lens_model[:, :2]
            if lens_model_xy.shape[0] != lens_obs_x_arcsec.shape[0]:
                raise SmokeTestError(
                    "BAGLE sanity: BAGLE lens astrometry length mismatch when comparing to Gulls lens columns."
                )

            lens_resid_e_mas = (lens_obs_x_arcsec - lens_model_xy[:, 0]) * MAS_PER_ARCSEC
            lens_resid_n_mas = (lens_obs_y_arcsec - lens_model_xy[:, 1]) * MAS_PER_ARCSEC
            lens_ast_rms_raw_mas = float(
                np.sqrt(np.mean(lens_resid_e_mas**2 + lens_resid_n_mas**2))
            )
            lens_resid_e_demean = lens_resid_e_mas - float(np.median(lens_resid_e_mas))
            lens_resid_n_demean = lens_resid_n_mas - float(np.median(lens_resid_n_mas))
            lens_ast_rms_demean_mas = float(
                np.sqrt(np.mean(lens_resid_e_demean**2 + lens_resid_n_demean**2))
            )
            if lens_ast_rms_demean_mas > lens_ast_rms_demean_mas_max:
                deferred_failures.append(
                    f"BAGLE sanity: primary-lens astrometry mismatch for {lc_file.name} "
                    f"(RMS_raw={lens_ast_rms_raw_mas:.3f} mas, "
                    f"RMS_after_xy_offset={lens_ast_rms_demean_mas:.3f} mas, "
                    f"limit={lens_ast_rms_demean_mas_max:.3f} mas, epochs={len(t_lens)})."
                    f"{lensing_context}"
                )
            warnings.append(
                f"{lc_file.name}: lens-track comparison vs BAGLE get_lens_astrometry "
                f"(RMS_raw={lens_ast_rms_raw_mas:.3f} mas, "
                f"RMS_after_xy_offset={lens_ast_rms_demean_mas:.3f} mas)."
            )
    else:
        available_lens_cols = [
            col for col in df.columns if col.startswith("lens") and ("_x" in col or "_y" in col)
        ]
        hint = ""
        if available_lens_cols:
            hint = (
                " Available lens-track-like columns: "
                + ", ".join(available_lens_cols[:8])
                + (" ..." if len(available_lens_cols) > 8 else "")
                + "."
            )
        if lens_frame in ("event_xy_thetae", "event_thetae", "event_frame_thetae"):
            warnings.append(
                f"{lc_file.name}: BAGLE lens-track comparison skipped because event-frame lens columns "
                f"({lens_col_x}/{lens_col_y}) were not present.{hint}"
            )
        else:
            deferred_failures.append(
                f"BAGLE sanity: {lc_file.name} missing required primary-lens astrometry columns "
                f"({lens_col_x}/{lens_col_y}; lens_frame={lens_frame}). "
                "Either publish these in the .lc output or update #Astrometry_BAGLE lens_columns/lens_frame metadata."
                + hint
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
    mu_ref = -mu_ref_raw
    warnings.append(
        "PM comparison converts .out murel_helio_* from lens-source convention to BAGLE muRel source-lens convention."
    )
    mu_amp_fit = float(np.hypot(mu_fit[0], mu_fit[1]))
    mu_amp_ref = float(np.hypot(mu_ref[0], mu_ref[1]))
    if mu_amp_ref <= 1.0e-6:
        raise SmokeTestError("BAGLE sanity: .out murel_helio vector is ~0, cannot compare direction.")
    mu_amp_frac = abs(mu_amp_fit - mu_amp_ref) / mu_amp_ref
    mu_dir_fit = math.degrees(math.atan2(mu_fit[1], mu_fit[0]))
    mu_dir_ref = math.degrees(math.atan2(mu_ref[1], mu_ref[0]))
    mu_dir_diff = abs(_angle_diff_deg(mu_dir_fit, mu_dir_ref))
    mu_mismatch_exceeds_tol = bool(
        mu_amp_frac > mu_amp_frac_tol or mu_dir_diff > mu_dir_tol_deg
    )
    if mu_mismatch_exceeds_tol:
        warnings.append(
            f"BAGLE sanity diagnostic: fitted muRel differs from .out heliocentric mu_rel "
            f"for {lc_file.name} "
            f"(fit_mu=({mu_fit[0]:.4f},{mu_fit[1]:.4f}) mas/yr, "
            f"out_mu_bagle=({mu_ref[0]:.4f},{mu_ref[1]:.4f}) mas/yr, "
            f"out_mu_raw=({mu_ref_raw[0]:.4f},{mu_ref_raw[1]:.4f}) mas/yr, "
            f"|Δamp|/amp={mu_amp_frac:.3f} (diag={mu_amp_frac_tol:.3f}), "
            f"Δdir={mu_dir_diff:.2f} deg (diag={mu_dir_tol_deg:.2f} deg))."
            f"{lensing_context}"
        )

    piE_fit = _vec_from_model_attr(best_model, "piE")
    if piE_fit is None:
        raise SmokeTestError(
            "BAGLE sanity: best-fit model has no finite piE vector; cannot compare parallax."
        )
    piE_amp_out = _safe_get(row, "piE")
    if not math.isfinite(piE_amp_out) or piE_amp_out <= 0.0:
        raise SmokeTestError(
            "BAGLE sanity: .out missing finite piE amplitude; cannot compare parallax."
        )
    mu_rel_amp_ref = math.hypot(mu_rel_e_ref, mu_rel_n_ref)
    if mu_rel_amp_ref < 1.0e-12:
        raise SmokeTestError(
            "BAGLE sanity: .out murel_helio vector is ~0; cannot derive equatorial piE direction."
        )
    piE_ref_raw = np.array([
        piE_amp_out * mu_rel_e_ref / mu_rel_amp_ref,
        piE_amp_out * mu_rel_n_ref / mu_rel_amp_ref,
    ], dtype=float)
    piE_ref = -piE_ref_raw
    warnings.append(
        "Parallax comparison derives equatorial piE from piE_amp * murel_helio direction, "
        "then converts from lens-source to BAGLE source-lens convention."
    )
    piE_amp_fit = float(np.hypot(piE_fit[0], piE_fit[1]))
    piE_amp_ref = float(np.hypot(piE_ref[0], piE_ref[1]))
    if piE_amp_ref <= 1.0e-6:
        raise SmokeTestError("BAGLE sanity: .out parallax vector is ~0, cannot compare direction.")
    piE_amp_frac = abs(piE_amp_fit - piE_amp_ref) / piE_amp_ref
    piE_dir_fit = math.degrees(math.atan2(piE_fit[1], piE_fit[0]))
    piE_dir_ref = math.degrees(math.atan2(piE_ref[1], piE_ref[0]))
    piE_dir_diff = abs(_angle_diff_deg(piE_dir_fit, piE_dir_ref))
    piE_mismatch_exceeds_tol = bool(
        piE_amp_frac > piE_amp_frac_tol or piE_dir_diff > piE_dir_tol_deg
    )
    if piE_mismatch_exceeds_tol:
        warnings.append(
            f"BAGLE sanity diagnostic: fitted piE differs from .out heliocentric piE "
            f"for {lc_file.name} "
            f"(fit_piE=({piE_fit[0]:.4f},{piE_fit[1]:.4f}), "
            f"out_piE_bagle=({piE_ref[0]:.4f},{piE_ref[1]:.4f}), "
            f"out_piE_raw=({piE_ref_raw[0]:.4f},{piE_ref_raw[1]:.4f}), "
            f"|Δamp|/amp={piE_amp_frac:.3f} (diag={piE_amp_frac_tol:.3f}), "
            f"Δdir={piE_dir_diff:.2f} deg (diag={piE_dir_tol_deg:.2f} deg))."
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
    lens_plot_candidate = fit_dir / f"event_{evt:06d}_lens_track.png"
    if lens_plot_candidate.exists():
        lens_plot_candidate.unlink()
    lens_plot_path: Path | None = None
    if (
        lens_t_plot is not None
        and lens_obs_x_plot_arcsec is not None
        and lens_obs_y_plot_arcsec is not None
    ):
        subtitle_parts: List[str] = []
        ra_evt = _safe_get(row, "ra_deg")
        dec_evt = _safe_get(row, "dec_deg")
        l_evt = _safe_get(row, "galactic_l")
        b_evt = _safe_get(row, "galactic_b")
        if math.isfinite(ra_evt) and math.isfinite(dec_evt):
            subtitle_parts.append(f"RA={ra_evt:.6f}° Dec={dec_evt:.6f}°")
        if math.isfinite(l_evt) and math.isfinite(b_evt):
            subtitle_parts.append(f"l={l_evt:.6f}° b={b_evt:.6f}°")
        mu_l = _safe_get(row, "Lens_mul")
        mu_b = _safe_get(row, "Lens_mub")
        if math.isfinite(mu_l) and math.isfinite(mu_b):
            subtitle_parts.append(f"mu_L(l,b)=({mu_l:.3f},{mu_b:.3f}) mas/yr")
        subtitle = " | ".join(subtitle_parts) if subtitle_parts else None

        lens_plot_path = lens_plot_candidate
        _plot_lens_track_diagnostics(
            lens_plot_path,
            lens_t_plot,
            lens_obs_x_plot_arcsec,
            lens_obs_y_plot_arcsec,
            lens_formula_x_plot_arcsec,
            lens_formula_y_plot_arcsec,
            best_model,
            t0_guess,
            title=(
                f"BAGLE vs Gulls primary-lens track: event {evt} "
                f"(single-lens chi2={out_chi2:.3f})"
            ),
            subtitle=subtitle,
            lens_rms_raw_mas=lens_ast_rms_raw_mas,
            lens_rms_demean_mas=lens_ast_rms_demean_mas,
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
        "obs_location": obs_location_used,
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
            "n_epochs": int(true_ast_epochs),
        },
        "proper_motion": {
            "fit": {"E": float(mu_fit[0]), "N": float(mu_fit[1]), "amp": mu_amp_fit},
            "out_helio_flipped_for_bagle_source_minus_lens_compare": {
                "E": float(mu_ref[0]),
                "N": float(mu_ref[1]),
                "amp": mu_amp_ref,
            },
            "out_helio_raw_gulls_lens_minus_source": {
                "E": float(mu_ref_raw[0]),
                "N": float(mu_ref_raw[1]),
                "amp": float(np.hypot(mu_ref_raw[0], mu_ref_raw[1])),
            },
            "amplitude_fraction_difference": mu_amp_frac,
            "direction_difference_deg": mu_dir_diff,
            "diagnostic_thresholds": {
                "amplitude_fraction": float(mu_amp_frac_tol),
                "direction_deg": float(mu_dir_tol_deg),
            },
            "exceeds_diagnostic_threshold": mu_mismatch_exceeds_tol,
        },
        "parallax": {
            "fit": {"E": float(piE_fit[0]), "N": float(piE_fit[1]), "amp": piE_amp_fit},
            "out_helio_flipped_for_bagle_source_minus_lens_compare": {
                "E": float(piE_ref[0]),
                "N": float(piE_ref[1]),
                "amp": piE_amp_ref,
            },
            "out_helio_raw_gulls_lens_minus_source": {
                "E": float(piE_ref_raw[0]),
                "N": float(piE_ref_raw[1]),
                "amp": float(np.hypot(piE_ref_raw[0], piE_ref_raw[1])),
            },
            "amplitude_fraction_difference": piE_amp_frac,
            "direction_difference_deg": piE_dir_diff,
            "diagnostic_thresholds": {
                "amplitude_fraction": float(piE_amp_frac_tol),
                "direction_deg": float(piE_dir_tol_deg),
            },
            "exceeds_diagnostic_threshold": piE_mismatch_exceeds_tol,
        },
        "best_fit_parameters": {name: float(best[name]) for name in needed if name in best},
        "best_fit_derived": {
            "tE_days": float(getattr(best_model, "tE", np.nan)),
            "thetaE_mas": float(getattr(best_model, "thetaE_amp", np.nan)),
            "piE_amp": float(getattr(best_model, "piE_amp", np.nan)),
        },
        "timing_crosscheck": {
            "timing_source_used": "SIMULATION_ZERO_TIME + Simulation_time",
            "bjd_minus_simulation_seconds": timing_bjd_diff_stats,
        },
        "lens_astrometry_comparison": {
            "rms_raw_mas": lens_ast_rms_raw_mas,
            "rms_after_xy_offset_mas": lens_ast_rms_demean_mas,
            "rms_after_xy_offset_limit_mas": float(lens_ast_rms_demean_mas_max),
        },
        "lens_plot_path": (str(lens_plot_path) if lens_plot_path is not None else None),
        "deferred_failures": deferred_failures,
        "warnings": warnings,
    }
    result_json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    if deferred_failures:
        plot_note = f" Plot: {plot_path}."
        raise SmokeTestError(
            "BAGLE sanity (deferred): diagnostics were generated before failing. "
            f"{deferred_failures[0]}{plot_note}"
        )

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
        obs_location_used=obs_location_used,
        lens_ast_rms_raw_mas=lens_ast_rms_raw_mas,
        lens_ast_rms_demean_mas=lens_ast_rms_demean_mas,
        plot_path=plot_path,
        lens_plot_path=lens_plot_path,
        result_json_path=result_json_path,
        warnings=warnings,
    )


__all__ = ["BagleJointFitSummary", "run_bagle_joint_fit_sanity"]
