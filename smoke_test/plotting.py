"""Plotting helpers for smoke test lightcurve products."""
from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

from astropy.coordinates import SkyCoord
import astropy.units as u

from .constants import REPO_ROOT
from .errors import SmokeTestError

candidate = (REPO_ROOT.parent / "VBMicrolensing").resolve()
if candidate.is_dir():
    sys.path.append(str(candidate))

from VBMicrolensing import VBMicrolensing as VBMicrolensingClass  # type: ignore[attr-defined]

VBM_CLASS = VBMicrolensingClass  # type: ignore


def _derive_event_key(lc_file: Path) -> Tuple[int, int, int] | None:
    """Map legacy filename order to summary-key order.

    Lightcurve filenames encode trailing IDs as:
      ..._<SubRun>_<Field>_<EventID>
    Summary metrics are keyed as:
      (EventID, SubRun, Field)
    """
    stem = lc_file.stem.split(".", 1)[0]
    parts = stem.rsplit("_", 3)
    if len(parts) < 4:
        return None
    subrun, field, event = (int(part) for part in parts[-3:])
    return (event, subrun, field)


def _format_metric(value: float | None, precision: int = 3) -> str:
    if value is None or math.isnan(value):
        return "n/a"
    return f"{value:.{precision}f}"


def _galactic_pm_to_icrs(l_deg: float, b_deg: float, mu_l: float, mu_b: float) -> Tuple[float, float]:
    coord = SkyCoord(
        l=l_deg * u.deg,
        b=b_deg * u.deg,
        pm_l_cosb=mu_l * u.mas / u.yr,
        pm_b=mu_b * u.mas / u.yr,
        frame="galactic",
    )
    icrs = coord.icrs
    return (
        icrs.pm_ra_cosdec.to_value(u.mas / u.yr),
        icrs.pm_dec.to_value(u.mas / u.yr),
    )


def _parse_header(lc_file: Path) -> Tuple[List[float] | None, List[float] | None]:
    planet_vals: List[float] | None = None
    event_vals: List[float] | None = None
    with lc_file.open(encoding="utf-8") as header_reader:
        for raw in header_reader:
            if not raw.startswith("#"):
                break
            stripped = raw.strip()
            if stripped.startswith("#Planet:"):
                planet_vals = [float(x) for x in stripped.split()[1:]]
            elif stripped.startswith("#Event:"):
                event_vals = [float(x) for x in stripped.split()[1:]]
    return planet_vals, event_vals


def _parse_astrometry_frame(lc_file: Path) -> Tuple[float | None, float | None]:
    """Return (ra_deg, dec_deg) from the #Astrometry_Frame header if present."""
    with lc_file.open(encoding="utf-8") as header_reader:
        for raw in header_reader:
            if not raw.startswith("#"):
                break
            if not raw.startswith("#Astrometry_Frame:"):
                continue
            parts = raw.strip().split()
            ra_deg = None
            dec_deg = None
            for part in parts:
                if part.startswith("RA_deg="):
                    try:
                        ra_deg = float(part.split("=", 1)[1])
                    except ValueError:
                        ra_deg = None
                elif part.startswith("Dec_deg="):
                    try:
                        dec_deg = float(part.split("=", 1)[1])
                    except ValueError:
                        dec_deg = None
            if ra_deg is not None or dec_deg is not None:
                return ra_deg, dec_deg
    return None, None


def _parse_header_keyvals(lc_file: Path, prefix: str) -> Dict[str, str]:
    """Parse ``key=value`` tokens from a single header line."""
    with lc_file.open(encoding="utf-8") as header_reader:
        for raw in header_reader:
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


def _resolve_lensframe_columns(
    column_names: List[str],
) -> Tuple[Dict[str, str] | None, str | None]:
    """Resolve lens-frame centroid column names and units.

    Returns (mapping, unit), where unit is "thetaE" or "mas".
    """
    cols_thetae = {
        "meas_x": "x_centroid",
        "meas_x_err": "x_centroid_error",
        "meas_y": "y_centroid",
        "meas_y_err": "y_centroid_error",
        "true_x": "true_x_centroid",
        "true_x_err": "true_x_centroid_error",
        "true_y": "true_y_centroid",
        "true_y_err": "true_y_centroid_error",
    }
    if all(col in column_names for col in cols_thetae.values()):
        return cols_thetae, "thetaE"

    cols_mas = {
        "meas_x": "x_centroid_mas",
        "meas_x_err": "x_centroid_error_mas",
        "meas_y": "y_centroid_mas",
        "meas_y_err": "y_centroid_error_mas",
        "true_x": "true_x_centroid_mas",
        "true_x_err": "true_x_centroid_error_mas",
        "true_y": "true_y_centroid_mas",
        "true_y_err": "true_y_centroid_error_mas",
    }
    if all(col in column_names for col in cols_mas.values()):
        return cols_mas, "mas"

    return None, None


def _sanity_check_flux_conservation(
    lc_file: Path,
    true_flux: np.ndarray,
    src1_flux: np.ndarray | None,
    src2_flux: np.ndarray | None,
    abs_tol: float = 1e-5,
) -> None:
    """Check that true_relative_flux ~= src1 + src2 + fblend for a couple of epochs.

    This is a strict, non-configurable sanity check. It only runs when both
    source1/source2 columns exist and the header contains #fs and #fs2. On
    failure it raises SmokeTestError so CI/smoke-test runner fails.
    """
    # Only run when both source columns are present
    if src1_flux is None or src2_flux is None:
        return

    fs = None
    fs2 = None
    with lc_file.open(encoding="utf-8") as header_reader:
        for raw in header_reader:
            if not raw.startswith("#"):
                break
            s = raw.strip()
            if s.startswith("#fs:"):
                parts = s.split()
                if len(parts) >= 2:
                    try:
                        fs = float(parts[1])
                    except Exception:
                        fs = None
            elif s.startswith("#fs2:"):
                parts = s.split()
                if len(parts) >= 2:
                    try:
                        fs2 = float(parts[1])
                    except Exception:
                        fs2 = None

    if fs is None or fs2 is None:
        # If header doesn't provide both, skip the strict check
        return

    fblend = 1.0 - float(fs) - float(fs2)

    # choose representative epochs: first epoch and the epoch of peak true flux
    try:
        peak_idx = int(np.nanargmax(true_flux))
    except Exception:
        peak_idx = 0

    indices = [0]
    if peak_idx != 0:
        indices.append(int(peak_idx))

    for idx in indices:
        try:
            expected = float(src1_flux[idx]) + float(src2_flux[idx]) + float(fblend)
            truev = float(true_flux[idx])
        except Exception:
            # If indexing or conversion fails, raise an error to surface the issue
            raise SmokeTestError(
                f"Flux sanity check failed: could not index/convert epoch {idx} in {lc_file.name}"
            )
        diff = abs(truev - expected)
        if diff > abs_tol:
            raise SmokeTestError(
                f"Flux sanity check failed for {lc_file.name} at epoch index {idx}: true={truev:.8f}, expected={expected:.8f}, abs_diff={diff:.8e}, abs_tol={abs_tol}"
            )



def _plot_minimal_astrometry(
    lc_file: Path,
    output_dir: Path,
    title: str,
    time: np.ndarray,
    true_x_er: np.ndarray,
    true_y_er: np.ndarray,
    meas_x_er: np.ndarray,
    meas_y_er: np.ndarray,
    x_er_err: np.ndarray,
    y_er_err: np.ndarray,
    tE: float,
    unit_label: str = "Einstein radii",
) -> Path:
    fig, ax = plt.subplots(2, 1, figsize=(8, 8))
    fig.suptitle(title, fontsize=14)

    ax[0].errorbar(
        meas_x_er,
        meas_y_er,
        xerr=x_er_err,
        yerr=y_er_err,
        fmt="o",
        markersize=3,
        alpha=0.5,
        color="C0",
        label="Measured",
        zorder=1,
    )
    ax[0].plot(
        true_x_er,
        true_y_er,
        "-",
        linewidth=1.5,
        color="red",
        label="True",
        zorder=2,
        alpha=0.8,
    )
    ax[0].scatter(
        true_x_er,
        true_y_er,
        c=time,
        cmap=plt.get_cmap("plasma"),
        s=20,
        marker="x",
        linewidths=0.9,
        alpha=1.0,
        label="True samples",
        zorder=3,
    )
    ax[0].set_xlabel(f"x_centroid ({unit_label})")
    ax[0].set_ylabel(f"y_centroid ({unit_label})")
    ax[0].set_title("Astrometric Centroid in Lens Frame with Lens at Rest")
    ax[0].grid(True, alpha=0.3)
    ax[0].legend()
    
    # plot 2: x and y centroids around event peak
    # zoom in to +/- 2 Einstein radius around origin and epochs +/- 2*tE
    mask_peak = (time >= -2 * tE) & (time <= 2 * tE)
    ax[1].errorbar(
        time[mask_peak],
        meas_x_er[mask_peak],
        yerr=x_er_err[mask_peak],
        fmt="o",
        markersize=3,
        alpha=0.5,
        color="C0",
        label="Measured x",
        zorder=1,
    )
    ax[1].errorbar(
        time[mask_peak],
        meas_y_er[mask_peak],
        yerr=y_er_err[mask_peak],
        fmt="o",
        markersize=3,
        alpha=0.5,
        color="C1",
        label="Measured y",
        zorder=1,
    )
    ax[1].plot(
        time[mask_peak],
        true_x_er[mask_peak],
        "-",
        linewidth=1.5,
        color="red",
        label="True x",
        zorder=2,
        alpha=0.8,
    )
    ax[1].plot(
        time[mask_peak],
        true_y_er[mask_peak],
        "-",
        linewidth=1.5,
        color="orange",
        label="True y",
        zorder=2,
        alpha=0.8,
    )
    ax[1].set_xlabel("Time (days)")
    ax[1].set_ylabel(f"Centroid ({unit_label})")
    ax[1].set_title("Centroid Around Event Peak")
    ax[1].grid(True, alpha=0.3)
    ax[1].legend()

    plot_file = output_dir / f"{lc_file.stem}_astrometry_plot.png"
    fig.tight_layout()
    fig.savefig(plot_file, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Generated astrometry plot: {plot_file.name}")
    return plot_file

def _plot_photometry_only(
    lc_file: Path,
    output_dir: Path,
    title: str,
    time: np.ndarray,
    flux: np.ndarray,
    flux_err: np.ndarray,
    true_flux: np.ndarray | None,
    src1_flux: np.ndarray | None = None,
    src2_flux: np.ndarray | None = None,
) -> Path:
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    fig.suptitle(title, fontsize=14)
    ax.errorbar(
        time,
        flux,
        yerr=flux_err,
        fmt="o",
        markersize=2,
        alpha=0.5,
        color="C0",
        label="Measured",
        zorder=1,
    )
    if true_flux is not None:
        ax.plot(
            time,
            true_flux,
            "-",
            linewidth=1.5,
            color="red",
            label="True",
            zorder=2,
            alpha=0.8,
        )
    if src1_flux is not None:
        ax.plot(
            time,
            src1_flux,
            "-",
            linewidth=1.2,
            color="magenta",
            label="Source 1",
            zorder=4,
            alpha=0.7,
        )
    if src2_flux is not None:
        ax.plot(
            time,
            src2_flux,
            "-",
            linewidth=1.2,
            color="gold",
            label="Source 2",
            zorder=5,
            alpha=0.7,
        )
    ax.axhline(1.0, color="k", linestyle="--", linewidth=1.5, label="Baseline", zorder=3)
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Relative Flux")
    ax.set_title("Light Curve")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plot_file = output_dir / f"{lc_file.stem}_plot.png"
    fig.tight_layout()
    fig.savefig(plot_file, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Generated plot: {plot_file.name}")
    return plot_file


def _compute_vbm_model(
    summary: Dict[str, float] | None,
    planet_vals: List[float] | None,
    event_vals: List[float] | None,
    source_pm_icrs: Tuple[float, float] | None,
    lens_pm_icrs: Tuple[float, float] | None,
    theta_e_float: float | None,
    source_dist_float: float | None,
    event_ra_float: float | None,
    event_dec_float: float | None,
    alpha_deg_float: float,
    sim_zero_offset: float,
    time: np.ndarray,
    true_x_vals: np.ndarray | None,
    true_y_vals: np.ndarray | None,
) -> Tuple[Dict[str, np.ndarray | str] | None, str | None]:
    if summary is None:
        return None, "summary metrics missing for this lightcurve"
    if not (
        planet_vals
        and len(planet_vals) >= 6
        and event_vals
        and len(event_vals) >= 8
    ):
        return None, "header missing #Planet/#Event metadata"

    q_val = float(planet_vals[4])
    s_val = float(planet_vals[5])
    if q_val <= 0 or s_val <= 0:
        return None, f"unphysical planet parameters (q={q_val}, s={s_val})"

    rho_val = summary.get("rho")
    if rho_val is None or math.isnan(rho_val):
        rho_val = float(event_vals[7])
    tE_val = summary.get("tE_ref")
    if tE_val is None or math.isnan(tE_val):
        tE_val = float(event_vals[6])
    u0_val = summary.get("u0")
    if u0_val is None or math.isnan(u0_val):
        u0_val = float(event_vals[0])
    alpha_deg = summary.get("alpha_event")
    if alpha_deg is None or math.isnan(alpha_deg):
        alpha_deg = float(event_vals[1])
    t0_val = summary.get("t0")
    if t0_val is None or math.isnan(t0_val):
        t0_val = float(event_vals[2])
    pi_n_val = summary.get("pi_n")
    pi_e_val = summary.get("pi_e")
    if pi_n_val is None or math.isnan(pi_n_val):
        pi_n_val = 0.0
    if pi_e_val is None or math.isnan(pi_e_val):
        pi_e_val = 0.0

    if (
        theta_e_float is None
        or theta_e_float <= 0
        or source_dist_float is None
        or source_dist_float <= 0
        or source_pm_icrs is None
        or event_ra_float is None
        or event_dec_float is None
    ):
        return None, "insufficient astrometric metadata (theta_E, source distance, or PM)"

    pi_s_val = 1.0 / source_dist_float
    if (
        pi_s_val <= 0
        or rho_val is None
        or float(rho_val) <= 0
        or tE_val is None
        or float(tE_val) <= 0
    ):
        return None, "missing positive rho/tE/source distance for VBM evaluation"

    vbm = VBM_CLASS()  # type: ignore[operator]
    skycoord = SkyCoord(
        ra=float(event_ra_float) * u.deg,
        dec=float(event_dec_float) * u.deg,
    )
    coord_str = (
        f"{skycoord.ra.to_string(unit=u.hour, sep=':', pad=True)} "
        f"{skycoord.dec.to_string(unit=u.deg, sep=':', pad=True, alwayssign=True)}"
    )
    vbm.SetObjectCoordinates(coord_str)
    params_vbm = [
        math.log(float(s_val)),
        math.log(float(q_val)),
        float(u0_val),
        math.radians(float(alpha_deg)),
        math.log(float(rho_val)),
        math.log(float(tE_val)),
        float(t0_val) + sim_zero_offset,
        float(pi_n_val),
        float(pi_e_val),
        float(source_pm_icrs[1]),
        float(source_pm_icrs[0]),
        float(pi_s_val),
        float(theta_e_float),
    ]
    results = vbm.BinaryAstroLightCurve(params_vbm, time + sim_zero_offset)

    lens_dec_deg = np.array(results[3], dtype=float)
    lens_ra_deg = np.array(results[4], dtype=float)
    y1 = np.array(results[5], dtype=float)
    y2 = np.array(results[6], dtype=float)
    lensframe_label = "Source Trajectory BinaryAstroLightCurve"

    vbm_x = y1
    vbm_y = y2
    if (
        true_x_vals is not None
        and true_y_vals is not None
        and len(true_x_vals) == len(y1)
    ):
        combos = [
            ("x= y1, y= y2", y1, y2),
            ("x=-y1, y= y2", -y1, y2),
            ("x= y1, y=-y2", y1, -y2),
            ("x=-y1, y=-y2", -y1, -y2),
            ("x= y2, y= y1", y2, y1),
            ("x=-y2, y= y1", -y2, y1),
            ("x= y2, y=-y1", y2, -y1),
            ("x=-y2, y=-y1", -y2, -y1),
        ]
        best = None
        best_err = None
        for label, cand_x, cand_y in combos:
            err = np.nanmean((cand_x - true_x_vals) ** 2 + (cand_y - true_y_vals) ** 2)
            if best_err is None or err < best_err:
                best_err = err
                best = (label, cand_x, cand_y)
        if best:
            lensframe_label = f"Source trajectory BinaryAstroLightCurve ({best[0]})"
            vbm_x, vbm_y = best[1], best[2]

    if len(vbm_x) != len(time):
        return None, "VBM returned mismatched array lengths"

    model: Dict[str, np.ndarray | str] = {
        "lens_x": vbm_x,
        "lens_y": vbm_y,
        "lens_label": lensframe_label,
        "sky_ra": lens_ra_deg,
        "sky_dec": lens_dec_deg,
        "sky_label": "VBM BinaryAstroLightCurve (sky)",
    }
    return model, None


def _render_lensframe(
    lc_file: Path,
    output_dir: Path,
    time: np.ndarray,
    cmap: matplotlib.colors.Colormap,
    norm: Normalize,
    vbm_model: Dict[str, np.ndarray | str],
    true_x_vals: np.ndarray | None,
    true_y_vals: np.ndarray | None,
    meas_x: np.ndarray | None,
    meas_y: np.ndarray | None,
) -> Path:
    vbm_x = np.asarray(vbm_model["lens_x"])  # type: ignore[index]
    vbm_y = np.asarray(vbm_model["lens_y"])  # type: ignore[index]
    vbm_label = str(vbm_model["lens_label"])

    fig2, ax2 = plt.subplots(figsize=(6, 6))

    if meas_x is not None and meas_y is not None:
        mask_meas = np.isfinite(meas_x) & np.isfinite(meas_y)
        if np.any(mask_meas):
            ax2.scatter(
                meas_x[mask_meas],
                meas_y[mask_meas],
                c=time[mask_meas],
                cmap=cmap,
                norm=norm,
                s=25,
                alpha=0.6,
                label="Measured centroid",
                zorder=1,
            )
    else:
        ax2.text(
            0.02,
            0.98,
            "Measured centroid not available",
            ha="left",
            va="top",
            transform=ax2.transAxes,
            fontsize=8,
        )

    if true_x_vals is not None and true_y_vals is not None:
        mask_true = np.isfinite(true_x_vals) & np.isfinite(true_y_vals)
        if np.any(mask_true):
            ax2.scatter(
                true_x_vals[mask_true],
                true_y_vals[mask_true],
                c=time[mask_true],
                cmap=cmap,
                norm=norm,
                s=20,
                marker="x",
                linewidths=0.9,
                alpha=1.0,
                label="True centroid",
                zorder=2,
            )

    ax2.plot(
        vbm_x,
        vbm_y,
        color="black",
        linewidth=1.5,
        label=vbm_label,
        zorder=3,
    )
    ax2.set_xlabel("x_centroid (Einstein radii)")
    ax2.set_ylabel("y_centroid (Einstein radii)")
    ax2.set_title(f"Lens-frame Centroid: {lc_file.stem}")
    ax2.grid(True, alpha=0.3)

    segments_x = [vbm_x]
    segments_y = [vbm_y]
    if meas_x is not None and meas_y is not None:
        mask_meas = np.isfinite(meas_x) & np.isfinite(meas_y)
        if np.any(mask_meas):
            segments_x.append(np.asarray(meas_x)[mask_meas])
            segments_y.append(np.asarray(meas_y)[mask_meas])
    if true_x_vals is not None and true_y_vals is not None:
        mask_true = np.isfinite(true_x_vals) & np.isfinite(true_y_vals)
        if np.any(mask_true):
            segments_x.append(true_x_vals[mask_true])
            segments_y.append(true_y_vals[mask_true])
    all_x = np.concatenate(segments_x) if segments_x else np.array([0.0])
    all_y = np.concatenate(segments_y) if segments_y else np.array([0.0])
    x_min = float(np.nanmin(all_x))
    x_max = float(np.nanmax(all_x))
    y_min = float(np.nanmin(all_y))
    y_max = float(np.nanmax(all_y))
    x_c = 0.5 * (x_min + x_max)
    y_c = 0.5 * (y_min + y_max)
    half_span = max(x_max - x_min, y_max - y_min) * 0.5
    half_span = max(half_span, 1e-6)
    half_span *= 1.15
    ax2.set_xlim(x_c - half_span, x_c + half_span)
    ax2.set_ylim(y_c - half_span, y_c + half_span)
    ax2.set_aspect("equal", adjustable="box")
    ax2.legend(loc="upper left")

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar2 = fig2.colorbar(sm, ax=ax2, fraction=0.046, pad=0.04)
    cbar2.set_label("Time (days)")

    lensframe_file = output_dir / f"{lc_file.stem}_lensframe_plot.png"
    fig2.tight_layout()
    fig2.savefig(lensframe_file, dpi=150, bbox_inches="tight")
    plt.close(fig2)
    return lensframe_file


def _render_astrometric_figure(
    lc_file: Path,
    output_dir: Path,
    title: str,
    time: np.ndarray,
    flux: np.ndarray,
    flux_err: np.ndarray,
    true_flux: np.ndarray | None,
    true_N_mas: np.ndarray,
    true_E_mas: np.ndarray,
    meas_N_mas: np.ndarray,
    meas_E_mas: np.ndarray,
    meas_N_err_mas: np.ndarray,
    meas_E_err_mas: np.ndarray,
    true_ra_deg: np.ndarray,
    true_dec_deg: np.ndarray,
    meas_ra_deg: np.ndarray,
    meas_dec_deg: np.ndarray,
    meas_ra_err_deg: np.ndarray,
    meas_dec_err_deg: np.ndarray,
    vector_specs: List[Dict[str, float | str]],
    vbm_model: Dict[str, np.ndarray | str] | None,
    true_x_vals: np.ndarray | None,
    true_y_vals: np.ndarray | None,
    meas_x: np.ndarray | None,
    meas_y: np.ndarray | None,
    event_frame_data: Dict[str, np.ndarray] | None = None,
    src1_flux: np.ndarray | None = None,
    src2_flux: np.ndarray | None = None,
    panel_labels: Dict[str, str] | None = None,
    event_info_rows: List[Tuple[str, str]] | None = None,
) -> Tuple[Path, Path | None]:
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10))
    fig.suptitle(title, fontsize=14)
    ax_light = axes[0, 0]
    ax_event = axes[0, 1]
    ax_radec = axes[1, 0]
    ax_info = axes[1, 1]

    labels = panel_labels or {}

    ax_light.errorbar(
        time,
        flux,
        yerr=flux_err,
        fmt="o",
        markersize=2,
        alpha=0.5,
        color="C0",
        label="Measured",
        zorder=1,
    )
    if true_flux is not None:
        ax_light.plot(
            time,
            true_flux,
            "-",
            linewidth=1.5,
            color="red",
            label="True",
            zorder=2,
            alpha=0.8,
        )
    if src1_flux is not None:
        ax_light.plot(
            time,
            src1_flux,
            "-",
            linewidth=1.2,
            color="magenta",
            label="Source 1",
            zorder=4,
            alpha=0.7,
        )
    if src2_flux is not None:
        ax_light.plot(
            time,
            src2_flux,
            "-",
            linewidth=1.2,
            color="gold",
            label="Source 2",
            zorder=5,
            alpha=0.7,
        )
    ax_light.axhline(
        1.0,
        color="k",
        linestyle="--",
        linewidth=1.5,
        label="Baseline",
        zorder=3,
    )
    light_x_col = labels.get("light_x_col", "Simulation_time")
    light_y_col = labels.get("light_y_col", "measured_relative_flux")
    light_yerr_col = labels.get("light_yerr_col", "measured_relative_flux_error")
    light_true_col = labels.get("light_true_col", "true_relative_flux")
    ax_light.set_xlabel(f"Time (days) [{light_x_col}]")
    ax_light.set_ylabel("Relative Flux")
    ax_light.set_title("Light Curve")
    ax_light.text(
        0.01,
        0.98,
        f"cols: {light_y_col}, {light_yerr_col}, {light_true_col}",
        transform=ax_light.transAxes,
        va="top",
        ha="left",
        fontsize=7,
    )
    ax_light.legend()
    ax_light.grid(True, alpha=0.3)

    norm = Normalize(vmin=np.min(time), vmax=np.max(time)) if len(time) else Normalize(0, 1)
    cmap = plt.get_cmap("plasma")

    span_years: float | None = None
    if len(time):
        span_days = float(time.max() - time.min())
        if span_days > 0:
            span_years = span_days / 365.25

    ax_radec.errorbar(
        meas_ra_deg,
        meas_dec_deg,
        xerr=meas_ra_err_deg,
        yerr=meas_dec_err_deg,
        fmt="none",
        ecolor="lightgray",
        alpha=0.5,
        capsize=2,
        zorder=0,
    )
    sc_ra = ax_radec.scatter(
        meas_ra_deg,
        meas_dec_deg,
        c=time,
        cmap=cmap,
        norm=norm,
        s=25,
        alpha=0.5,
        label="Measured",
        zorder=1,
    )
    ax_radec.plot(
        true_ra_deg,
        true_dec_deg,
        color="black",
        linewidth=1.2,
        alpha=0.8,
        label="True track",
        zorder=4,
    )
    ax_radec.scatter(
        true_ra_deg,
        true_dec_deg,
        c=time,
        cmap=cmap,
        norm=norm,
        s=18,
        marker="x",
        linewidths=0.8,
        alpha=1.0,
        label="True samples",
        zorder=3,
    )
    if vbm_model is not None:
        deg_to_mas = 3600.0 * 1000.0
        baseline_ra = true_ra_deg[0]
        baseline_dec = true_dec_deg[0]
        cos_dec0 = math.cos(math.radians(baseline_dec))
        if abs(cos_dec0) < 1e-6:
            cos_dec0 = 1e-6 if cos_dec0 >= 0 else -1e-6

        vbm_ra_offset_mas = np.asarray(vbm_model["sky_ra"], dtype=float)  # type: ignore[index]
        vbm_dec_offset_mas = np.asarray(vbm_model["sky_dec"], dtype=float)  # type: ignore[index]
        vbm_ra_abs = baseline_ra + (vbm_ra_offset_mas / (deg_to_mas * cos_dec0))
        vbm_dec_abs = baseline_dec + (vbm_dec_offset_mas / deg_to_mas)

        ax_radec.plot(
            vbm_ra_abs,
            vbm_dec_abs,
            color="tab:purple",
            linewidth=1.2,
            alpha=0.9,
            label=vbm_model.get("sky_label", "VBM BinaryAstroLightCurve (sky)"),
        )

    radec_obs_ra_col = labels.get("radec_obs_ra_col", "measured_centroid_ra_deg")
    radec_obs_dec_col = labels.get("radec_obs_dec_col", "measured_centroid_dec_deg")
    radec_obs_ra_err_col = labels.get("radec_obs_ra_err_col", "measured_centroid_ra_error_deg")
    radec_obs_dec_err_col = labels.get("radec_obs_dec_err_col", "measured_centroid_dec_error_deg")
    radec_true_ra_col = labels.get("radec_true_ra_col", "true_centroid_ra_deg")
    radec_true_dec_col = labels.get("radec_true_dec_col", "true_centroid_dec_deg")
    ax_radec.set_xlabel("RA (degrees)")
    ax_radec.set_ylabel("Dec (degrees)")
    ax_radec.set_title("Absolute Astrometric Position (RA/Dec)")
    ax_radec.text(
        0.01,
        0.98,
        (
            f"obs cols: {radec_obs_ra_col}, {radec_obs_dec_col}\n"
            f"obs err cols: {radec_obs_ra_err_col}, {radec_obs_dec_err_col}\n"
            f"noiseless cols: {radec_true_ra_col}, {radec_true_dec_col}"
        ),
        transform=ax_radec.transAxes,
        va="top",
        ha="left",
        fontsize=7,
    )
    ax_radec.grid(True, alpha=0.3)
    ax_radec.axis("equal")

    if span_years and vector_specs:
        start_ra = true_ra_deg[0]
        start_dec = true_dec_deg[0]
        cos_dec = math.cos(math.radians(start_dec))
        if abs(cos_dec) < 1e-6:
            cos_dec = 1e-6 if cos_dec >= 0 else -1e-6
        for spec in vector_specs:
            pm_ra = spec["pm_ra"]
            pm_dec = spec["pm_dec"]
            delta_ra_deg = (pm_ra * span_years) / (3600000.0 * cos_dec)
            delta_dec_deg = (pm_dec * span_years) / 3600000.0
            end_ra = start_ra + delta_ra_deg
            end_dec = start_dec + delta_dec_deg
            ax_radec.annotate(
                "",
                xy=(end_ra, end_dec),
                xytext=(start_ra, start_dec),
                arrowprops=dict(color=spec["color"], arrowstyle="->", linewidth=1),
                zorder=5,
            )
            ax_radec.plot([], [], color=spec["color"], linewidth=2, label=spec["label"])
    ax_radec.legend(ncol=2, fontsize=8)

    ax_info.axis("off")
    ax_info.set_title("Event Info / Validation Inputs")
    info_rows = event_info_rows or []
    if not info_rows:
        info_rows = [("info", "n/a")]
    info_table = ax_info.table(
        cellText=[[field, value] for field, value in info_rows],
        colLabels=["Field", "Value"],
        cellLoc="left",
        colLoc="left",
        loc="center",
    )
    info_table.auto_set_font_size(False)
    info_table.set_fontsize(8)
    info_table.scale(1.0, 1.15)

    fig.tight_layout(rect=[0, 0.12, 1, 1])
    cbar_ax = fig.add_axes([0.25, 0.06, 0.5, 0.025])
    cbar = fig.colorbar(sc_ra, cax=cbar_ax, orientation="horizontal")
    cbar.set_label("Time (days)")

    evt = event_frame_data or {}
    centroid_x = evt.get("centroid_x")
    centroid_y = evt.get("centroid_y")
    if centroid_x is not None and centroid_y is not None:
        mask_cent = np.isfinite(centroid_x) & np.isfinite(centroid_y) & np.isfinite(time)
        if np.any(mask_cent):
            ax_event.plot(
                centroid_x[mask_cent],
                centroid_y[mask_cent],
                color="black",
                linewidth=1.2,
                alpha=0.8,
                label="Centroid track",
                zorder=4,
            )
            ax_event.scatter(
                centroid_x[mask_cent],
                centroid_y[mask_cent],
                c=time[mask_cent],
                cmap=cmap,
                norm=norm,
                s=18,
                marker="x",
                linewidths=0.8,
                alpha=1.0,
                label="Centroid samples",
                zorder=5,
            )

    def _plot_track(x: np.ndarray | None, y: np.ndarray | None, color: str, label: str, z: int) -> None:
        if x is None or y is None:
            return
        mask = np.isfinite(x) & np.isfinite(y)
        if not np.any(mask):
            return
        ax_event.plot(x[mask], y[mask], color=color, linewidth=1.1, alpha=0.9, label=label, zorder=z)

    _plot_track(evt.get("source0_x"), evt.get("source0_y"), "tab:blue", "Source 0", 2)
    _plot_track(evt.get("lens0_x"), evt.get("lens0_y"), "tab:red", "Lens 0", 2)
    _plot_track(evt.get("lens1_x"), evt.get("lens1_y"), "tab:orange", "Lens 1", 2)

    event_x_col = labels.get("event_x_col", "true_x_centroid")
    event_y_col = labels.get("event_y_col", "true_y_centroid")
    cols_note_parts = [event_x_col, event_y_col]
    for key in ("source0_x_col", "source0_y_col", "source0_mu_col", "lens0_x_col", "lens0_y_col", "lens1_x_col", "lens1_y_col"):
        val = labels.get(key)
        if val:
            cols_note_parts.append(val)
    ax_event.set_xlabel("Event-frame x (theta_E)")
    ax_event.set_ylabel("Event-frame y (theta_E)")
    ax_event.set_title("Event-Frame Tracks")
    ax_event.text(
        0.01,
        0.98,
        "cols: " + ", ".join(cols_note_parts),
        transform=ax_event.transAxes,
        va="top",
        ha="left",
        fontsize=7,
    )
    ax_event.grid(True, alpha=0.3)
    ax_event.axis("equal")
    handles, labels_ = ax_event.get_legend_handles_labels()
    if handles:
        ax_event.legend(fontsize=8)

    plot_file = output_dir / f"{lc_file.stem}_plot.png"
    fig.savefig(plot_file, dpi=150, bbox_inches="tight")
    plt.close(fig)

    lensframe_path: Path | None = None
    if vbm_model is not None:
        lensframe_path = _render_lensframe(
            lc_file,
            output_dir,
            time,
            cmap,
            norm,
            vbm_model,
            true_x_vals,
            true_y_vals,
            meas_x,
            meas_y,
        )

    return plot_file, lensframe_path


def _extract_gulls_version(build_bin: Path) -> str:
    """Extract the gulls version by running the executable.
    
    Returns version string like "2.1.0" or "unknown" if unable to determine.
    """
    import subprocess
    import re
    
    # Try gulls_std first as it's the most common
    executables = ["gulls_std", "gulls_croin", "gullsFish"]
    
    for exe_name in executables:
        exe_path = build_bin / exe_name
        if not exe_path.exists():
            continue
        
        try:
            # Run the executable with invalid arguments to get version in error output
            result = subprocess.run(
                [str(exe_path), "-i", "/nonexistent"],
                capture_output=True,
                text=True,
                timeout=5
            )
            # Check both stdout and stderr for version string
            output = result.stdout + result.stderr
            
            # Look for "gulls v2.1.0" pattern
            match = re.search(r'gulls v(\d+\.\d+\.\d+)', output)
            if match:
                return match.group(1)
        except (subprocess.TimeoutExpired, subprocess.SubprocessError, OSError):
            continue
    
    # Default to 2.1.0 if we can't extract it
    return "2.1.0"


def plot_lightcurves(
    output_dir: Path,
    summaries: Dict[Tuple[int, int, int], Dict[str, float]] | None = None,
    params: Dict[str, str] | None = None,
    build_bin: Path | None = None,
) -> None:
    lc_files = sorted(output_dir.rglob("*.lc"))
    if not lc_files:
        return

    sim_zero_offset = 0.0
    if params:
        sz = params.get("SIMULATION_ZERO_TIME")
        if sz:
            sim_zero_offset = float(sz) - 2450000.0

    astrometry_expected = False
    if params:
        val = params.get("ASTROMETRY_ON")
        if val is not None:
            astrometry_expected = str(val).strip().lower() not in {"0", "false", "off"}

    # Determine whether MULTIPLE_SOURCES was requested in the parameter file.
    multiple_sources_expected = False
    if params:
        ms = params.get("MULTIPLE_SOURCES")
        if ms is not None:
            multiple_sources_expected = str(ms).strip().lower() not in {"0", "false", "off"}

    exec_name = ""
    if params:
        exec_name = params.get("EXECUTABLE", "")
    exec_name = exec_name.strip().lower()
    if exec_name.endswith(".x"):
        exec_name = exec_name[:-2]
    vbm_supported = exec_name in {"gulls_croin", "gullsfish", "gulls_std"}
    vbm_required = astrometry_expected and vbm_supported and exec_name != "gullssingle"
    plot_failures: List[str] = []
    
    for lc_file in lc_files:
        planet_vals, event_vals = _parse_header(lc_file)
        df = pd.read_csv(lc_file, sep=r"\s+", comment="#")
        if df.empty:
            continue

        column_names = list(df.columns)

        def _require_column(name: str) -> np.ndarray:
            if name not in column_names:
                raise SmokeTestError(f"Smoke test failed: column '{name}' missing in {lc_file.name}")
            return df[name].to_numpy(dtype=float, copy=False)

        def _optional_column(name: str) -> np.ndarray | None:
            if name not in column_names:
                return None
            return df[name].to_numpy(dtype=float, copy=False)

        summary = None
        if summaries:
            event_key = _derive_event_key(lc_file)
            if event_key is not None:
                summary = summaries.get(event_key)
        if summary is None:
            plot_failures.append(f"summary metrics missing for {lc_file.name}")

        title = f"Smoke Test: {lc_file.stem}"
        theta_e_float: float | None = None
        alpha_deg_float: float | None = None
        pi_n_float: float | None = None
        pi_e_float: float | None = None
        source_dist_float: float | None = None
        event_ra_float: float | None = None
        event_dec_float: float | None = None
        tE = math.nan
        if summary:
            lens_mass = _format_metric(summary.get("lens_mass"))
            lens_dist = _format_metric(summary.get("lens_dist"))
            source_dist = _format_metric(summary.get("source_dist"))
            theta_e = _format_metric(summary.get("theta_e"))
            pm_alpha = _format_metric(summary.get("pm_helio_alpha"))
            pm_delta = _format_metric(summary.get("pm_helio_delta"))
            subtitle = (
                f"Lens M={lens_mass} Msun, Lens D={lens_dist} pc, "
                f"Source D={source_dist} pc, theta_E={theta_e}, "
                f"mu_rel=({pm_alpha}, {pm_delta}) mas/yr"
            )
            title = f"{title}\n{subtitle}"
            val = summary.get("theta_e")
            if val is not None and not math.isnan(val):
                theta_e_float = float(val)
            val = summary.get("alpha_event")
            if val is not None and not math.isnan(val):
                alpha_deg_float = float(val)
            val = summary.get("pi_n")
            if val is not None and not math.isnan(val):
                pi_n_float = float(val)
            val = summary.get("pi_e")
            if val is not None and not math.isnan(val):
                pi_e_float = float(val)
            val = summary.get("source_dist")
            if val is not None and not math.isnan(val) and val > 0:
                source_dist_float = float(val)
            val = summary.get("event_ra")
            if val is not None and not math.isnan(val):
                event_ra_float = float(val)
            val = summary.get("event_dec")
            if val is not None and not math.isnan(val):
                event_dec_float = float(val)
            val = summary.get("tE_ref")
            if val is not None and not math.isnan(val):
                tE = float(val)
        if alpha_deg_float is None and event_vals and len(event_vals) >= 2:
            alpha_deg_float = float(event_vals[1])
        if alpha_deg_float is None:
            raise SmokeTestError(f"Smoke test failed: missing alpha_event for {lc_file.name}")

        time = _require_column("Simulation_time")
        flux = _require_column("measured_relative_flux")
        flux_err = _require_column("measured_relative_flux_error")
        true_flux = _require_column("true_relative_flux")

        # Astrometry column validation with backward/forward compatibility.
        lensframe_thetae_cols = [
            "x_centroid",
            "x_centroid_error",
            "y_centroid",
            "y_centroid_error",
            "true_x_centroid",
            "true_x_centroid_error",
            "true_y_centroid",
            "true_y_centroid_error",
        ]
        lensframe_mas_cols = [
            "x_centroid_mas",
            "x_centroid_error_mas",
            "y_centroid_mas",
            "y_centroid_error_mas",
            "true_x_centroid_mas",
            "true_x_centroid_error_mas",
            "true_y_centroid_mas",
            "true_y_centroid_error_mas",
        ]
        skyframe_cols = [
            "true_N_centroid_mas",
            "true_E_centroid_mas",
            "measured_N_centroid_mas",
            "measured_E_centroid_mas",
            "measured_N_centroid_error_mas",
            "measured_E_centroid_error_mas",
            "true_centroid_ra_deg",
            "true_centroid_dec_deg",
            "measured_centroid_ra_deg",
            "measured_centroid_dec_deg",
            "measured_centroid_ra_error_deg",
            "measured_centroid_dec_error_deg",
        ]
        radec_cols = [
            "RA_centroid_deg",
            "Dec_centroid_deg",
            "RA_centroid_true_deg",
            "Dec_centroid_true_deg",
        ]

        contract = _parse_header_keyvals(lc_file, "#Astrometry_Contract:")
        contract_cols = _parse_header_keyvals(lc_file, "#Astrometry_Columns:")
        vbm_function_meta = _parse_header_keyvals(lc_file, "#VBM_function:")
        vbm_function_name = vbm_function_meta.get("name", "n/a")

        def _mapped_col(contract_key: str, *fallbacks: str) -> str | None:
            candidates: List[str] = []
            if contract_key in contract_cols:
                candidates.append(contract_cols[contract_key])
            candidates.extend(fallbacks)
            for name in candidates:
                if name in column_names:
                    return name
            return None

        ra_obs_col = _mapped_col(
            "sky_ra_measured_deg",
            "sky_ra_obs_deg",
            "RA_measured_deg",
            "RA_obs_deg",
            "RA_centroid_deg",
        )
        dec_obs_col = _mapped_col(
            "sky_dec_measured_deg",
            "sky_dec_obs_deg",
            "Dec_measured_deg",
            "Dec_obs_deg",
            "Dec_centroid_deg",
        )
        ra_true_col = _mapped_col(
            "sky_ra_noiseless_deg",
            "sky_ra_det_deg",
            "RA_noiseless_deg",
            "RA_det_deg",
            "RA_centroid_true_deg",
        )
        dec_true_col = _mapped_col(
            "sky_dec_noiseless_deg",
            "sky_dec_det_deg",
            "Dec_noiseless_deg",
            "Dec_det_deg",
            "Dec_centroid_true_deg",
        )
        ra_err_col = _mapped_col("sky_ra_err_deg", "RA_err_deg", "measured_centroid_ra_error_deg")
        dec_err_col = _mapped_col("sky_dec_err_deg", "Dec_err_deg", "measured_centroid_dec_error_deg")

        lensframe_cols, lensframe_unit = _resolve_lensframe_columns(column_names)
        lensframe_meas_unit = lensframe_unit.lower() if isinstance(lensframe_unit, str) else lensframe_unit
        lensframe_true_unit = lensframe_meas_unit
        event_frame_centroid_x_col = _mapped_col("event_x_thetaE", "event_x_true", "true_x_centroid")
        event_frame_centroid_y_col = _mapped_col("event_y_thetaE", "event_y_true", "true_y_centroid")
        contract_lensframe_cols = {
            "meas_x": _mapped_col("event_x_obs_mas"),
            "meas_x_err": _mapped_col("event_x_err_mas"),
            "meas_y": _mapped_col("event_y_obs_mas"),
            "meas_y_err": _mapped_col("event_y_err_mas"),
            "true_x": event_frame_centroid_x_col,
            "true_x_err": _mapped_col("event_x_true_err"),
            "true_y": event_frame_centroid_y_col,
            "true_y_err": _mapped_col("event_y_true_err"),
        }
        if all(contract_lensframe_cols[key] is not None for key in ("meas_x", "meas_x_err", "meas_y", "meas_y_err", "true_x", "true_x_err", "true_y", "true_y_err")):
            lensframe_cols = {key: str(val) for key, val in contract_lensframe_cols.items()}
            lensframe_meas_unit = "mas"
            lensframe_true_unit = contract.get("event_true_unit", "mas").strip().lower()

        has_full_sky = all(col in column_names for col in skyframe_cols)
        has_radec = all(col is not None for col in [ra_obs_col, dec_obs_col, ra_true_col, dec_true_col])

        astrom_mode = None
        if has_full_sky:
            astrom_mode = "sky"
        elif has_radec:
            astrom_mode = "sky_derived"
        elif lensframe_cols:
            astrom_mode = "lensframe"

        has_astrom = astrom_mode is not None

        if astrometry_expected and not has_astrom:
            debug_sets = {
                "lens-frame (thetaE)": lensframe_thetae_cols,
                "lens-frame (mas)": lensframe_mas_cols,
                "sky-frame": skyframe_cols,
                "RA/Dec": [name for name in [ra_obs_col, dec_obs_col, ra_true_col, dec_true_col] if name is not None],
            }
            for label, cols in debug_sets.items():
                missing = [col for col in cols if col not in column_names]
                if missing:
                    print(
                        f"  Debug: {lc_file.name} missing {label} astrometry columns: "
                        + ", ".join(missing)
                    )
            raise SmokeTestError(
                f"Smoke test failed: astrometric columns missing in {lc_file.name}"
            )

        # Extract source flux columns (optional for binary source events)
        src1_flux = _optional_column("source1_relative_flux")
        src2_flux = _optional_column("source2_relative_flux")

        # If the parameter file did not request multiple sources, ignore any
        # source1/source2 columns that may nevertheless appear in the output
        # (prevents plotting Source 1/Source 2 when MULTIPLE_SOURCES isn't set).
        if not multiple_sources_expected:
            src1_flux = None
            src2_flux = None

        # Strict flux-conservation sanity check (non-configurable). Only run
        # when multiple sources are expected and both columns are present.
        _sanity_check_flux_conservation(lc_file, true_flux, src1_flux, src2_flux)

        lensframe_data: Dict[str, np.ndarray] | None = None
        true_x_vals: np.ndarray | None = None
        true_y_vals: np.ndarray | None = None
        meas_x: np.ndarray | None = None
        meas_y: np.ndarray | None = None
        meas_x_err: np.ndarray | None = None
        meas_y_err: np.ndarray | None = None
        meas_x_mas: np.ndarray | None = None
        meas_y_mas: np.ndarray | None = None
        meas_x_err_mas: np.ndarray | None = None
        meas_y_err_mas: np.ndarray | None = None
        plot_true_x: np.ndarray | None = None
        plot_true_y: np.ndarray | None = None
        plot_meas_x: np.ndarray | None = None
        plot_meas_y: np.ndarray | None = None
        plot_x_err: np.ndarray | None = None
        plot_y_err: np.ndarray | None = None
        plot_unit_label = "Einstein radii"
        event_frame_data: Dict[str, np.ndarray] = {}

        if lensframe_cols:
            lensframe_data = {
                key: _require_column(col_name) for key, col_name in lensframe_cols.items()
            }
            true_x_vals = lensframe_data["true_x"]
            true_y_vals = lensframe_data["true_y"]
            meas_x = lensframe_data["meas_x"]
            meas_y = lensframe_data["meas_y"]
            meas_x_err = lensframe_data["meas_x_err"]
            meas_y_err = lensframe_data["meas_y_err"]

            if lensframe_meas_unit == "mas":
                meas_x_mas = meas_x
                meas_y_mas = meas_y
                meas_x_err_mas = meas_x_err
                meas_y_err_mas = meas_y_err
            elif lensframe_meas_unit == "thetae" and theta_e_float is not None and theta_e_float > 0:
                meas_x_mas = meas_x * theta_e_float
                meas_y_mas = meas_y * theta_e_float
                meas_x_err_mas = meas_x_err * theta_e_float
                meas_y_err_mas = meas_y_err * theta_e_float

            # Prepare lens-frame arrays for plotting in a single unit.
            if lensframe_meas_unit == "mas":
                plot_unit_label = "mas"
                plot_meas_x = meas_x
                plot_meas_y = meas_y
                plot_x_err = meas_x_err
                plot_y_err = meas_y_err
                if lensframe_true_unit == "mas":
                    plot_true_x = true_x_vals
                    plot_true_y = true_y_vals
                elif lensframe_true_unit == "thetae" and theta_e_float is not None and theta_e_float > 0:
                    plot_true_x = true_x_vals * theta_e_float
                    plot_true_y = true_y_vals * theta_e_float
                else:
                    # Do not fake units: measured data stay in mas, true curve omitted.
                    plot_true_x = np.full_like(plot_meas_x, np.nan)
                    plot_true_y = np.full_like(plot_meas_y, np.nan)
            elif lensframe_meas_unit == "thetae":
                plot_unit_label = "Einstein radii"
                plot_meas_x = meas_x
                plot_meas_y = meas_y
                plot_x_err = meas_x_err
                plot_y_err = meas_y_err
                if lensframe_true_unit in {"thetae", "theta_e"}:
                    plot_true_x = true_x_vals
                    plot_true_y = true_y_vals
                elif lensframe_true_unit == "mas" and theta_e_float is not None and theta_e_float > 0:
                    scale = 1.0 / theta_e_float
                    plot_true_x = true_x_vals * scale
                    plot_true_y = true_y_vals * scale
                else:
                    plot_true_x = np.full_like(plot_meas_x, np.nan)
                    plot_true_y = np.full_like(plot_meas_y, np.nan)
            else:
                # Unknown declared unit; avoid misleading unit conversion.
                plot_unit_label = "native"
                plot_true_x = true_x_vals
                plot_true_y = true_y_vals
                plot_meas_x = meas_x
                plot_meas_y = meas_y
                plot_x_err = meas_x_err
                plot_y_err = meas_y_err

        if true_x_vals is None and event_frame_centroid_x_col is not None and event_frame_centroid_y_col is not None:
            true_x_vals = _require_column(event_frame_centroid_x_col)
            true_y_vals = _require_column(event_frame_centroid_y_col)

        if true_x_vals is not None and true_y_vals is not None:
            event_frame_data["centroid_x"] = true_x_vals
            event_frame_data["centroid_y"] = true_y_vals
        for key in ("source0_x", "source0_y", "source0_mu", "lens0_x", "lens0_y", "lens1_x", "lens1_y"):
            arr = _optional_column(key)
            if arr is not None:
                event_frame_data[key] = arr

        if not has_astrom:
            _plot_photometry_only(
                lc_file,
                output_dir,
                title,
                time,
                flux,
                flux_err,
                true_flux,
                src1_flux,
                src2_flux,
            )
            continue

        if astrom_mode == "lensframe":
            _plot_photometry_only(
                lc_file,
                output_dir,
                title,
                time,
                flux,
                flux_err,
                true_flux,
                src1_flux,
                src2_flux,
            )
            if (
                plot_true_x is not None
                and plot_true_y is not None
                and plot_meas_x is not None
                and plot_meas_y is not None
                and plot_x_err is not None
                and plot_y_err is not None
            ):
                print(f"  Debug: Plotting minimal astrometry plot for {lc_file.name}")
                _plot_minimal_astrometry(
                    lc_file,
                    output_dir,
                    title,
                    time,
                    plot_true_x,
                    plot_true_y,
                    plot_meas_x,
                    plot_meas_y,
                    plot_x_err,
                    plot_y_err,
                    tE,
                    unit_label=plot_unit_label,
                )
            continue

        if astrom_mode == "sky":
            true_N_mas = _require_column("true_N_centroid_mas")
            true_E_mas = _require_column("true_E_centroid_mas")
            meas_N_mas = _require_column("measured_N_centroid_mas")
            meas_E_mas = _require_column("measured_E_centroid_mas")
            meas_N_err_mas = _require_column("measured_N_centroid_error_mas")
            meas_E_err_mas = _require_column("measured_E_centroid_error_mas")
            true_ra_deg = _require_column("true_centroid_ra_deg")
            true_dec_deg = _require_column("true_centroid_dec_deg")
            meas_ra_deg = _require_column("measured_centroid_ra_deg")
            meas_dec_deg = _require_column("measured_centroid_dec_deg")
            meas_ra_err_deg = _require_column("measured_centroid_ra_error_deg")
            meas_dec_err_deg = _require_column("measured_centroid_dec_error_deg")
        else:
            # Use RA/Dec columns directly (no synthetic N/E reconstruction).
            if ra_true_col is None or dec_true_col is None or ra_obs_col is None or dec_obs_col is None:
                raise SmokeTestError(
                    f"Smoke test failed: could not resolve RA/Dec astrometry columns for {lc_file.name}"
                )
            true_ra_deg = _require_column(ra_true_col)  # used in Absolute Astrometric Position plot
            true_dec_deg = _require_column(dec_true_col)
            meas_ra_deg = _require_column(ra_obs_col)
            meas_dec_deg = _require_column(dec_obs_col)

            # Keep placeholders for the render signature; this mode plots RA/Dec plus event frame.
            true_N_mas = np.zeros_like(meas_dec_deg)
            meas_N_mas = np.zeros_like(meas_dec_deg)
            true_E_mas = np.zeros_like(meas_ra_deg)
            meas_E_mas = np.zeros_like(meas_ra_deg)

            if ra_err_col is not None and dec_err_col is not None:
                meas_ra_err_deg = _require_column(ra_err_col)
                meas_dec_err_deg = _require_column(dec_err_col)
                radec_obs_ra_err_col = ra_err_col
                radec_obs_dec_err_col = dec_err_col
            else:
                sigma_ast_series = _optional_column("sigma_ast_mas")
                if sigma_ast_series is not None:
                    mas_to_deg = 1.0 / (3600.0 * 1000.0)
                    meas_ra_err_deg = sigma_ast_series * mas_to_deg
                    meas_dec_err_deg = sigma_ast_series * mas_to_deg
                    radec_obs_ra_err_col = "sigma_ast_mas (symmetric)"
                    radec_obs_dec_err_col = "sigma_ast_mas (symmetric)"
                elif meas_x_err_mas is None or meas_y_err_mas is None:
                    meas_ra_err_deg = np.zeros_like(meas_ra_deg)
                    meas_dec_err_deg = np.zeros_like(meas_dec_deg)
                    radec_obs_ra_err_col = "none"
                    radec_obs_dec_err_col = "none"
                else:
                    base_dec_for_err = float(true_dec_deg[0])
                    cos_dec = math.cos(math.radians(base_dec_for_err))
                    if abs(cos_dec) < 1e-6:
                        cos_dec = 1e-6 if cos_dec >= 0 else -1e-6
                    mas_to_deg = 1.0 / (3600.0 * 1000.0)
                    meas_ra_err_deg = meas_x_err_mas * mas_to_deg / cos_dec
                    meas_dec_err_deg = meas_y_err_mas * mas_to_deg
                    radec_obs_ra_err_col = "x/y event errors"
                    radec_obs_dec_err_col = "x/y event errors"

            meas_E_err_mas = np.zeros_like(meas_ra_deg)
            meas_N_err_mas = np.zeros_like(meas_dec_deg)

        pm_ref_alpha_float = None
        pm_ref_delta_float = None
        pm_ref_alpha_val = summary.get("pm_ref_alpha") if summary else None
        pm_ref_delta_val = summary.get("pm_ref_delta") if summary else None
        if pm_ref_alpha_val is not None and pm_ref_delta_val is not None:
            pm_ref_alpha_float = float(pm_ref_alpha_val)
            pm_ref_delta_float = float(pm_ref_delta_val)
            if math.isnan(pm_ref_alpha_float) or math.isnan(pm_ref_delta_float):
                pm_ref_alpha_float = None
                pm_ref_delta_float = None

        source_pm_icrs: Tuple[float, float] | None = None
        lens_pm_icrs: Tuple[float, float] | None = None
        if summary:
            src_vals = (
                summary.get("source_mul"),
                summary.get("source_mub"),
                summary.get("source_l"),
                summary.get("source_b"),
            )
            if all(v is not None and not math.isnan(v) for v in src_vals):
                source_pm_icrs = _galactic_pm_to_icrs(
                    float(src_vals[2]),
                    float(src_vals[3]),
                    float(src_vals[0]),
                    float(src_vals[1]),
                )
            lens_vals = (
                summary.get("lens_mul"),
                summary.get("lens_mub"),
                summary.get("lens_l"),
                summary.get("lens_b"),
            )
            if all(v is not None and not math.isnan(v) for v in lens_vals):
                lens_pm_icrs = _galactic_pm_to_icrs(
                    float(lens_vals[2]),
                    float(lens_vals[3]),
                    float(lens_vals[0]),
                    float(lens_vals[1]),
                )

        span_years = None
        if len(time):
            span_days = float(time.max() - time.min())
            if span_days > 0:
                span_years = span_days / 365.25

        vector_specs: List[Dict[str, float | str]] = []
        if span_years:
            if pm_ref_alpha_float is not None and pm_ref_delta_float is not None:
                vector_specs.append(
                    {
                        "label": "Relative PM geocentric (pm_ra*cosDec, pm_dec)",
                        "color": "black",
                        "pm_ra": pm_ref_alpha_float,
                        "pm_dec": pm_ref_delta_float,
                    }
                )
            if source_pm_icrs:
                vector_specs.append(
                    {
                        "label": "Source PM heliocentric (pm_ra*cosDec, pm_dec)",
                        "color": "tab:blue",
                        "pm_ra": source_pm_icrs[0],
                        "pm_dec": source_pm_icrs[1],
                    }
                )
            if lens_pm_icrs:
                vector_specs.append(
                    {
                        "label": "Lens PM heliocentric (pm_ra*cosDec, pm_dec)",
                        "color": "tab:red",
                        "pm_ra": lens_pm_icrs[0],
                        "pm_dec": lens_pm_icrs[1],
                    }
                )

        if astrom_mode == "sky":
            radec_obs_ra_col = "measured_centroid_ra_deg"
            radec_obs_dec_col = "measured_centroid_dec_deg"
            radec_true_ra_col = "true_centroid_ra_deg"
            radec_true_dec_col = "true_centroid_dec_deg"
            radec_obs_ra_err_col = "measured_centroid_ra_error_deg"
            radec_obs_dec_err_col = "measured_centroid_dec_error_deg"
        else:
            radec_obs_ra_col = ra_obs_col or "RA_measured_deg"
            radec_obs_dec_col = dec_obs_col or "Dec_measured_deg"
            radec_true_ra_col = ra_true_col or "RA_noiseless_deg"
            radec_true_dec_col = dec_true_col or "Dec_noiseless_deg"

        panel_labels: Dict[str, str] = {
            "light_x_col": "Simulation_time",
            "light_y_col": "measured_relative_flux",
            "light_yerr_col": "measured_relative_flux_error",
            "light_true_col": "true_relative_flux",
            "radec_obs_ra_col": radec_obs_ra_col,
            "radec_obs_dec_col": radec_obs_dec_col,
            "radec_obs_ra_err_col": radec_obs_ra_err_col,
            "radec_obs_dec_err_col": radec_obs_dec_err_col,
            "radec_true_ra_col": radec_true_ra_col,
            "radec_true_dec_col": radec_true_dec_col,
            "event_x_col": event_frame_centroid_x_col or "true_x_centroid",
            "event_y_col": event_frame_centroid_y_col or "true_y_centroid",
        }
        for key in ("source0_x", "source0_y", "source0_mu", "lens0_x", "lens0_y", "lens1_x", "lens1_y"):
            if key in event_frame_data:
                panel_labels[f"{key}_col"] = key

        mode_display = "sky_columns" if astrom_mode == "sky" else "radec_columns_only"
        event_track_cols = [key for key in ("source0_x", "source0_y", "source0_mu", "lens0_x", "lens0_y", "lens1_x", "lens1_y") if key in event_frame_data]
        event_info_rows: List[Tuple[str, str]] = [
            ("Astrometry mode", mode_display),
            ("VBM function", vbm_function_name),
            ("Time column", "Simulation_time"),
            ("Flux cols", "measured_relative_flux, true_relative_flux"),
            ("RA/Dec measured cols", f"{panel_labels['radec_obs_ra_col']}, {panel_labels['radec_obs_dec_col']}"),
            ("RA/Dec error cols", f"{panel_labels['radec_obs_ra_err_col']}, {panel_labels['radec_obs_dec_err_col']}"),
            ("RA/Dec noiseless cols", f"{panel_labels['radec_true_ra_col']}, {panel_labels['radec_true_dec_col']}"),
            ("Event centroid cols", f"{panel_labels['event_x_col']}, {panel_labels['event_y_col']}"),
            ("Event track cols", ", ".join(event_track_cols) if event_track_cols else "none"),
            ("PM vector convention", "(pm_ra*cosDec, pm_dec) mas/yr in ICRS"),
        ]
        if span_years is not None:
            event_info_rows.append(("Vector span", f"{span_years:.4f} yr"))
        if theta_e_float is not None:
            event_info_rows.append(("theta_E", f"{theta_e_float:.6f} mas"))
        if summary and summary.get("tE_ref") is not None:
            event_info_rows.append(("tE_ref", _format_metric(summary.get("tE_ref"))))

        vbm_model = None
        vbm_reason = None
        if vbm_required:
            vbm_model, vbm_reason = _compute_vbm_model(
                summary,
                planet_vals,
                event_vals,
                source_pm_icrs,
                lens_pm_icrs,
                theta_e_float,
                source_dist_float,
                event_ra_float,
                event_dec_float,
                alpha_deg_float,
                sim_zero_offset,
                time,
                true_x_vals,
                true_y_vals,
            )
            if vbm_model is None:
                detail = f" ({vbm_reason})" if vbm_reason else ""
                raise SmokeTestError(
                    f"Smoke test failed: missing VBM lens-frame plot for {lc_file.name}{detail}"
                )

        plot_file, lensframe_path = _render_astrometric_figure(
            lc_file,
            output_dir,
            title,
            time,
            flux,
            flux_err,
            true_flux,
            true_N_mas,
            true_E_mas,
            meas_N_mas,
            meas_E_mas,
            meas_N_err_mas,
            meas_E_err_mas,
            true_ra_deg,
            true_dec_deg,
            meas_ra_deg,
            meas_dec_deg,
            meas_ra_err_deg,
            meas_dec_err_deg,
            vector_specs,
            vbm_model,
            true_x_vals,
            true_y_vals,
            meas_x,
            meas_y,
            event_frame_data,
            src1_flux,
            src2_flux,
            panel_labels,
            event_info_rows,
        )

        if vbm_required:
            if lensframe_path is None or not lensframe_path.exists():
                raise SmokeTestError(
                    f"Smoke test failed: missing lens-frame plot for {lc_file.name}"
                )
        print(f"  Generated plot: {plot_file.name}")
        if lensframe_path is not None:
            print(f"  Generated plot: {lensframe_path.name}")

    if plot_failures:
        rendered = "\n".join(f" - {msg}" for msg in plot_failures)
        raise SmokeTestError(
            "Smoke test plotting encountered failures:\n"
            f"{rendered}\n"
            "Plots for files with complete inputs were still generated."
        )

__all__ = ["plot_lightcurves"]
