"""Direct BAGLE forward-model diagnostics for explicit 1s1l Gulls events."""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Sequence

import numpy as np
import pandas as pd

from .bagle_fit_sanity import _parse_header_keyvals
from .constants import REPO_ROOT
from .errors import SmokeTestError

MAS_PER_ARCSEC = 1000.0


@dataclass
class BagleForwardSummary:
    run_dir: Path
    lc_file: Path
    event_id: int
    bagle_model_name: str
    astrometry_rms_mas: float
    photometry_rms_relative_flux: float
    photometry_max_abs_relative_flux_diff: float
    summary_json_path: Path
    plot_path: Path
    warnings: List[str] = field(default_factory=list)


def _find_matching_lc(run_dir: Path, event_id: int) -> Path:
    candidates = sorted(run_dir.glob(f"*_{event_id}.all.lc"))
    if not candidates:
        raise SmokeTestError(
            f"BAGLE forward sanity: could not find .all.lc for EventID={event_id} in {run_dir}"
        )
    return candidates[0]


def _parse_scalar_header(lc_file: Path, key: str) -> float:
    prefix = f"#{key}:"
    for raw in lc_file.read_text(encoding="utf-8").splitlines():
        if raw.startswith(prefix):
            parts = raw.split(":", 1)[1].split()
            if not parts:
                raise SmokeTestError(
                    f"BAGLE forward sanity: header {prefix} is empty in {lc_file.name}"
                )
            return float(parts[0])
    raise SmokeTestError(
        f"BAGLE forward sanity: missing header {prefix} in {lc_file.name}"
    )


def _header_magnitudes(lc_file: Path) -> tuple[float, float]:
    """Return the published source and lens magnitudes from a Gulls lightcurve header."""
    return _parse_scalar_header(lc_file, "Obssrcmag"), _parse_scalar_header(lc_file, "Obslensmag")


def _infer_simulation_zero_jd(lc_file: Path, df: pd.DataFrame) -> float:
    """Recover SIMULATION_ZERO_TIME from the published BJD and Simulation_time columns."""
    if "BJD" not in df.columns:
        raise SmokeTestError(
            f"BAGLE forward sanity: {lc_file.name} is missing BJD, so "
            "SIMULATION_ZERO_TIME cannot be inferred."
        )
    if "Simulation_time" not in df.columns:
        raise SmokeTestError(
            f"BAGLE forward sanity: {lc_file.name} is missing Simulation_time."
        )

    bjd = df["BJD"].to_numpy(dtype=float, copy=False)
    sim_time = df["Simulation_time"].to_numpy(dtype=float, copy=False)
    if not np.all(np.isfinite(bjd)) or not np.all(np.isfinite(sim_time)):
        raise SmokeTestError(
            f"BAGLE forward sanity: {lc_file.name} has non-finite BJD or Simulation_time values."
        )

    sim_zero = bjd - sim_time
    sim_zero_ref = float(sim_zero[0])
    max_diff_days = float(np.max(np.abs(sim_zero - sim_zero_ref)))
    max_diff_sec = max_diff_days * 86400.0
    if max_diff_sec > 1.0e-4:
        raise SmokeTestError(
            f"BAGLE forward sanity: inferred SIMULATION_ZERO_TIME from {lc_file.name} is inconsistent "
            f"across epochs (max spread {max_diff_sec:.6e} s)."
        )
    return sim_zero_ref


def _load_out_table(run_dir: Path) -> pd.DataFrame:
    out_candidates = sorted(run_dir.glob("*.out"))
    if not out_candidates:
        raise SmokeTestError(f"BAGLE forward sanity: no .out file found in {run_dir}")
    return pd.read_csv(out_candidates[0], sep=r"\s+")


def _select_event_row(df: pd.DataFrame, event_id: int | None) -> pd.Series:
    if event_id is None:
        matches = df[(df["NSource"] == 1) & (df["NLens"] == 1) & (df["NPlanets"] == 0)]
        if matches.empty:
            raise SmokeTestError(
                "BAGLE forward sanity: no explicit 1s1l event found in the .out file."
            )
        return matches.iloc[0]

    matches = df[df["EventID"] == event_id]
    if matches.empty:
        raise SmokeTestError(f"BAGLE forward sanity: EventID={event_id} not found in .out")
    return matches.iloc[0]


def _plot_forward_comparison(
    output_path: Path,
    *,
    event_id: int,
    bagle_model_name: str,
    sim_time_days: np.ndarray,
    t0_event_days: float,
    true_flux: np.ndarray,
    bagle_flux: np.ndarray,
    x_true_arcsec: np.ndarray,
    y_true_arcsec: np.ndarray,
    bagle_ast_arcsec: np.ndarray,
    astrometry_rms_mas: float,
    photometry_rms_relative_flux: float,
    photometry_max_abs_relative_flux_diff: float,
) -> None:
    """Plot direct BAGLE-vs-Gulls noiseless photometry and astrometry for a 1s1l event."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - optional dependency guard
        raise SmokeTestError(
            "BAGLE forward sanity: matplotlib is required to generate forward-model plots."
        ) from exc

    sim_time_days = np.asarray(sim_time_days, dtype=float)
    true_flux = np.asarray(true_flux, dtype=float)
    bagle_flux = np.asarray(bagle_flux, dtype=float)
    x_true_arcsec = np.asarray(x_true_arcsec, dtype=float)
    y_true_arcsec = np.asarray(y_true_arcsec, dtype=float)
    bagle_ast_arcsec = np.asarray(bagle_ast_arcsec, dtype=float)

    fig, axes = plt.subplots(3, 1, figsize=(11, 12), constrained_layout=True)

    axes[0].plot(
        sim_time_days,
        true_flux,
        ".",
        ms=2.0,
        alpha=0.40,
        color="tab:blue",
        label="Gulls noiseless",
        zorder=6,
    )
    axes[0].plot(
        sim_time_days,
        bagle_flux,
        "--",
        lw=2.0,
        color="tab:red",
        label="BAGLE forward",
        zorder=20,
    )
    axes[0].axvline(t0_event_days, color="0.55", lw=1.0, ls=":", zorder=1)
    axes[0].set_ylabel("Relative Flux")
    axes[0].set_title(
        (
            f"BAGLE Forward vs Gulls Noiseless: event {event_id} ({bagle_model_name})\n"
            f"astrometry RMS={astrometry_rms_mas:.4f} mas, "
            f"photometry RMS={photometry_rms_relative_flux:.5f}, "
            f"max |dflux|={photometry_max_abs_relative_flux_diff:.5f}"
        )
    )
    axes[0].legend(loc="best", fontsize=9)

    axes[1].plot(
        sim_time_days,
        x_true_arcsec * MAS_PER_ARCSEC,
        ".",
        ms=2.0,
        alpha=0.40,
        color="tab:olive",
        label="E Gulls",
        zorder=6,
    )
    axes[1].plot(
        sim_time_days,
        y_true_arcsec * MAS_PER_ARCSEC,
        ".",
        ms=2.0,
        alpha=0.40,
        color="tab:brown",
        label="N Gulls",
        zorder=6,
    )
    axes[1].plot(
        sim_time_days,
        bagle_ast_arcsec[:, 0] * MAS_PER_ARCSEC,
        "--",
        lw=2.0,
        color="tab:green",
        label="E BAGLE",
        zorder=20,
    )
    axes[1].plot(
        sim_time_days,
        bagle_ast_arcsec[:, 1] * MAS_PER_ARCSEC,
        "--",
        lw=2.0,
        color="tab:purple",
        label="N BAGLE",
        zorder=20,
    )
    axes[1].axvline(t0_event_days, color="0.55", lw=1.0, ls=":", zorder=1)
    axes[1].set_ylabel("Centroid (mas)")
    axes[1].legend(loc="best", fontsize=9)

    axes[2].plot(
        x_true_arcsec * MAS_PER_ARCSEC,
        y_true_arcsec * MAS_PER_ARCSEC,
        ".",
        ms=2.0,
        alpha=0.40,
        color="tab:orange",
        label="Gulls noiseless",
        zorder=6,
    )
    axes[2].plot(
        bagle_ast_arcsec[:, 0] * MAS_PER_ARCSEC,
        bagle_ast_arcsec[:, 1] * MAS_PER_ARCSEC,
        "--",
        lw=2.0,
        color="tab:red",
        label="BAGLE forward",
        zorder=20,
    )
    axes[2].set_xlabel("E (mas)")
    axes[2].set_ylabel("N (mas)")
    axes[2].set_aspect("equal", adjustable="box")
    axes[2].legend(loc="best", fontsize=9)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def run_bagle_forward_1s1l_sanity(
    run_dir: Path,
    *,
    event_id: int | None = None,
    obs_location: str = "earth",
) -> BagleForwardSummary:
    """Compare a direct BAGLE geocentric-projected forward model against Gulls noiseless 1s1l outputs."""
    try:
        from bagle import model
    except Exception as exc:  # pragma: no cover - environment guard
        raise SmokeTestError(
            "BAGLE forward sanity requires the bagle package to be importable."
        ) from exc

    run_dir = run_dir.resolve()
    out_df = _load_out_table(run_dir)
    row = _select_event_row(out_df, event_id)

    evt = int(float(row["EventID"]))
    if int(float(row["NSource"])) != 1 or int(float(row["NLens"])) != 1 or int(float(row["NPlanets"])) != 0:
        raise SmokeTestError(
            f"BAGLE forward sanity: EventID={evt} is not 1s1l "
            f"(NSource={row['NSource']}, NLens={row['NLens']}, NPlanets={row['NPlanets']})."
        )

    lc_file = _find_matching_lc(run_dir, evt)
    df = pd.read_csv(lc_file, sep=r"\s+", comment="#")
    if df.empty:
        raise SmokeTestError(f"BAGLE forward sanity: {lc_file.name} contains no data rows.")

    if "Simulation_time" not in df.columns:
        raise SmokeTestError(
            f"BAGLE forward sanity: {lc_file.name} is missing Simulation_time."
        )
    for required_col in ("true_relative_flux", "true_RA_deg", "true_Dec_deg"):
        if required_col not in df.columns:
            raise SmokeTestError(
                f"BAGLE forward sanity: {lc_file.name} is missing required column {required_col}."
            )

    sim_zero_jd = _infer_simulation_zero_jd(lc_file, df)

    ra_deg = float(row["ra_deg"])
    dec_deg = float(row["dec_deg"])
    cos_dec = math.cos(math.radians(dec_deg))
    if abs(cos_dec) < 1.0e-8:
        raise SmokeTestError(
            f"BAGLE forward sanity: cos(dec) too small at event pointing (dec={dec_deg:.8f} deg)."
        )

    sim_time = df["Simulation_time"].to_numpy(dtype=float, copy=False)
    t_mjd = sim_zero_jd + sim_time - 2400000.5
    true_flux = df["true_relative_flux"].to_numpy(dtype=float, copy=False)
    dra_deg = (df["true_RA_deg"].to_numpy(dtype=float, copy=False) - ra_deg + 180.0) % 360.0 - 180.0
    x_true_arcsec = dra_deg * cos_dec * 3600.0
    y_true_arcsec = (df["true_Dec_deg"].to_numpy(dtype=float, copy=False) - dec_deg) * 3600.0

    thetaE_mas = float(row["thetaE"])
    if not math.isfinite(thetaE_mas) or thetaE_mas <= 0.0:
        raise SmokeTestError(f"BAGLE forward sanity: invalid thetaE for EventID={evt}.")
    dS_kpc = float(row["Source_Dist"])
    if not math.isfinite(dS_kpc) or dS_kpc <= 0.0:
        raise SmokeTestError(f"BAGLE forward sanity: invalid Source_Dist for EventID={evt}.")
    obs_fs = float(row["Obs_0_fs"])
    if not math.isfinite(obs_fs) or not (0.0 < obs_fs <= 1.0):
        raise SmokeTestError(f"BAGLE forward sanity: invalid Obs_0_fs for EventID={evt}.")

    source_mag, _ = _header_magnitudes(lc_file)
    mag_base = source_mag + 2.5 * math.log10(obs_fs)
    t0_event_days = float(row["t0lens1"])
    t0par_event_days = float(row["tref"])
    t0_mjd = sim_zero_jd + t0_event_days - 2400000.5
    t0par_mjd = sim_zero_jd + t0par_event_days - 2400000.5
    reference_idx = int(np.argmin(np.abs(t_mjd - t0_mjd)))
    xS0_guess = float(x_true_arcsec[reference_idx])
    yS0_guess = float(y_true_arcsec[reference_idx])

    bagle_model_name = "pspl_param4_geoproj"
    piE_amp = float(row["piE"])
    mu_rel_ref_E = float(row["murel_ref_alpha"])
    mu_rel_ref_N = float(row["murel_ref_delta"])
    mu_rel_ref_amp = math.hypot(mu_rel_ref_E, mu_rel_ref_N)
    if not math.isfinite(mu_rel_ref_amp) or mu_rel_ref_amp <= 0.0:
        raise SmokeTestError(
            f"BAGLE forward sanity: invalid reference-frame equatorial mu_rel for EventID={evt}."
        )
    # BAGLE expects geocentric-projected piE in on-sky East/North coordinates.
    # Reuse the already-published Gulls reference-frame equatorial mu_rel direction
    # rather than the ecliptic piEE/piEN pair.
    piE_E_geotr = piE_amp * mu_rel_ref_E / mu_rel_ref_amp
    piE_N_geotr = piE_amp * mu_rel_ref_N / mu_rel_ref_amp

    # GULLS defines u0_hat as the 90-degree CCW rotation of tau_hat in the
    # ecliptic event frame.  BAGLE's geoproj model uses coord_in='tb'
    # internally, where u0>0 corresponds to the CW rotation of tau_hat —
    # the opposite handedness.  Negating u0 maps from GULLS convention to
    # BAGLE's tau-beta convention.
    u0_for_bagle = -float(row["u0lens1"])

    muS_E = float(row["mu_source_helio_alpha"])
    muS_N = float(row["mu_source_helio_delta"])
    bagle_model = model.PSPL_PhotAstrom_Par_Param4_geoproj(
        t0_mjd,
        u0_for_bagle,
        float(row["tE_ref"]),
        thetaE_mas,
        1.0 / dS_kpc,
        piE_E_geotr,
        piE_N_geotr,
        xS0_guess,
        yS0_guess,
        muS_E,
        muS_N,
        [obs_fs],
        [mag_base],
        t0par_mjd,
        raL=ra_deg,
        decL=dec_deg,
        obsLocation=obs_location,
    )

    warnings: List[str] = [
        (
            f"{lc_file.name}: BAGLE forward model uses geocentric-projected photometric tuple "
            "from Gulls: t0lens1, tE_ref, and piE projected into equatorial E/N "
            "from the published murel_ref_alpha/delta direction, with t0par=tref.  "
            "u0 is negated (GULLS CCW u0_hat convention -> BAGLE tb CW convention)."
        ),
        (
            f"{lc_file.name}: astrometric source position and proper motion remain in the SSB/heliocentric "
            "form BAGLE expects for the geoproj class."
        ),
    ]
    truth_mapping: Dict[str, Any] = {
        "t0_geotr_mjd_from_gulls_t0lens1": float(t0_mjd),
        "t0par_mjd_from_gulls_tref": float(t0par_mjd),
        "u0_amp_geotr_from_gulls_u0lens1": float(row["u0lens1"]),
        "u0_amp_geotr_negated_for_bagle_tb": float(u0_for_bagle),
        "tE_geotr_days_from_gulls_tE_ref": float(row["tE_ref"]),
        "thetaE_mas_from_gulls_thetaE": float(thetaE_mas),
        "piS_mas_from_gulls_source_dist_kpc": float(1.0 / dS_kpc),
        "piE_geotr_lens_minus_source_from_gulls_reference_eq": {
            "E": float(piE_E_geotr),
            "N": float(piE_N_geotr),
        },
        "piE_ref_ecliptic_from_gulls_raw": {
            "E_lambda": float(row["piEE"]),
            "N_beta": float(row["piEN"]),
        },
        "muRel_ref_eq": {"E": float(mu_rel_ref_E), "N": float(mu_rel_ref_N)},
        "muS_helio_eq": {"E": float(muS_E), "N": float(muS_N)},
        "blend_source_flux_fraction": float(obs_fs),
        "source_mag": float(source_mag),
        "mag_base_from_gulls_source_mag_and_fs": float(mag_base),
        "xS0_arcsec_initial": {"E": float(xS0_guess), "N": float(yS0_guess)},
    }

    ast_guess = np.asarray(bagle_model.get_astrometry(t_mjd), dtype=float)
    translation_arcsec = np.array(
        [
            x_true_arcsec[reference_idx] - ast_guess[reference_idx, 0],
            y_true_arcsec[reference_idx] - ast_guess[reference_idx, 1],
        ],
        dtype=float,
    )
    bagle_ast = ast_guess + translation_arcsec
    bagle_mag = np.asarray(bagle_model.get_photometry(t_mjd), dtype=float)
    mag_base = float(bagle_model.mag_base[0])
    bagle_flux_rel = 10.0 ** (-0.4 * (bagle_mag - mag_base))

    astrometry_rms_mas = float(
        np.sqrt(
            np.mean(
                ((bagle_ast[:, 0] - x_true_arcsec) * MAS_PER_ARCSEC) ** 2
                + ((bagle_ast[:, 1] - y_true_arcsec) * MAS_PER_ARCSEC) ** 2
            )
        )
    )
    photometry_rms_relative_flux = float(
        np.sqrt(np.mean((bagle_flux_rel - true_flux) ** 2))
    )
    photometry_max_abs_relative_flux_diff = float(
        np.max(np.abs(bagle_flux_rel - true_flux))
    )

    fit_dir = run_dir / "bagle_forward_sanity"
    fit_dir.mkdir(parents=True, exist_ok=True)
    model_suffix = bagle_model_name.lower()
    summary_json_path = fit_dir / f"event_{evt:06d}_{model_suffix}_summary.json"
    plot_path = fit_dir / f"event_{evt:06d}_{model_suffix}_forward_model.png"

    _plot_forward_comparison(
        plot_path,
        event_id=evt,
        bagle_model_name=bagle_model_name,
        sim_time_days=sim_time,
        t0_event_days=t0_event_days,
        true_flux=true_flux,
        bagle_flux=bagle_flux_rel,
        x_true_arcsec=x_true_arcsec,
        y_true_arcsec=y_true_arcsec,
        bagle_ast_arcsec=bagle_ast,
        astrometry_rms_mas=astrometry_rms_mas,
        photometry_rms_relative_flux=photometry_rms_relative_flux,
        photometry_max_abs_relative_flux_diff=photometry_max_abs_relative_flux_diff,
    )

    payload: Dict[str, Any] = {
        "event_id": evt,
        "lc_file": str(lc_file),
        "obs_location": obs_location,
        "plot_path": str(plot_path),
        "bagle_model": bagle_model_name,
        "simulation_zero_jd_inferred": float(sim_zero_jd),
        "truth_mapping": truth_mapping,
        "astrometry": {
            "rms_mas": astrometry_rms_mas,
            "x_translation_matched_at_reference_epoch_arcsec": float(translation_arcsec[0]),
            "y_translation_matched_at_reference_epoch_arcsec": float(translation_arcsec[1]),
            "reference_epoch_simulation_day": float(t0_event_days),
        },
        "photometry": {
            "rms_relative_flux": photometry_rms_relative_flux,
            "max_abs_relative_flux_diff": photometry_max_abs_relative_flux_diff,
            "peak_simulation_day_bagle": float(sim_time[int(np.argmax(bagle_flux_rel))]),
            "peak_simulation_day_gulls": float(sim_time[int(np.argmax(true_flux))]),
            "peak_offset_days": float(
                sim_time[int(np.argmax(bagle_flux_rel))] - sim_time[int(np.argmax(true_flux))]
            ),
        },
        "warnings": warnings,
    }
    summary_json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    return BagleForwardSummary(
        run_dir=run_dir,
        lc_file=lc_file,
        event_id=evt,
        bagle_model_name=bagle_model_name,
        astrometry_rms_mas=astrometry_rms_mas,
        photometry_rms_relative_flux=photometry_rms_relative_flux,
        photometry_max_abs_relative_flux_diff=photometry_max_abs_relative_flux_diff,
        summary_json_path=summary_json_path,
        plot_path=plot_path,
        warnings=warnings,
    )
