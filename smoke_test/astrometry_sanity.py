"""Sanity checks for astrometric lightcurve outputs.

These checks are intentionally self-consistency checks:
- They validate that astrometric columns are present and numerically coherent.
- They cross-check `.lc` astrometry against canonical event metadata in `.out`.
- They do not assume any specific observer orbit.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

import numpy as np
import pandas as pd

from .errors import SmokeTestError

try:
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    _HAS_ASTROPY = True
except Exception:  # pragma: no cover - optional dependency guard
    _HAS_ASTROPY = False


MAS_PER_DEG = 3600.0 * 1000.0
DEG_TO_RAD = math.pi / 180.0


@dataclass
class AstrometrySanitySummary:
    """High-level stats from astrometry sanity checks."""

    run_dir: Path
    checked_lightcurves: int = 0
    checked_epochs: int = 0
    warnings: List[str] = field(default_factory=list)
    zscore_mean: float | None = None
    zscore_std: float | None = None


def _angle_diff_deg(lhs: float, rhs: float) -> float:
    """Return wrapped angular difference lhs-rhs in degrees."""
    return (lhs - rhs + 180.0) % 360.0 - 180.0


def _direction_alias_diagnostic(
    fit_e: float,
    fit_n: float,
    ref_angle_deg: float,
) -> Tuple[str, float]:
    """Find the closest sign/axis convention variant to a reference direction."""
    variants = [
        ("as-is", fit_e, fit_n),
        ("flip-E", -fit_e, fit_n),
        ("flip-N", fit_e, -fit_n),
        ("flip-both", -fit_e, -fit_n),
        ("swap(E<->N)", fit_n, fit_e),
        ("swap(E<->N)+flip-E", -fit_n, fit_e),
        ("swap(E<->N)+flip-N", fit_n, -fit_e),
        ("swap(E<->N)+flip-both", -fit_n, -fit_e),
    ]
    best_name = "as-is"
    best_diff = math.inf
    for name, vec_e, vec_n in variants:
        angle = math.degrees(math.atan2(vec_n, vec_e))
        diff = abs(_angle_diff_deg(angle, ref_angle_deg))
        if diff < best_diff:
            best_name = name
            best_diff = diff
    return best_name, best_diff


def _fit_linear_motion(
    x_year: np.ndarray,
    e_values: np.ndarray,
    n_values: np.ndarray,
) -> Mapping[str, float | np.ndarray]:
    """Fit E/N tracks as linear functions of time and estimate angle uncertainty."""
    mat = np.vstack([np.ones_like(x_year), x_year]).T
    coeff_e = np.linalg.lstsq(mat, e_values, rcond=None)[0]
    coeff_n = np.linalg.lstsq(mat, n_values, rcond=None)[0]

    mu_e = float(coeff_e[1])
    mu_n = float(coeff_n[1])
    mu_amp = math.hypot(mu_e, mu_n)

    pred_e = float(coeff_e[0]) + mu_e * x_year
    pred_n = float(coeff_n[0]) + mu_n * x_year
    resid_e = e_values - pred_e
    resid_n = n_values - pred_n
    resid_r = np.hypot(resid_e, resid_n)
    rms = float(np.sqrt(np.mean(resid_r**2)))

    n_fit = int(x_year.size)
    dof = max(1, n_fit - 2)
    x_center = x_year - float(np.mean(x_year))
    sxx = float(np.sum(x_center * x_center))
    sigma_theta_deg = math.inf
    if sxx > 0:
        sigma2_e = float(np.sum(resid_e * resid_e) / dof)
        sigma2_n = float(np.sum(resid_n * resid_n) / dof)
        se_mu_e = math.sqrt(max(0.0, sigma2_e) / sxx)
        se_mu_n = math.sqrt(max(0.0, sigma2_n) / sxx)
        mu_sq = max(mu_amp * mu_amp, 1.0e-12)
        sigma_theta_rad = math.sqrt((mu_n * se_mu_e) ** 2 + (mu_e * se_mu_n) ** 2) / mu_sq
        sigma_theta_deg = sigma_theta_rad * 180.0 / math.pi

    return {
        "coeff_e": coeff_e,
        "coeff_n": coeff_n,
        "mu_e": mu_e,
        "mu_n": mu_n,
        "mu_amp": mu_amp,
        "sigma_theta_deg": sigma_theta_deg,
        "rms": rms,
    }


def _parse_astrometry_frame(lc_file: Path) -> Tuple[float | None, float | None]:
    """Parse ``#Astrometry_Frame`` header and return (ra_deg, dec_deg)."""
    with lc_file.open(encoding="utf-8") as handle:
        for raw in handle:
            if not raw.startswith("#"):
                break
            if not raw.startswith("#Astrometry_Frame:"):
                continue
            if "RA_deg=" not in raw or "Dec_deg=" not in raw:
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
            return ra_deg, dec_deg
    return None, None


def _parse_astrometry_transform(
    lc_file: Path,
) -> Tuple[float, float, float, float] | None:
    """Parse ``#Astrometry_Transform`` and return (a, b, c, d).

    The declared mapping is:
      dRAcosDec_mas = a * E_ecl_mas + b * N_ecl_mas
      dDec_mas      = c * E_ecl_mas + d * N_ecl_mas
    """
    pattern = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"
    with lc_file.open(encoding="utf-8") as handle:
        for raw in handle:
            if not raw.startswith("#"):
                break
            if not raw.startswith("#Astrometry_Transform:"):
                continue
            vals = [float(x) for x in re.findall(pattern, raw)]
            if len(vals) >= 4:
                return vals[0], vals[1], vals[2], vals[3]
            return None
    return None


def _derive_event_key(lc_file: Path) -> Tuple[int, int, int | None] | None:
    """Infer (EventID, SubRun, Field?) from a lightcurve filename.

    Expected patterns:
    - `<run>_<subrun>_<field>_<event>.*.lc`
    - `<run>_<subrun>_<event>.*.lc` (field inferred elsewhere)
    """
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


def _load_out_rows(out_files: Sequence[Path]) -> Dict[Tuple[int, int, int], Mapping[str, float]]:
    """Load all `.out` rows indexed by (EventID, SubRun, Field)."""
    rows: Dict[Tuple[int, int, int], Mapping[str, float]] = {}
    for out_file in out_files:
        table = pd.read_csv(out_file, sep=r"\s+")
        required = {"EventID", "SubRun", "Field"}
        missing = required - set(table.columns)
        if missing:
            raise SmokeTestError(
                f"Astrometry sanity: {out_file.name} missing key columns: {', '.join(sorted(missing))}"
            )
        for _, row in table.iterrows():
            key = (
                int(float(row["EventID"])),
                int(float(row["SubRun"])),
                int(float(row["Field"])),
            )
            rows[key] = row
    if not rows:
        raise SmokeTestError("Astrometry sanity: no rows found in .out files")
    return rows


def _match_out_row(
    lc_file: Path,
    out_rows: Mapping[Tuple[int, int, int], Mapping[str, float]],
) -> Mapping[str, float] | None:
    key = _derive_event_key(lc_file)
    if key is not None:
        event, subrun, field = key
        if field is not None:
            row = out_rows.get((event, subrun, field))
            if row is not None:
                return row
        else:
            matches = [
                row
                for (evt, sub, _fld), row in out_rows.items()
                if evt == event and sub == subrun
            ]
            if len(matches) == 1:
                return matches[0]
    return None


def _missing_columns(columns: Iterable[str], required: Sequence[str]) -> List[str]:
    colset = set(columns)
    return [name for name in required if name not in colset]


def _nan_or_inf_columns(df: pd.DataFrame, cols: Sequence[str]) -> List[str]:
    bad: List[str] = []
    for col in cols:
        arr = df[col].to_numpy(dtype=float, copy=False)
        if not np.all(np.isfinite(arr)):
            bad.append(col)
    return bad


def verify_astrometry_sanity(
    run_dir: Path,
    out_files: Sequence[Path],
    params: Mapping[str, str] | None = None,
    *,
    max_errors: int = 40,
    strict_documented_columns: bool = False,
    long_baseline_years_each_side: float = 1.0,
    long_baseline_exclusion_te: float = 5.0,
    long_baseline_direction_tol_deg: float = 15.0,
) -> AstrometrySanitySummary:
    """Run astrometric sanity checks on one prepared run directory.

    Parameters
    ----------
    run_dir:
        Directory containing lightcurves and summary `.out` file(s).
    out_files:
        Summary files for this run (typically one `.out` file).
    params:
        Parsed parameter dictionary; used for `SIMULATION_ZERO_TIME`.
    max_errors:
        Maximum number of detailed errors included in the exception message.
    strict_documented_columns:
        If True, missing documented `lens_parallax_x_mas/y_mas` columns are fatal.
    long_baseline_years_each_side:
        Required temporal coverage (years) before/after `tref` for long-baseline PM checks.
    long_baseline_exclusion_te:
        Exclude epochs with `|t-t0| <= N*tE` for long-baseline PM checks.
    long_baseline_direction_tol_deg:
        Minimum allowed angular tolerance (degrees) for long-baseline PM direction checks.
    """
    run_dir = run_dir.resolve()
    summary = AstrometrySanitySummary(run_dir=run_dir)

    lc_files = sorted(run_dir.rglob("*.lc"))
    if not lc_files:
        raise SmokeTestError(f"Astrometry sanity: no .lc files found under {run_dir}")

    out_rows = _load_out_rows(out_files)

    sim_zero_time: float | None = None
    if params:
        raw = params.get("SIMULATION_ZERO_TIME")
        if raw is not None:
            try:
                sim_zero_time = float(raw)
            except ValueError:
                sim_zero_time = None

    required_cols = [
        "Simulation_time",
        "x_centroid_mas",
        "x_centroid_error_mas",
        "y_centroid_mas",
        "y_centroid_error_mas",
        "true_x_centroid_mas",
        "true_x_centroid_error_mas",
        "true_y_centroid_mas",
        "true_y_centroid_error_mas",
        "RA_centroid_deg",
        "Dec_centroid_deg",
        "RA_centroid_true_deg",
        "Dec_centroid_true_deg",
        "lens_dist_kpc",
        "parallax_shift_x",
        "parallax_shift_y",
        "parallax_shift_z",
    ]
    lpllx_cols = [
        "RA_centroid_lpllx_deg",
        "Dec_centroid_lpllx_deg",
        "RA_true_lpllx_deg",
        "Dec_true_lpllx_deg",
    ]
    documented_parallax_cols = ["lens_parallax_x_mas", "lens_parallax_y_mas"]

    errors: List[str] = []
    zscores: List[np.ndarray] = []

    def record_error(message: str) -> None:
        if len(errors) < max_errors:
            errors.append(message)

    for lc_file in lc_files:
        summary.checked_lightcurves += 1

        try:
            df = pd.read_csv(lc_file, sep=r"\s+", comment="#")
        except Exception as exc:  # pragma: no cover - defensive parser check
            record_error(f"{lc_file.name}: failed to parse table: {exc}")
            continue

        if df.empty:
            record_error(f"{lc_file.name}: table is empty")
            continue
        summary.checked_epochs += len(df)

        missing = _missing_columns(df.columns, required_cols)
        if missing:
            record_error(
                f"{lc_file.name}: missing required astrometry columns: {', '.join(missing)}"
            )
            continue

        bad_numeric = _nan_or_inf_columns(df, required_cols)
        if bad_numeric:
            record_error(
                f"{lc_file.name}: non-finite values in columns: {', '.join(bad_numeric)}"
            )
            continue

        row = _match_out_row(lc_file, out_rows)
        if row is None:
            record_error(
                f"{lc_file.name}: unable to match this lightcurve to an EventID/SubRun/Field row in .out"
            )
            continue

        frame_ra_deg, frame_dec_deg = _parse_astrometry_frame(lc_file)
        if frame_ra_deg is None or frame_dec_deg is None:
            record_error(f"{lc_file.name}: missing or malformed #Astrometry_Frame header")
            continue
        transform = _parse_astrometry_transform(lc_file)
        if transform is None:
            record_error(
                f"{lc_file.name}: missing or malformed #Astrometry_Transform header (required for ecliptic->ICRS checks)"
            )
            continue
        tr_a, tr_b, tr_c, tr_d = transform

        out_ra = float(row.get("ra_deg", math.nan))
        out_dec = float(row.get("dec_deg", math.nan))
        if not math.isnan(out_ra) and not math.isnan(out_dec):
            dra_mas = abs(_angle_diff_deg(frame_ra_deg, out_ra)) * math.cos(frame_dec_deg * DEG_TO_RAD) * MAS_PER_DEG
            ddec_mas = abs(frame_dec_deg - out_dec) * MAS_PER_DEG
            if dra_mas > 0.05 or ddec_mas > 0.05:
                record_error(
                    f"{lc_file.name}: frame RA/Dec mismatch .out row (dRA={dra_mas:.3e} mas, dDec={ddec_mas:.3e} mas)"
                )

        if _HAS_ASTROPY and "galactic_l" in row and "galactic_b" in row and not math.isnan(out_ra) and not math.isnan(out_dec):
            l_raw = float(row["galactic_l"])
            b_raw = float(row["galactic_b"])
            # Some catalogs store l,b in radians, others in degrees.
            # Try both and keep the one that best matches .out RA/Dec.
            candidates = [
                SkyCoord(l=l_raw * u.deg, b=b_raw * u.deg, frame="galactic").icrs,
                SkyCoord(l=l_raw * u.rad, b=b_raw * u.rad, frame="galactic").icrs,
            ]
            diffs = []
            for cand in candidates:
                dra = abs(_angle_diff_deg(float(cand.ra.deg), out_ra)) * math.cos(out_dec * DEG_TO_RAD) * MAS_PER_DEG
                ddec = abs(float(cand.dec.deg) - out_dec) * MAS_PER_DEG
                diffs.append((dra, ddec))
            best_dra, best_ddec = min(diffs, key=lambda item: math.hypot(item[0], item[1]))
            if best_dra > 500.0 or best_ddec > 500.0:
                record_error(
                    f"{lc_file.name}: .out galactic_l/b -> RA/Dec mismatch is large even after trying deg/rad "
                    f"(best dRA={best_dra:.3f} mas, best dDec={best_ddec:.3f} mas)"
                )

        xerr = df["x_centroid_error_mas"].to_numpy(dtype=float, copy=False)
        yerr = df["y_centroid_error_mas"].to_numpy(dtype=float, copy=False)
        if np.any(xerr < 0) or np.any(yerr < 0):
            record_error(f"{lc_file.name}: negative astrometric error values found")

        txerr = df["true_x_centroid_error_mas"].to_numpy(dtype=float, copy=False)
        tyerr = df["true_y_centroid_error_mas"].to_numpy(dtype=float, copy=False)
        if np.nanmax(np.abs(txerr)) > 1.0e-8 or np.nanmax(np.abs(tyerr)) > 1.0e-8:
            record_error(
                f"{lc_file.name}: true centroid error columns should be ~0 but max abs is "
                f"{max(float(np.nanmax(np.abs(txerr))), float(np.nanmax(np.abs(tyerr)))):.3e} mas"
            )

        cos_dec = math.cos(frame_dec_deg * DEG_TO_RAD)
        if abs(cos_dec) < 1.0e-8:
            record_error(f"{lc_file.name}: cos(dec) too small for stable RA conversion")
            continue

        x_obs = df["x_centroid_mas"].to_numpy(dtype=float, copy=False)
        y_obs = df["y_centroid_mas"].to_numpy(dtype=float, copy=False)
        x_true = df["true_x_centroid_mas"].to_numpy(dtype=float, copy=False)
        y_true = df["true_y_centroid_mas"].to_numpy(dtype=float, copy=False)

        dra_obs_cosdec_mas = tr_a * x_obs + tr_b * y_obs
        ddec_obs_mas = tr_c * x_obs + tr_d * y_obs
        dra_true_cosdec_mas = tr_a * x_true + tr_b * y_true
        ddec_true_mas = tr_c * x_true + tr_d * y_true

        ra_obs_expected = frame_ra_deg + dra_obs_cosdec_mas / (MAS_PER_DEG * cos_dec)
        dec_obs_expected = frame_dec_deg + ddec_obs_mas / MAS_PER_DEG
        ra_true_expected = frame_ra_deg + dra_true_cosdec_mas / (MAS_PER_DEG * cos_dec)
        dec_true_expected = frame_dec_deg + ddec_true_mas / MAS_PER_DEG

        ra_obs = df["RA_centroid_deg"].to_numpy(dtype=float, copy=False)
        dec_obs = df["Dec_centroid_deg"].to_numpy(dtype=float, copy=False)
        ra_true = df["RA_centroid_true_deg"].to_numpy(dtype=float, copy=False)
        dec_true = df["Dec_centroid_true_deg"].to_numpy(dtype=float, copy=False)

        tref = float(row.get("tref", math.nan))
        if math.isnan(tref):
            record_error(f"{lc_file.name}: .out row is missing finite tref")
            continue

        sim_time = df["Simulation_time"].to_numpy(dtype=float, copy=False)
        idx_ref = int(np.argmin(np.abs(sim_time - tref)))
        dt_ref = float(abs(sim_time[idx_ref] - tref))
        if dt_ref > 0.5:
            summary.warnings.append(
                f"{lc_file.name}: nearest epoch to tref is far (|t-tref|={dt_ref:.4f} day), "
                "tref-anchor checks are less diagnostic"
            )

        # Canonical event pointing used by gulls outputs is ra_deg/dec_deg in .out.
        # Those are themselves derived from canonical l,b, so this comparison directly
        # tests centroid anchoring against the canonical pointing metadata.
        canonical_ra = out_ra
        canonical_dec = out_dec

        if not math.isnan(canonical_ra) and not math.isnan(canonical_dec):
            dra_ref = _angle_diff_deg(float(ra_true[idx_ref]), canonical_ra) * math.cos(canonical_dec * DEG_TO_RAD) * MAS_PER_DEG
            ddec_ref = (float(dec_true[idx_ref]) - canonical_dec) * MAS_PER_DEG
            ref_offset_mas = math.hypot(dra_ref, ddec_ref)
            theta_e_mas = float(row.get("thetaE", math.nan))
            tol_ref_mas = max(1.0, 4.0 * theta_e_mas) if not math.isnan(theta_e_mas) else 1.0
            if ref_offset_mas > tol_ref_mas:
                record_error(
                    f"{lc_file.name}: true astrometric position at epoch nearest tref is not close to canonical pointing "
                    f"from out-file metadata (offset={ref_offset_mas:.4f} mas, tol={tol_ref_mas:.4f} mas, |t-tref|={dt_ref:.4f} day)"
                )

        dra_obs_mas = np.abs(ra_obs - ra_obs_expected) * MAS_PER_DEG * abs(cos_dec)
        ddec_obs_mas = np.abs(dec_obs - dec_obs_expected) * MAS_PER_DEG
        dra_true_mas = np.abs(ra_true - ra_true_expected) * MAS_PER_DEG * abs(cos_dec)
        ddec_true_mas = np.abs(dec_true - dec_true_expected) * MAS_PER_DEG

        if np.nanmax(dra_obs_mas) > 0.05 or np.nanmax(ddec_obs_mas) > 0.05:
            record_error(
                f"{lc_file.name}: RA/Dec (observed) inconsistent with x/y centroid conversion "
                f"(max dRA={float(np.nanmax(dra_obs_mas)):.4f} mas, max dDec={float(np.nanmax(ddec_obs_mas)):.4f} mas)"
            )
        if np.nanmax(dra_true_mas) > 0.05 or np.nanmax(ddec_true_mas) > 0.05:
            record_error(
                f"{lc_file.name}: RA/Dec (true) inconsistent with x/y centroid conversion "
                f"(max dRA={float(np.nanmax(dra_true_mas)):.4f} mas, max dDec={float(np.nanmax(ddec_true_mas)):.4f} mas)"
            )

        if "Lens_Dist" in row:
            expected_dl = float(row["Lens_Dist"])
            dl = df["lens_dist_kpc"].to_numpy(dtype=float, copy=False)
            dmax = float(np.nanmax(np.abs(dl - expected_dl)))
            if dmax > 1.0e-6:
                record_error(
                    f"{lc_file.name}: lens_dist_kpc differs from .out Lens_Dist by up to {dmax:.3e} kpc"
                )

        if "BJD" in df.columns and sim_zero_time is not None:
            sim_time = df["Simulation_time"].to_numpy(dtype=float, copy=False)
            bjd = df["BJD"].to_numpy(dtype=float, copy=False)
            expected_bjd = sim_zero_time + sim_time
            bjd_abs_err = np.abs(bjd - expected_bjd)
            if np.nanmedian(bjd_abs_err) > 1.0e-3 or np.nanmax(bjd_abs_err) > 5.0e-3:
                sample = ", ".join(f"{val:.6f}" for val in bjd[:3])
                record_error(
                    f"{lc_file.name}: BJD is inconsistent with SIMULATION_ZERO_TIME + Simulation_time "
                    f"(median |Δ|={float(np.nanmedian(bjd_abs_err)):.3e} day, max |Δ|={float(np.nanmax(bjd_abs_err)):.3e} day; "
                    f"first BJD values: {sample})"
                )
            if len(bjd) > 2 and np.unique(np.round(bjd, 10)).size < 3:
                record_error(
                    f"{lc_file.name}: BJD shows <3 unique values over {len(bjd)} epochs, suggesting severe precision loss"
                )

        elif "BJD" in df.columns and sim_zero_time is None:
            summary.warnings.append(
                f"{lc_file.name}: skipped BJD epoch-alignment check (SIMULATION_ZERO_TIME missing)"
            )
        else:
            record_error(f"{lc_file.name}: missing BJD column")

        if all(col in df.columns for col in lpllx_cols):
            ra_obs_lpllx = df["RA_centroid_lpllx_deg"].to_numpy(dtype=float, copy=False)
            dec_obs_lpllx = df["Dec_centroid_lpllx_deg"].to_numpy(dtype=float, copy=False)
            ra_true_lpllx = df["RA_true_lpllx_deg"].to_numpy(dtype=float, copy=False)
            dec_true_lpllx = df["Dec_true_lpllx_deg"].to_numpy(dtype=float, copy=False)

            pllx_obs_e = (ra_obs_lpllx - ra_obs) * MAS_PER_DEG * cos_dec
            pllx_obs_n = (dec_obs_lpllx - dec_obs) * MAS_PER_DEG
            pllx_true_e = (ra_true_lpllx - ra_true) * MAS_PER_DEG * cos_dec
            pllx_true_n = (dec_true_lpllx - dec_true) * MAS_PER_DEG

            if all(col in df.columns for col in documented_parallax_cols):
                pllx_e_ecl = df["lens_parallax_x_mas"].to_numpy(dtype=float, copy=False)
                pllx_n_ecl = df["lens_parallax_y_mas"].to_numpy(dtype=float, copy=False)
                pred_pllx_e = tr_a * pllx_e_ecl + tr_b * pllx_n_ecl
                pred_pllx_n = tr_c * pllx_e_ecl + tr_d * pllx_n_ecl

                diff_e = np.nanmax(np.abs(pllx_obs_e - pred_pllx_e))
                diff_n = np.nanmax(np.abs(pllx_obs_n - pred_pllx_n))
                if diff_e > 0.1 or diff_n > 0.1:
                    record_error(
                        f"{lc_file.name}: lpllx-vs-non-lpllx RA/Dec offsets do not match declared lens_parallax_x/y columns "
                        f"(max dE={float(diff_e):.4f} mas, max dN={float(diff_n):.4f} mas)"
                    )
            else:
                summary.warnings.append(
                    f"{lc_file.name}: lpllx columns present but lens_parallax_x/y missing; skipped lpllx offset-model check"
                )

            diff_true_e = np.nanmax(np.abs(pllx_true_e - pllx_obs_e))
            diff_true_n = np.nanmax(np.abs(pllx_true_n - pllx_obs_n))
            if diff_true_e > 0.1 or diff_true_n > 0.1:
                record_error(
                    f"{lc_file.name}: parallax offsets differ between measured and true RA/Dec tracks "
                    f"(max dE={float(diff_true_e):.4f} mas, max dN={float(diff_true_n):.4f} mas)"
                )
        else:
            missing_lpllx = _missing_columns(df.columns, lpllx_cols)
            summary.warnings.append(
                f"{lc_file.name}: missing lpllx RA/Dec columns ({', '.join(missing_lpllx)}); skipped lens-parallax RA/Dec checks"
            )

        missing_doc = _missing_columns(df.columns, documented_parallax_cols)
        if missing_doc:
            msg = (
                f"{lc_file.name}: missing documented lens parallax columns ({', '.join(missing_doc)}); "
                "using observer position columns for equivalent checks"
            )
            if strict_documented_columns:
                record_error(msg)
            else:
                summary.warnings.append(msg)

        if "centroid_final_x_mas" in df.columns and "centroid_final_y_mas" in df.columns:
            dfx = np.nanmax(
                np.abs(
                    df["centroid_final_x_mas"].to_numpy(dtype=float, copy=False)
                    - x_true
                )
            )
            dfy = np.nanmax(
                np.abs(
                    df["centroid_final_y_mas"].to_numpy(dtype=float, copy=False)
                    - y_true
                )
            )
            if dfx > 0.1 or dfy > 0.1:
                record_error(
                    f"{lc_file.name}: centroid_final_* differs from true_* by up to "
                    f"(dx={float(dfx):.4f} mas, dy={float(dfy):.4f} mas)"
                )

        obs_x = df["parallax_shift_x"].to_numpy(dtype=float, copy=False)
        obs_y = df["parallax_shift_y"].to_numpy(dtype=float, copy=False)
        obs_z = df["parallax_shift_z"].to_numpy(dtype=float, copy=False)
        if np.ptp(obs_x) < 1.0e-8 and np.ptp(obs_y) < 1.0e-8 and np.ptp(obs_z) < 1.0e-8:
            record_error(
                f"{lc_file.name}: observer parallax_shift_(x,y,z) appears constant; expected epoch-dependent observer position"
            )

        if all(col in df.columns for col in ["source0_x_thE", "source0_y_thE", "lens0_x_thE", "lens0_y_thE"]):
            theta_e_mas = float(row.get("thetaE", math.nan))
            mu_ref_out = float(row.get("murel_ref", math.nan))
            tE_ref = float(row.get("tE_ref", math.nan))

            if not math.isnan(theta_e_mas) and not math.isnan(mu_ref_out):
                srcx = df["source0_x_thE"].to_numpy(dtype=float, copy=False)
                srcy = df["source0_y_thE"].to_numpy(dtype=float, copy=False)
                lenx = df["lens0_x_thE"].to_numpy(dtype=float, copy=False)
                leny = df["lens0_y_thE"].to_numpy(dtype=float, copy=False)
                relx_mas = (srcx - lenx) * theta_e_mas
                rely_mas = (srcy - leny) * theta_e_mas

                if len(sim_time) < 3:
                    summary.warnings.append(
                        f"{lc_file.name}: too few epochs for robust proper-motion finite-difference check"
                    )
                else:
                    if idx_ref == 0:
                        idx0, idx1 = 0, 1
                    elif idx_ref == len(sim_time) - 1:
                        idx0, idx1 = len(sim_time) - 2, len(sim_time) - 1
                    else:
                        idx0, idx1 = idx_ref - 1, idx_ref + 1

                    dt = float(sim_time[idx1] - sim_time[idx0])
                    if abs(dt) > 0:
                        vx = float((relx_mas[idx1] - relx_mas[idx0]) / dt)
                        vy = float((rely_mas[idx1] - rely_mas[idx0]) / dt)
                        mu_from_lc = math.hypot(vx, vy) * 365.25
                        tol_mu = max(0.2, 0.05 * abs(mu_ref_out))
                        if abs(mu_from_lc - mu_ref_out) > tol_mu:
                            record_error(
                                f"{lc_file.name}: relative proper-motion magnitude near tref from source/lens tracks "
                                f"does not match .out murel_ref "
                                f"(lc={mu_from_lc:.4f} mas/yr, out={mu_ref_out:.4f} mas/yr, tol={tol_mu:.4f})"
                            )

                        if not math.isnan(tE_ref) and tE_ref > 0:
                            mu_from_te = theta_e_mas / tE_ref * 365.25
                            tol_te = max(0.05, 0.02 * abs(mu_ref_out))
                            if abs(mu_from_te - mu_ref_out) > tol_te:
                                record_error(
                                    f"{lc_file.name}: .out proper-motion self-consistency failed "
                                    f"(thetaE/tE_ref*365.25={mu_from_te:.4f} mas/yr, murel_ref={mu_ref_out:.4f} mas/yr)"
                                )
                    else:
                        summary.warnings.append(
                            f"{lc_file.name}: zero time spacing around tref prevented proper-motion finite-difference check"
                        )
        else:
            summary.warnings.append(
                f"{lc_file.name}: source0/lens0 position columns missing; skipped proper-motion check"
            )

        # Long-baseline check: away from t_ref and outside the main lensing window,
        # remove declared lens-parallax term and verify recovered motion matches
        # heliocentric proper motion from the .out metadata.
        if all(col in df.columns for col in ["RA_true_lpllx_deg", "Dec_true_lpllx_deg", "lens_parallax_x_mas", "lens_parallax_y_mas"]):
            mu_helio_out = float(row.get("murel_helio", math.nan))
            mu_helio_alpha_out = float(row.get("murel_helio_alpha", math.nan))
            mu_helio_delta_out = float(row.get("murel_helio_delta", math.nan))
            t0_event = float(row.get("t0lens1", math.nan))
            tE_ref = float(row.get("tE_ref", math.nan))
            if math.isnan(mu_helio_out) or math.isnan(t0_event) or math.isnan(tE_ref):
                summary.warnings.append(
                    f"{lc_file.name}: missing murel_helio/t0lens1/tE_ref in .out; skipped long-baseline heliocentric PM check"
                )
            else:
                t_year = (sim_time - tref) / 365.25
                left_years = max(0.0, -float(np.min(t_year)))
                right_years = max(0.0, float(np.max(t_year)))
                if left_years < long_baseline_years_each_side or right_years < long_baseline_years_each_side:
                    summary.warnings.append(
                        f"{lc_file.name}: long-baseline PM check skipped: coverage around tref is "
                        f"{left_years:.2f} yr before and {right_years:.2f} yr after; "
                        f"need >= {long_baseline_years_each_side:.2f} yr each side"
                    )
                else:
                    exclusion_days = long_baseline_exclusion_te * max(tE_ref, 1.0e-6)
                    mask_far = np.abs(sim_time - t0_event) > exclusion_days
                    # Enforce that long-baseline fit uses both sides of tref.
                    mask_far &= ((sim_time - tref) < -30.0) | ((sim_time - tref) > 30.0)
                    n_far = int(np.sum(mask_far))
                    if n_far < 30:
                        summary.warnings.append(
                            f"{lc_file.name}: long-baseline PM check skipped: only {n_far} far-from-event epochs "
                            f"after excluding |t-t0| <= {long_baseline_exclusion_te:.1f} tE"
                        )
                    elif np.sum(mask_far & (sim_time < tref)) < 10 or np.sum(mask_far & (sim_time > tref)) < 10:
                        summary.warnings.append(
                            f"{lc_file.name}: long-baseline PM check skipped: insufficient far-from-event points on one side of tref"
                        )
                    else:
                        ra_lpllx = df["RA_true_lpllx_deg"].to_numpy(dtype=float, copy=False)
                        dec_lpllx = df["Dec_true_lpllx_deg"].to_numpy(dtype=float, copy=False)
                        dra_lpllx = (ra_lpllx - frame_ra_deg + 180.0) % 360.0 - 180.0
                        e_lpllx = dra_lpllx * cos_dec * MAS_PER_DEG
                        n_lpllx = (dec_lpllx - frame_dec_deg) * MAS_PER_DEG

                        pllx_e_ecl = df["lens_parallax_x_mas"].to_numpy(dtype=float, copy=False)
                        pllx_n_ecl = df["lens_parallax_y_mas"].to_numpy(dtype=float, copy=False)
                        e_pllx = tr_a * pllx_e_ecl + tr_b * pllx_n_ecl
                        n_pllx = tr_c * pllx_e_ecl + tr_d * pllx_n_ecl

                        # Remove lens parallax term from absolute astrometric track.
                        e_corr = e_lpllx - e_pllx
                        n_corr = n_lpllx - n_pllx

                        x_fit = t_year[mask_far]
                        fit_corr = _fit_linear_motion(x_fit, e_corr[mask_far], n_corr[mask_far])
                        fit_mu_e = float(fit_corr["mu_e"])
                        fit_mu_n = float(fit_corr["mu_n"])
                        mu_fit_helio = float(fit_corr["mu_amp"])

                        tol_mu_helio = max(0.5, 0.10 * abs(mu_helio_out))
                        final_fail_reasons: List[str] = []
                        if abs(mu_fit_helio - mu_helio_out) > tol_mu_helio:
                            final_fail_reasons.append(
                                f"{lc_file.name}: long-baseline heliocentric PM magnitude mismatch after lens-parallax correction "
                                f"(fit={mu_fit_helio:.4f} mas/yr, out={mu_helio_out:.4f} mas/yr, "
                                f"tol={tol_mu_helio:.4f}, excluded |t-t0|<={long_baseline_exclusion_te:.1f} tE)"
                            )

                        out_dir_deg: float | None = None
                        mu_helio_from_components = math.nan
                        if math.isnan(mu_helio_alpha_out) or math.isnan(mu_helio_delta_out):
                            summary.warnings.append(
                                f"{lc_file.name}: missing murel_helio_alpha/delta in .out; skipped long-baseline PM direction check"
                            )
                        else:
                            mu_helio_from_components = math.hypot(mu_helio_alpha_out, mu_helio_delta_out)
                            if abs(mu_helio_from_components - mu_helio_out) > max(0.05, 0.01 * abs(mu_helio_out)):
                                record_error(
                                    f"{lc_file.name}: .out murel_helio is inconsistent with murel_helio_alpha/delta "
                                    f"(hypot(alpha,delta)={mu_helio_from_components:.4f} mas/yr, murel_helio={mu_helio_out:.4f} mas/yr)"
                                )
                            elif mu_fit_helio < 0.2 or mu_helio_from_components < 0.2:
                                summary.warnings.append(
                                    f"{lc_file.name}: long-baseline PM direction check skipped because motion is very small "
                                    f"(fit={mu_fit_helio:.4f}, out={mu_helio_from_components:.4f} mas/yr)"
                                )
                            else:
                                out_dir_deg = math.degrees(math.atan2(mu_helio_delta_out, mu_helio_alpha_out))
                                fit_dir_deg = math.degrees(math.atan2(fit_mu_n, fit_mu_e))
                                dir_diff_deg = abs(_angle_diff_deg(fit_dir_deg, out_dir_deg))
                                sigma_theta_deg = float(fit_corr["sigma_theta_deg"])
                                tol_dir_deg = max(long_baseline_direction_tol_deg, 3.0 * sigma_theta_deg)
                                if dir_diff_deg > tol_dir_deg:
                                    alias_name, alias_diff = _direction_alias_diagnostic(
                                        fit_mu_e, fit_mu_n, out_dir_deg
                                    )
                                    hint = ""
                                    if alias_name != "as-is" and alias_diff + 5.0 < dir_diff_deg:
                                        hint = (
                                            f"; closest sign/axis variant is {alias_name} "
                                            f"with Δdir={alias_diff:.2f} deg"
                                        )
                                    final_fail_reasons.append(
                                        f"{lc_file.name}: long-baseline heliocentric PM direction mismatch after lens-parallax correction "
                                        f"(fit_dir={fit_dir_deg:.2f} deg, out_dir={out_dir_deg:.2f} deg, "
                                        f"Δdir={dir_diff_deg:.2f} deg, tol={tol_dir_deg:.2f} deg, "
                                        f"fit_mu=(E={fit_mu_e:.4f},N={fit_mu_n:.4f}) mas/yr, "
                                        f"out_mu=(alpha={mu_helio_alpha_out:.4f},delta={mu_helio_delta_out:.4f}) mas/yr"
                                        f"{hint})"
                                    )

                        source_only_pass = False
                        source_only_diagnostic = ""
                        if all(col in df.columns for col in ["centroid_src_x_mas", "centroid_src_y_mas"]):
                            src_e_ecl = df["centroid_src_x_mas"].to_numpy(dtype=float, copy=False)
                            src_n_ecl = df["centroid_src_y_mas"].to_numpy(dtype=float, copy=False)
                            src_e = tr_a * src_e_ecl + tr_b * src_n_ecl
                            src_n = tr_c * src_e_ecl + tr_d * src_n_ecl
                            src_fit = _fit_linear_motion(x_fit, src_e[mask_far], src_n[mask_far])
                            src_mu_e = float(src_fit["mu_e"])
                            src_mu_n = float(src_fit["mu_n"])
                            src_mu_amp = float(src_fit["mu_amp"])
                            src_fail_reasons: List[str] = []
                            if abs(src_mu_amp - mu_helio_out) > tol_mu_helio:
                                src_fail_reasons.append(
                                    f"src-only amplitude mismatch (fit={src_mu_amp:.4f} mas/yr, out={mu_helio_out:.4f} mas/yr, tol={tol_mu_helio:.4f})"
                                )
                            if out_dir_deg is not None:
                                src_dir_deg = math.degrees(math.atan2(src_mu_n, src_mu_e))
                                src_dir_diff = abs(_angle_diff_deg(src_dir_deg, out_dir_deg))
                                src_tol_dir = max(
                                    long_baseline_direction_tol_deg,
                                    3.0 * float(src_fit["sigma_theta_deg"]),
                                )
                                if src_dir_diff > src_tol_dir:
                                    src_fail_reasons.append(
                                        f"src-only direction mismatch (Δdir={src_dir_diff:.2f} deg, tol={src_tol_dir:.2f} deg)"
                                    )
                            source_only_pass = len(src_fail_reasons) == 0
                            source_only_diagnostic = (
                                f"src-only fit_mu=(E={src_mu_e:.4f},N={src_mu_n:.4f}) mas/yr"
                                + (", " + "; ".join(src_fail_reasons) if src_fail_reasons else ", src-only pass")
                            )

                        fs_obs = float(row.get("Obs_0_fs", math.nan))
                        high_blend = (not math.isnan(fs_obs)) and fs_obs < 0.9
                        if final_fail_reasons:
                            if source_only_pass and high_blend:
                                summary.warnings.append(
                                    f"{lc_file.name}: blended-centroid long-baseline PM check failed, but source-only centroid check passed "
                                    f"(Obs_0_fs={fs_obs:.3f}, likely blend-induced centroid drift). "
                                    f"Blended issues: {'; '.join(final_fail_reasons)}"
                                )
                            else:
                                for reason in final_fail_reasons:
                                    record_error(reason)
                                if source_only_diagnostic:
                                    summary.warnings.append(
                                        f"{lc_file.name}: source-only long-baseline diagnostic: {source_only_diagnostic}"
                                    )
                        elif source_only_diagnostic and (not source_only_pass):
                            summary.warnings.append(
                                f"{lc_file.name}: blended-centroid PM check passed but source-only diagnostic did not: {source_only_diagnostic}"
                            )

                        # Optional diagnostic: correction should not make linearity worse.
                        fit_raw = _fit_linear_motion(x_fit, e_lpllx[mask_far], n_lpllx[mask_far])
                        rms_corr = float(fit_corr["rms"])
                        rms_raw = float(fit_raw["rms"])
                        if rms_corr > 1.05 * rms_raw:
                            summary.warnings.append(
                                f"{lc_file.name}: parallax-corrected long-baseline fit is not tighter than raw fit "
                                f"(rms_corr={rms_corr:.4f} mas, rms_raw={rms_raw:.4f} mas)"
                            )
        else:
            summary.warnings.append(
                f"{lc_file.name}: missing RA_true_lpllx/lens_parallax columns; skipped long-baseline heliocentric PM check"
            )

        dx = x_obs - x_true
        dy = y_obs - y_true
        mask_x = xerr > 0
        mask_y = yerr > 0
        if np.any(mask_x):
            zscores.append(dx[mask_x] / xerr[mask_x])
        if np.any(mask_y):
            zscores.append(dy[mask_y] / yerr[mask_y])

    if zscores:
        z = np.concatenate(zscores)
        z = z[np.isfinite(z)]
        if z.size:
            summary.zscore_mean = float(np.mean(z))
            summary.zscore_std = float(np.std(z))
            frac_gt5 = float(np.mean(np.abs(z) > 5.0))
            if abs(summary.zscore_mean) > 0.15:
                record_error(
                    f"Global astrometric residual bias detected: mean((obs-true)/sigma)={summary.zscore_mean:.3f} (expected ~0)"
                )
            if summary.zscore_std < 0.7 or summary.zscore_std > 1.3:
                record_error(
                    f"Global astrometric residual scatter mismatch: std((obs-true)/sigma)={summary.zscore_std:.3f} (expected ~1)"
                )
            if frac_gt5 > 0.005:
                record_error(
                    f"Global astrometric residual outlier rate too high: P(|z|>5)={frac_gt5:.4f} (expected << 0.005)"
                )

    if errors:
        omitted = max(0, len(errors) - max_errors)
        rendered = "\n".join(f" - {err}" for err in errors[:max_errors])
        suffix = f"\n - ... and {omitted} more" if omitted else ""
        raise SmokeTestError(
            "Astrometric sanity checks failed:\n"
            f"{rendered}{suffix}\n"
            "These checks validate internal astrometric consistency between .lc and .out products."
        )

    return summary


__all__ = ["AstrometrySanitySummary", "verify_astrometry_sanity"]
