"""Command-line entry point for the gulls smoke test."""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

from .constants import (
    BUILD_BIN_DEFAULT,
    CASE_CHOICES,
    CASE_EXEC_MAP,
    CASE_LABELS,
    CASE_LOOKUP,
    CASES,
    REPO_ROOT,
)
from .astrometry_sanity import verify_astrometry_sanity
from .bagle_forward_sanity import run_bagle_forward_1s1l_sanity
from .bagle_fit_sanity import run_bagle_joint_fit_sanity
from .errors import SmokeTestError
from .execution import run_command
from .metrics import gather_case_metrics
from .plotting import plot_lightcurves
from .prep import PreparedCase, prepare_cases
from .validation import (
    verify_binary_source_columns,
    verify_catalog_alignment,
    verify_catalog_columns,
    verify_input_files_exist,
    verify_planet_file_schema,
    verify_psf_files,
    verify_nfilters_matches_catalogs,
    verify_outputs,
    verify_rates_file,
    verify_sequence_has_observations,
    verify_source_lens_compatibility,
    verify_weather_coverage,
)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--build-bin",
        type=Path,
        default=BUILD_BIN_DEFAULT,
        help="Directory containing the gulls executables (default: %(default)s)",
    )
    parser.add_argument(
        "--keep-output",
        action="store_true",
        help="Skip cleaning existing smoke_test/output directories before running.",
    )
    parser.add_argument(
        "--cases",
        nargs="*",
        choices=list(CASE_CHOICES),
        help="Subset of cases to run (accepts case labels or executable names; default: all).",
    )
    parser.add_argument(
        "--instance",
        default="0",
        help="Instance identifier passed via -s (default: %(default)s).",
    )
    parser.add_argument(
        "--field",
        type=int,
        default=0,
        help="Field index passed via -f (default: %(default)s). Use --field -1 to let gulls auto-select.",
    )
    parser.add_argument(
        "--exec-timeout",
        type=float,
        default=180.0,
        help="Seconds to wait for each executable before aborting (<=0 disables).",
    )
    parser.add_argument(
        "--ci",
        action="store_true",
        help="Run CI-optimized subset of tests (faster, essential cases only).",
    )
    parser.add_argument(
        "--debug",
        type=int,
        default=0,
        help="The debug verbosity level of the output"
    )
    parser.add_argument(
        "--astrometry-sanity",
        action="store_true",
        help=(
            "Run strict astrometry self-consistency checks on generated .lc/.out products. "
            "Useful for validating pre-generated astrometric outputs."
        ),
    )
    parser.add_argument(
        "--astrometry-strict-documented-columns",
        action="store_true",
        help=(
            "With --astrometry-sanity, fail if documented lens_parallax_x_mas/y_mas "
            "columns are missing."
        ),
    )
    parser.add_argument(
        "--astrometry-long-baseline-years",
        type=float,
        default=1.0,
        help=(
            "Required years of coverage on each side of tref for long-baseline "
            "heliocentric proper-motion checks when --astrometry-sanity is enabled "
            "(default: %(default)s)."
        ),
    )
    parser.add_argument(
        "--astrometry-long-baseline-exclusion-te",
        type=float,
        default=5.0,
        help=(
            "Exclude epochs within this many tE of t0 for long-baseline heliocentric "
            "proper-motion checks when --astrometry-sanity is enabled (default: %(default)s)."
        ),
    )
    parser.add_argument(
        "--astrometry-long-baseline-direction-tol-deg",
        type=float,
        default=15.0,
        help=(
            "Minimum angular tolerance (deg) for long-baseline heliocentric "
            "proper-motion direction checks when --astrometry-sanity is enabled "
            "(default: %(default)s)."
        ),
    )
    parser.add_argument(
        "--bagle-joint-fit-sanity",
        action="store_true",
        help=(
            "Run an additional BAGLE joint photometric+astrometric fit sanity check "
            "on one single-source, single-lens-like event in each case output."
        ),
    )
    parser.add_argument(
        "--bagle-forward-sanity",
        action="store_true",
        help=(
            "Run an additional BAGLE forward-model sanity check "
            "(no fitting) on one explicit 1s1l event in each case output."
        ),
    )
    parser.add_argument(
        "--bagle-event-id",
        type=int,
        default=None,
        help=(
            "With --bagle-joint-fit-sanity, force BAGLE fit to use this EventID "
            "(default: auto-select using single-lens chi2)."
        ),
    )
    parser.add_argument(
        "--bagle-chi2-max",
        type=float,
        default=100.0,
        help=(
            "With --bagle-joint-fit-sanity, maximum .out ObsGroup_0_chi2 for event "
            "selection (default: %(default)s)."
        ),
    )
    parser.add_argument(
        "--bagle-forward-event-id",
        type=int,
        default=None,
        help=(
            "With --bagle-forward-sanity, force BAGLE forward comparison to use "
            "this EventID (default: auto-select first explicit 1s1l event)."
        ),
    )
    parser.add_argument(
        "--bagle-forward-obs-location",
        default="earth",
        help=(
            "With --bagle-forward-sanity, BAGLE observer location (default: %(default)s)."
        ),
    )
    parser.add_argument(
        "--bagle-n-live-points",
        type=int,
        default=200,
        help=(
            "With --bagle-joint-fit-sanity, nested-sampling live points "
            "(default: %(default)s)."
        ),
    )
    return parser.parse_args(argv)


def _resolve_case_selection(raw_choices: Sequence[str] | None, ci_mode: bool = False) -> Tuple[Tuple[str, str, str], ...]:
    if not raw_choices:
        if ci_mode:
            # CI subset: essential tests only (std, binary source validation)
            return (
                ("smoke_std", "gulls_std.x", "smoke_std.prm"),
                ("smoke_std_binary", "gulls_std.x", "smoke_std_binary.prm"),
                ("smoke_fish", "gullsFish.x", "smoke_fish.prm"),
                # For the fish binary CI case we should use the fish binary
                # parameter file so outputs land under the fish/ output tree.
                ("smoke_fish_binary", "gullsFish.x", "smoke_fish_binary.prm"),
                ("smoke_croin", "gulls_croin.x", "smoke_croin.prm"),
                ("smoke_croin_binary", "gulls_croin.x", "smoke_croin_binary.prm"),
                ("smoke_general", "gulls_general.x", "smoke_general.prm"),
                ("smoke_general_binary", "gulls_general.x", "smoke_general_binary.prm"),
            )
        return CASES

    ordered: List[Tuple[str, str, str]] = []
    for choice in raw_choices:
        if choice in CASE_LABELS:
            ordered.append(CASE_LOOKUP[choice])
        else:
            ordered.extend(CASE_EXEC_MAP.get(choice, []))

    seen: set[str] = set()
    deduped: List[Tuple[str, str, str]] = []
    for case in ordered:
        if case[0] in seen:
            continue
        deduped.append(case)
        seen.add(case[0])
    return tuple(deduped)


def _prepare_environment(keep_output: bool, prepared_cases: Sequence[PreparedCase]) -> None:
    output_roots = {case.output_root for case in prepared_cases}
    if not keep_output:
        for root in output_roots:
            if root.exists():
                shutil.rmtree(root)
            root.mkdir(parents=True, exist_ok=True)
    else:
        for root in output_roots:
            root.mkdir(parents=True, exist_ok=True)


def _generate_psf_files(build_bin: Path) -> None:
    """Generate PSF files needed for smoke tests."""
    psf_dir = REPO_ROOT / "smoke_test" / "assets" / "observatories"
    psf_binary_file = psf_dir / "WFI_PSF.psf"
    
    # Check if we already have a valid PSF file (should be ~68MB for subpixel sampling)
    if psf_binary_file.exists() and psf_binary_file.stat().st_size > 10_000_000:  # > 10MB
        print(f"Using existing PSF file: {psf_binary_file} ({psf_binary_file.stat().st_size:,} bytes)")
        return
    
    # Generate PSF using PSF class with Moffat function
    generate_moffat_psf = build_bin / "generateMoffatPSF"
    
    if not generate_moffat_psf.exists():
        raise SmokeTestError(f"PSF generator not found: {generate_moffat_psf}")
    
    # Generate binary PSF with subpixel sampling
    print("Generating PSF with subpixel sampling using PSF class...")
    cmd = [
        str(generate_moffat_psf),
        "0.2",    # fwhm (arcsec) - matches PSFFWHM
        "0.11",   # pixel_scale (arcsec) - matches PIXELSCALE
        str(psf_binary_file)
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=REPO_ROOT)
    if result.returncode != 0:
        raise SmokeTestError(f"PSF generation failed: {result.stderr}")
    
    if not psf_binary_file.exists():
        raise SmokeTestError(f"PSF file was not created: {psf_binary_file}")
    
    print(f"Generated PSF file: {psf_binary_file} ({psf_binary_file.stat().st_size:,} bytes)")


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    build_bin = args.build_bin.resolve()
    if not build_bin.is_dir():
        raise SmokeTestError(f"Build directory not found: {build_bin}")

    selected = _resolve_case_selection(args.cases, args.ci)
    if not selected:
        print("No cases selected", file=sys.stderr)
        return 1

    espl_table = REPO_ROOT / "src" / "ESPL.tbl"
    if not espl_table.is_file():
        raise SmokeTestError(f"Missing ESPL.tbl at {espl_table}; copy it before running the smoke test.")

    # Generate PSF files needed for smoke tests
    _generate_psf_files(build_bin)

    prepared_cases, prep_failures = prepare_cases(build_bin, selected)
    if prep_failures:
        print("Smoke test setup issues:")
        for message in prep_failures:
            print(f" - {message}")

    if not prepared_cases:
        return 1
    
    # Validate catalogs: columns, compatibility, and binary-specific requirements
    for case in prepared_cases:
        try:
            verify_input_files_exist(case.params)
            verify_psf_files(case.params)
            verify_catalog_columns(case.params)
            verify_source_lens_compatibility(case.params)
            verify_binary_source_columns(case.params)
            verify_nfilters_matches_catalogs(case.params)
            verify_weather_coverage(case.params)
            verify_rates_file(case.params)
            verify_sequence_has_observations(case.params)
            verify_planet_file_schema(case.params, args.field if args.field is not None else -1, args.instance)
        except SmokeTestError as err:
            print(f"Validation failed for {case.label}:")
            print(f" - {err}")
            return 1

    env = os.environ.copy()
    base_dir = REPO_ROOT.as_posix() + "/"
    env["GULLS_BASE_DIR"] = base_dir
    env.setdefault("GULLS_STARS_DIR", base_dir)

    _prepare_environment(args.keep_output, prepared_cases)

    failures: List[str] = []

    for case in prepared_cases:
        if case.output_dir.exists() and not args.keep_output:
            shutil.rmtree(case.output_dir)
        case.output_dir.mkdir(parents=True, exist_ok=True)

        cmd = [str(case.exe_path), "-i", str(case.exec_param_path), "-s", args.instance] #, "-d","-d","-d"]
        if args.field is not None:
            cmd.extend(["-f", str(args.field)])
        for i in range(args.debug):
            cmd.extend(["-d"])

        print(f"\n=== Running {case.exe_name} ({case.label}) with {case.param_path.name} ===")
        result = run_command(cmd, env, args.exec_timeout, cwd=REPO_ROOT)

        print(result.stdout)
        if result.returncode != 0:
            failures.append(f"{case.label} exited with {result.returncode}")
            continue

        out_files = verify_outputs(case.output_dir)
        verify_catalog_alignment(out_files, case.params)
        if args.astrometry_sanity:
            summary = verify_astrometry_sanity(
                case.output_dir,
                out_files,
                case.params,
                strict_documented_columns=args.astrometry_strict_documented_columns,
                long_baseline_years_each_side=args.astrometry_long_baseline_years,
                long_baseline_exclusion_te=args.astrometry_long_baseline_exclusion_te,
                long_baseline_direction_tol_deg=args.astrometry_long_baseline_direction_tol_deg,
            )
            print(
                "Astrometry sanity passed: "
                f"{summary.checked_lightcurves} lightcurves, {summary.checked_epochs} epochs"
            )
            if summary.zscore_mean is not None and summary.zscore_std is not None:
                print(
                    "Astrometry residual z-score stats: "
                    f"mean={summary.zscore_mean:.4f}, std={summary.zscore_std:.4f}"
                )
            for report in summary.pm_conversion_reports:
                print(report)
            max_warn = 10
            for warning in summary.warnings[:max_warn]:
                print(f"Astrometry warning: {warning}")
            if len(summary.warnings) > max_warn:
                print(
                    "Astrometry warning: "
                    f"... and {len(summary.warnings) - max_warn} more warnings"
                )
        if args.bagle_joint_fit_sanity:
            bagle_summary = run_bagle_joint_fit_sanity(
                case.output_dir,
                out_files,
                case.params,
                chi2_max=args.bagle_chi2_max,
                event_id=args.bagle_event_id,
                n_live_points=args.bagle_n_live_points,
            )
            print(
                "BAGLE joint-fit sanity passed: "
                f"event={bagle_summary.event_id} (SubRun={bagle_summary.subrun}, Field={bagle_summary.field}), "
                f"fit reduced chi2={bagle_summary.fit_reduced_chi2:.4f}"
            )
            print(
                "BAGLE PM comparison: "
                f"out(E,N)=({bagle_summary.mu_ref_e:.4f},{bagle_summary.mu_ref_n:.4f}) mas/yr, "
                f"fit(E,N)=({bagle_summary.mu_fit_e:.4f},{bagle_summary.mu_fit_n:.4f}) mas/yr, "
                f"dirΔ={bagle_summary.mu_dir_diff_deg:.2f} deg"
            )
            print(
                "BAGLE parallax comparison: "
                f"out(E,N)=({bagle_summary.piE_ref_e:.4f},{bagle_summary.piE_ref_n:.4f}), "
                f"fit(E,N)=({bagle_summary.piE_fit_e:.4f},{bagle_summary.piE_fit_n:.4f}), "
                f"dirΔ={bagle_summary.piE_dir_diff_deg:.2f} deg"
            )
            for warning in bagle_summary.warnings:
                print(f"BAGLE warning: {warning}")
            print(f"BAGLE plot: {bagle_summary.plot_path}")
            print(f"BAGLE summary: {bagle_summary.result_json_path}")
        if args.bagle_forward_sanity:
            bagle_forward_summary = run_bagle_forward_1s1l_sanity(
                case.output_dir,
                event_id=args.bagle_forward_event_id,
                obs_location=args.bagle_forward_obs_location,
            )
            print(
                "BAGLE forward sanity passed: "
                f"event={bagle_forward_summary.event_id}, "
                f"model={bagle_forward_summary.bagle_model_name}"
            )
            print(
                "BAGLE forward astrometry RMS: "
                f"{bagle_forward_summary.astrometry_rms_mas:.6f} mas"
            )
            print(
                "BAGLE forward photometry RMS/max|Δflux|: "
                f"{bagle_forward_summary.photometry_rms_relative_flux:.6e} / "
                f"{bagle_forward_summary.photometry_max_abs_relative_flux_diff:.6e}"
            )
            for warning in bagle_forward_summary.warnings:
                print(f"BAGLE forward warning: {warning}")
            print(f"BAGLE forward plot: {bagle_forward_summary.plot_path}")
            print(f"BAGLE forward summary: {bagle_forward_summary.summary_json_path}")
        summaries = gather_case_metrics(out_files)
        plot_lightcurves(case.output_dir, summaries, case.params, build_bin)

    if failures:
        print("\nSmoke test failed:")
        for message in failures:
            print(f" - {message}")
        return 1

    print("\nAll smoke test cases completed successfully.")
    return 0


__all__ = ["main", "parse_args"]
