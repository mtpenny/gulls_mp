#!/usr/bin/env python3
"""Run BAGLE joint photometric+astrometric sanity fits on prepared outputs."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Sequence

if __package__ is None or __package__ == "":
    repo_root = Path(__file__).resolve().parents[1]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

from smoke_test.bagle_fit_sanity import run_bagle_joint_fit_sanity  # pylint: disable=wrong-import-position
from smoke_test.constants import REPO_ROOT  # pylint: disable=wrong-import-position
from smoke_test.errors import SmokeTestError  # pylint: disable=wrong-import-position
from smoke_test.prep import parse_parameter_file  # pylint: disable=wrong-import-position


def _discover_run_dirs(root: Path) -> List[Path]:
    runs: List[Path] = []
    for path in sorted(root.rglob("*")):
        if not path.is_dir():
            continue
        if list(path.glob("*.out")) and list(path.rglob("*.lc")):
            runs.append(path)
    return runs


def _resolve_run_dirs(raw_dirs: Sequence[Path]) -> List[Path]:
    if raw_dirs:
        return [path.resolve() for path in raw_dirs]
    return _discover_run_dirs((REPO_ROOT / "smoke_test" / "output").resolve())


def _infer_param_file(run_dir: Path) -> Path | None:
    param_dir = (REPO_ROOT / "smoke_test" / "parameterfiles").resolve()
    names = [run_dir.name]
    if "_" in run_dir.name:
        parts = run_dir.name.split("_")
        if len(parts) >= 2:
            names.append("_".join(parts[:2]))
    for name in names:
        candidate = param_dir / f"{name}.prm"
        if candidate.is_file():
            return candidate
    return None


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "run_dirs",
        nargs="*",
        type=Path,
        help=(
            "Run directory/directories containing .out and .lc files "
            "(default: auto-discover under smoke_test/output)"
        ),
    )
    parser.add_argument(
        "--params",
        type=Path,
        help=(
            "Optional parameter file used to provide SIMULATION_ZERO_TIME "
            "(if omitted, inferred from run directory name)."
        ),
    )
    parser.add_argument(
        "--event-id",
        type=int,
        default=None,
        help=(
            "Optional explicit EventID to fit. If omitted, the script chooses one "
            "single-source event with ObsGroup_0_chi2 < --chi2-max."
        ),
    )
    parser.add_argument(
        "--chi2-max",
        type=float,
        default=20.0,
        help="Upper bound on raw single-lens chi2 delta in .out (default: %(default)s).",
    )
    parser.add_argument(
        "--n-live-points",
        type=int,
        default=200,
        help="Nested-sampling live points for BAGLE MicrolensSolver (default: %(default)s).",
    )
    parser.add_argument(
        "--obs-location",
        type=str,
        default="earth",
        help="Observer location alias passed to BAGLE obsLocation (default: %(default)s).",
    )
    parser.add_argument(
        "--max-phot-points",
        type=int,
        default=0,
        help="Maximum photometric epochs used in BAGLE fit; <=0 uses all epochs (default: %(default)s).",
    )
    parser.add_argument(
        "--max-ast-points",
        type=int,
        default=0,
        help="Maximum astrometric epochs used in BAGLE fit; <=0 uses all epochs (default: %(default)s).",
    )
    parser.add_argument(
        "--fit-true-astrometry",
        action="store_true",
        help=(
            "Fit BAGLE to noiseless astrometry columns instead of noisy astrometry columns "
            "(still validates against noiseless astrometry)."
        ),
    )
    parser.add_argument(
        "--true-ast-err-mas",
        type=float,
        default=0.01,
        help="When --fit-true-astrometry is set, use this fixed astrometric uncertainty in mas.",
    )
    parser.add_argument(
        "--fit-reduced-chi2-max",
        type=float,
        default=3.0,
        help="Fail if BAGLE reduced chi2 exceeds this value (default: %(default)s).",
    )
    parser.add_argument(
        "--mu-amp-frac-tol",
        type=float,
        default=0.20,
        help="Diagnostic fractional proper-motion amplitude threshold recorded in the BAGLE summary (default: %(default)s).",
    )
    parser.add_argument(
        "--mu-dir-tol-deg",
        type=float,
        default=10.0,
        help="Diagnostic proper-motion direction threshold in degrees recorded in the BAGLE summary (default: %(default)s).",
    )
    parser.add_argument(
        "--pie-amp-frac-tol",
        type=float,
        default=0.35,
        help="Diagnostic fractional parallax amplitude threshold recorded in the BAGLE summary (default: %(default)s).",
    )
    parser.add_argument(
        "--pie-dir-tol-deg",
        type=float,
        default=15.0,
        help="Diagnostic parallax direction threshold in degrees recorded in the BAGLE summary (default: %(default)s).",
    )
    parser.add_argument(
        "--true-ast-rms-mas-max",
        type=float,
        default=0.05,
        help="Fail if BAGLE-vs-noiseless astrometric RMS exceeds this (mas).",
    )
    parser.add_argument(
        "--true-ast-sigma-max",
        type=float,
        default=1.5,
        help=(
            "Fail if BAGLE-vs-noiseless astrometric RMS exceeds this many median astrometric "
            "error sigmas."
        ),
    )
    parser.add_argument(
        "--lens-ast-rms-demean-mas-max",
        type=float,
        default=0.05,
        help=(
            "Fail if BAGLE get_lens_astrometry disagrees with GULLS primary-lens astrometry by more than this "
            "RMS (mas) after removing constant x/y offsets."
        ),
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    run_dirs = _resolve_run_dirs(args.run_dirs)
    if not run_dirs:
        print("No run directories found with both .out and .lc files.")
        return 1

    failures = 0
    for run_dir in run_dirs:
        out_files = sorted(run_dir.glob("*.out"))
        if not out_files:
            print(f"[SKIP] {run_dir}: no .out files")
            continue

        if args.params is not None:
            param_path = args.params.resolve()
        else:
            inferred = _infer_param_file(run_dir)
            if inferred is None:
                failures += 1
                print(f"\n[FAIL] {run_dir}")
                print(
                    "Could not infer parameter file for this run. "
                    "Pass --params to provide SIMULATION_ZERO_TIME."
                )
                continue
            param_path = inferred

        params: Dict[str, str] = parse_parameter_file(param_path)
        print(f"\n=== BAGLE joint-fit sanity: {run_dir} ===")
        print(f"Using params: {param_path}")

        try:
            summary = run_bagle_joint_fit_sanity(
                run_dir,
                out_files,
                params,
                chi2_max=args.chi2_max,
                event_id=args.event_id,
                max_phot_points=args.max_phot_points,
                max_ast_points=args.max_ast_points,
                n_live_points=args.n_live_points,
                fit_true_astrometry=args.fit_true_astrometry,
                true_ast_err_mas=args.true_ast_err_mas,
                obs_location=args.obs_location,
                fit_reduced_chi2_max=args.fit_reduced_chi2_max,
                mu_amp_frac_tol=args.mu_amp_frac_tol,
                mu_dir_tol_deg=args.mu_dir_tol_deg,
                piE_amp_frac_tol=args.pie_amp_frac_tol,
                piE_dir_tol_deg=args.pie_dir_tol_deg,
                true_ast_rms_mas_max=args.true_ast_rms_mas_max,
                true_ast_sigma_max=args.true_ast_sigma_max,
                lens_ast_rms_demean_mas_max=args.lens_ast_rms_demean_mas_max,
            )
        except SmokeTestError as err:
            failures += 1
            print(f"[FAIL] {run_dir}")
            print(err)
            continue

        print(
            f"[PASS] event={summary.event_id} (SubRun={summary.subrun}, Field={summary.field}) "
            f"from {summary.lc_file.name}"
        )
        print(
            "       astrometry fit mode={}{}".format(
                "noiseless" if args.fit_true_astrometry else "noisy",
                (
                    f" (fixed err={args.true_ast_err_mas:.4g} mas)"
                    if args.fit_true_astrometry
                    else ""
                ),
            )
        )
        print(
            f"       selected single-lens chi2={summary.out_chi2_single_lens:.4f}, "
            f"fit reduced chi2={summary.fit_reduced_chi2:.4f}"
        )
        print(
            "       mu_rel out=(E={:.4f},N={:.4f}) fit=(E={:.4f},N={:.4f}) "
            "amp(out/fit)=({:.4f}/{:.4f}) dirΔ={:.2f} deg".format(
                summary.mu_ref_e,
                summary.mu_ref_n,
                summary.mu_fit_e,
                summary.mu_fit_n,
                summary.mu_amp_ref,
                summary.mu_amp_fit,
                summary.mu_dir_diff_deg,
            )
        )
        print(
            "       piE out=(E={:.4f},N={:.4f}) fit=(E={:.4f},N={:.4f}) "
            "amp(out/fit)=({:.4f}/{:.4f}) dirΔ={:.2f} deg".format(
                summary.piE_ref_e,
                summary.piE_ref_n,
                summary.piE_fit_e,
                summary.piE_fit_n,
                summary.piE_amp_ref,
                summary.piE_amp_fit,
                summary.piE_dir_diff_deg,
            )
        )
        print(
            "       true-astrometry: RMS={:.4f} mas, sigma_equiv={:.2f}".format(
                summary.true_ast_rms_mas,
                summary.true_ast_sigma_equiv,
            )
        )
        print(f"       BAGLE obsLocation={summary.obs_location_used!r}")
        if summary.lens_ast_rms_demean_mas is not None:
            print(
                "       lens-track: RMS_raw={:.4f} mas, RMS_after_xy_offset={:.4f} mas".format(
                    summary.lens_ast_rms_raw_mas,
                    summary.lens_ast_rms_demean_mas,
                )
            )
        if summary.lens_plot_path is not None:
            print(f"       lens plot: {summary.lens_plot_path}")
        print(f"       plot: {summary.plot_path}")
        print(f"       summary: {summary.result_json_path}")
        for warning in summary.warnings:
            print(f"  [WARN] {warning}")

    if failures:
        print(f"\nBAGLE joint-fit sanity failed for {failures} run(s).")
        return 1

    print("\nBAGLE joint-fit sanity passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
