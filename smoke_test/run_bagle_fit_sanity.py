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
        default=100.0,
        help="Selection threshold for single-lens chi2 in .out (default: %(default)s).",
    )
    parser.add_argument(
        "--n-live-points",
        type=int,
        default=200,
        help="Nested-sampling live points for BAGLE MicrolensSolver (default: %(default)s).",
    )
    parser.add_argument(
        "--max-phot-points",
        type=int,
        default=500,
        help="Maximum photometric epochs used in BAGLE fit (default: %(default)s).",
    )
    parser.add_argument(
        "--max-ast-points",
        type=int,
        default=500,
        help="Maximum astrometric epochs used in BAGLE fit (default: %(default)s).",
    )
    parser.add_argument(
        "--fit-reduced-chi2-max",
        type=float,
        default=6.0,
        help="Fail if BAGLE reduced chi2 exceeds this value (default: %(default)s).",
    )
    parser.add_argument(
        "--mu-amp-frac-tol",
        type=float,
        default=0.35,
        help="Allowed fractional proper-motion amplitude mismatch (default: %(default)s).",
    )
    parser.add_argument(
        "--mu-dir-tol-deg",
        type=float,
        default=20.0,
        help="Allowed proper-motion direction mismatch in degrees (default: %(default)s).",
    )
    parser.add_argument(
        "--pie-amp-frac-tol",
        type=float,
        default=0.70,
        help="Allowed fractional parallax amplitude mismatch (default: %(default)s).",
    )
    parser.add_argument(
        "--pie-dir-tol-deg",
        type=float,
        default=30.0,
        help="Allowed parallax direction mismatch in degrees (default: %(default)s).",
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
                fit_reduced_chi2_max=args.fit_reduced_chi2_max,
                mu_amp_frac_tol=args.mu_amp_frac_tol,
                mu_dir_tol_deg=args.mu_dir_tol_deg,
                piE_amp_frac_tol=args.pie_amp_frac_tol,
                piE_dir_tol_deg=args.pie_dir_tol_deg,
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

