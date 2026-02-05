#!/usr/bin/env python3
"""Run astrometric sanity checks against pre-generated gulls outputs."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Sequence

if __package__ is None or __package__ == "":
    repo_root = Path(__file__).resolve().parents[1]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

from smoke_test.astrometry_sanity import verify_astrometry_sanity  # pylint: disable=wrong-import-position
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
        help="Optional parameter file used to provide SIMULATION_ZERO_TIME for BJD checks.",
    )
    parser.add_argument(
        "--max-errors",
        type=int,
        default=40,
        help="Maximum number of per-file failures shown before truncation (default: %(default)s).",
    )
    parser.add_argument(
        "--strict-documented-columns",
        action="store_true",
        help="Fail if documented lens_parallax_x_mas/y_mas columns are missing.",
    )
    parser.add_argument(
        "--long-baseline-years",
        type=float,
        default=1.0,
        help=(
            "Required coverage (years) on each side of tref for long-baseline "
            "heliocentric PM checks (default: %(default)s)."
        ),
    )
    parser.add_argument(
        "--long-baseline-exclusion-te",
        type=float,
        default=5.0,
        help=(
            "Exclude epochs within this many tE of t0 when fitting long-baseline "
            "heliocentric PM (default: %(default)s)."
        ),
    )
    parser.add_argument(
        "--long-baseline-direction-tol-deg",
        type=float,
        default=15.0,
        help=(
            "Minimum angular tolerance (deg) for long-baseline heliocentric PM "
            "direction checks (default: %(default)s)."
        ),
    )
    return parser.parse_args(argv)


def _resolve_run_dirs(raw_dirs: Sequence[Path]) -> List[Path]:
    if raw_dirs:
        return [path.resolve() for path in raw_dirs]
    return _discover_run_dirs((REPO_ROOT / "smoke_test" / "output").resolve())


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    run_dirs = _resolve_run_dirs(args.run_dirs)
    if not run_dirs:
        print("No run directories found with both .out and .lc files.")
        return 1

    params: Dict[str, str] | None = None
    if args.params is not None:
        params = parse_parameter_file(args.params.resolve())

    failures = 0
    for run_dir in run_dirs:
        out_files = sorted(run_dir.glob("*.out"))
        if not out_files:
            print(f"[SKIP] {run_dir}: no .out files")
            continue
        if not list(run_dir.rglob("*.lc")):
            print(f"[SKIP] {run_dir}: no .lc files")
            continue

        print(f"\n=== Astrometry sanity: {run_dir} ===")
        try:
            summary = verify_astrometry_sanity(
                run_dir,
                out_files,
                params,
                max_errors=args.max_errors,
                strict_documented_columns=args.strict_documented_columns,
                long_baseline_years_each_side=args.long_baseline_years,
                long_baseline_exclusion_te=args.long_baseline_exclusion_te,
                long_baseline_direction_tol_deg=args.long_baseline_direction_tol_deg,
            )
        except SmokeTestError as err:
            failures += 1
            print(f"[FAIL] {run_dir}")
            print(err)
            continue

        print(
            f"[PASS] checked {summary.checked_lightcurves} lightcurves, "
            f"{summary.checked_epochs} epochs"
        )
        if summary.zscore_mean is not None and summary.zscore_std is not None:
            print(
                f"       residual z-score mean={summary.zscore_mean:.4f}, "
                f"std={summary.zscore_std:.4f}"
            )
        max_warn = 20
        for warning in summary.warnings[:max_warn]:
            print(f"  [WARN] {warning}")
        if len(summary.warnings) > max_warn:
            print(f"  [WARN] ... and {len(summary.warnings) - max_warn} more warning(s)")

    if failures:
        print(f"\nAstrometry sanity failed for {failures} run(s).")
        return 1

    print("\nAstrometry sanity checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
