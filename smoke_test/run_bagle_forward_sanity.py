"""Command-line entry point for the explicit 1s1l BAGLE forward-model sanity check."""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

if __package__ in (None, ""):
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from smoke_test.bagle_forward_sanity import run_bagle_forward_1s1l_sanity
    from smoke_test.errors import SmokeTestError
else:
    from .bagle_forward_sanity import run_bagle_forward_1s1l_sanity
    from .errors import SmokeTestError


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path, help="Smoke-test run directory to inspect.")
    parser.add_argument(
        "--event-id",
        type=int,
        default=None,
        help="Explicit EventID to compare (default: first 1s1l event in the .out file).",
    )
    parser.add_argument(
        "--obs-location",
        type=str,
        default="earth",
        help="Observer location passed to BAGLE obsLocation (default: %(default)s).",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        summary = run_bagle_forward_1s1l_sanity(
            args.run_dir,
            event_id=args.event_id,
            obs_location=args.obs_location,
        )
    except SmokeTestError as exc:
        print(exc)
        return 1

    print(f"[PASS] event={summary.event_id} from {summary.lc_file.name}")
    print(f"       BAGLE model = {summary.bagle_model_name}")
    print(f"       BAGLE forward astrometry RMS = {summary.astrometry_rms_mas:.6f} mas")
    print(
        "       BAGLE forward photometry RMS = "
        f"{summary.photometry_rms_relative_flux:.6f} relflux "
        f"(max abs diff {summary.photometry_max_abs_relative_flux_diff:.6f})"
    )
    print(f"       summary: {summary.summary_json_path}")
    print(f"       plot: {summary.plot_path}")
    for warning in summary.warnings:
        print(f"  [WARN] {warning}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
