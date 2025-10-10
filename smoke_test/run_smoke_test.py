#!/usr/bin/env python3
"""Minimal end-to-end smoke test runner for gulls executables.

The script executes one or more gulls binaries using the synthetic inputs that
live under `smoke_test/`. Each selected executable must finish without errors
and emit both a `.out` report and at least one `.lc` light-curve file in its
configured output directory.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

try:
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    PLOTTING_AVAILABLE = True
except ImportError:
    PLOTTING_AVAILABLE = False

REPO_ROOT = Path(__file__).resolve().parents[1]
BUILD_BIN_DEFAULT = REPO_ROOT / "bin"
PARAM_DIR = REPO_ROOT / "smoke_test" / "parameterfiles"

CASES: Tuple[Tuple[str, str], ...] = (
    ("gulls_std.x", "smoke_std.prm"),
    ("gulls_croin.x", "smoke_croin.prm"),
    ("gullsFish.x", "smoke_fish.prm"),
)


class SmokeTestError(RuntimeError):
    """Raised when the smoke test encounters a setup or runtime failure."""


@dataclass
class PreparedCase:
    name: str
    exe_path: Path
    param_path: Path
    run_name: str
    output_root: Path
    output_dir: Path


def parse_args() -> argparse.Namespace:
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
        choices=[name for name, _ in CASES],
        help="Subset of executables to run (default: all).",
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
    return parser.parse_args()


def ensure_executable(path: Path) -> None:
    if not path.is_file():
        raise SmokeTestError(f"Missing executable: {path}")
    if not os.access(path, os.X_OK):
        raise SmokeTestError(f"Executable is not runnable: {path}")


def run_command(cmd: Sequence[str], env: Dict[str, str], timeout: float | None) -> subprocess.CompletedProcess[str]:
    effective_timeout = None if timeout is None or timeout <= 0 else timeout
    return subprocess.run(
        cmd,
        cwd=REPO_ROOT,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
        timeout=effective_timeout,
    )


def parse_parameter_file(path: Path) -> Dict[str, str]:
    params: Dict[str, str] = {}
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                continue
            key, value = line.split("=", 1)
            params[key.strip()] = value.strip()
    return params


def prepare_cases(build_bin: Path, selected: Sequence[Tuple[str, str]]) -> Tuple[List[PreparedCase], List[str]]:
    prepared: List[PreparedCase] = []
    failures: List[str] = []

    for name, prm_filename in selected:
        exe_path = build_bin / name
        try:
            ensure_executable(exe_path)
        except SmokeTestError as exc:
            failures.append(str(exc))
            continue

        param_path = PARAM_DIR / prm_filename
        if not param_path.is_file():
            failures.append(f"Parameter file not found: {param_path}")
            continue

        params = parse_parameter_file(param_path)
        run_name = params.get("RUN_NAME", param_path.stem)
        output_dir_value = params.get("OUTPUT_DIR")
        if not output_dir_value:
            failures.append(f"{name}: OUTPUT_DIR missing in {param_path}")
            continue

        output_root = (REPO_ROOT / output_dir_value).resolve()
        output_dir = output_root / run_name
        prepared.append(PreparedCase(name, exe_path, param_path, run_name, output_root, output_dir))

    return prepared, failures


def verify_outputs(output_dir: Path) -> None:
    out_files = sorted(output_dir.glob("*.out"))
    lc_files = list(output_dir.rglob("*.lc"))
    if not out_files:
        raise SmokeTestError(f"No .out files found in {output_dir}")
    if not lc_files:
        raise SmokeTestError(f"No .lc files found under {output_dir}")
    if out_files[0].stat().st_size == 0:
        raise SmokeTestError(f"Summary file is empty: {out_files[0]}")


def plot_lightcurves(output_dir: Path) -> None:
    """Generate photometry and astrometry plots from lightcurve files."""
    if not PLOTTING_AVAILABLE:
        return
    
    lc_files = list(output_dir.rglob("*.lc"))
    if not lc_files:
        return
    
    for lc_file in lc_files:
        try:
            # Read the file with pandas, using whitespace delimiter and skipping comment lines
            df = pd.read_csv(lc_file, sep=r'\s+', comment='#')
            
            if df.empty:
                continue
            
            # Extract time and photometry columns by name
            time = df['Simulation_time'].values
            flux = df['measured_relative_flux'].values
            flux_err = df['measured_relative_flux_error'].values
            true_flux = df['true_relative_flux'].values if 'true_relative_flux' in df.columns else None
            
            # Check if astrometry columns exist
            astrom_cols = ['true_N_centroid_mas', 'true_E_centroid_mas', 
                          'measured_N_centroid_mas', 'measured_E_centroid_mas',
                          'measured_N_centroid_error_mas', 'measured_E_centroid_error_mas',
                          'true_centroid_ra_deg', 'true_centroid_dec_deg',
                          'measured_centroid_ra_deg', 'measured_centroid_dec_deg',
                          'measured_centroid_ra_error_deg', 'measured_centroid_dec_error_deg']
            has_astrom = all(col in df.columns for col in astrom_cols)
            
            if has_astrom:
                # Extract astrometry columns by name
                true_N_mas = df['true_N_centroid_mas'].values
                true_E_mas = df['true_E_centroid_mas'].values
                meas_N_mas = df['measured_N_centroid_mas'].values
                meas_E_mas = df['measured_E_centroid_mas'].values
                meas_N_err_mas = df['measured_N_centroid_error_mas'].values
                meas_E_err_mas = df['measured_E_centroid_error_mas'].values
                true_ra_deg = df['true_centroid_ra_deg'].values
                true_dec_deg = df['true_centroid_dec_deg'].values
                meas_ra_deg = df['measured_centroid_ra_deg'].values
                meas_dec_deg = df['measured_centroid_dec_deg'].values
                meas_ra_err_deg = df['measured_centroid_ra_error_deg'].values
                meas_dec_err_deg = df['measured_centroid_dec_error_deg'].values
            
            # Create figure with subplots
            if has_astrom:
                fig, axes = plt.subplots(2, 2, figsize=(12, 10))
                fig.suptitle(f'Smoke Test: {lc_file.stem}', fontsize=14)
                
                # Photometry plot (top-left)
                ax = axes[0, 0]
                # Plot measured with errorbars (bottom layer)
                ax.errorbar(time, flux, yerr=flux_err, fmt='o', markersize=2, 
                           alpha=0.5, color='C0', label='Measured', zorder=1)
                # Plot true flux as connected line (middle layer)
                if true_flux is not None:
                    ax.plot(time, true_flux, '-', linewidth=1.5, color='red', 
                           label='True', zorder=2, alpha=0.8)
                # Baseline on top
                ax.axhline(1.0, color='k', linestyle='--', linewidth=1.5, 
                          label='Baseline', zorder=3)
                ax.set_xlabel('Time (days)')
                ax.set_ylabel('Relative Flux')
                ax.set_title('Light Curve')
                ax.legend()
                ax.grid(True, alpha=0.3)
                
                # RA/Dec position plot (top-right)
                ax = axes[0, 1]
                # Convert to milliarcsec for better visibility
                meas_ra_mas = (meas_ra_deg - true_ra_deg[0]) * 3600 * 1000
                meas_dec_mas = (meas_dec_deg - true_dec_deg[0]) * 3600 * 1000
                true_ra_mas = (true_ra_deg - true_ra_deg[0]) * 3600 * 1000
                true_dec_mas = (true_dec_deg - true_dec_deg[0]) * 3600 * 1000
                meas_ra_err_mas = meas_ra_err_deg * 3600 * 1000
                meas_dec_err_mas = meas_dec_err_deg * 3600 * 1000
                
                # Plot with error bars and connecting line
                ax.errorbar(meas_ra_mas, meas_dec_mas, 
                           xerr=meas_ra_err_mas, yerr=meas_dec_err_mas,
                           fmt='o-', markersize=2, alpha=0.5, color='red',
                           linewidth=1, label='Measured', capsize=2)
                ax.plot(true_ra_mas, true_dec_mas, 'b-', linewidth=2, 
                       alpha=0.7, label='True')
                ax.set_xlabel('ΔRA (mas)')
                ax.set_ylabel('ΔDec (mas)')
                ax.set_title('Astrometric Position (RA/Dec)')
                ax.legend()
                ax.grid(True, alpha=0.3)
                ax.axis('equal')
                
                # Astrometry: North-East trajectory (bottom-left)
                ax = axes[1, 0]
                ax.plot(true_E_mas, true_N_mas, 'b-', label='True', linewidth=2, alpha=0.7)
                ax.plot(meas_E_mas, meas_N_mas, 'r.', label='Measured', markersize=3, alpha=0.5)
                ax.set_xlabel('East (mas)')
                ax.set_ylabel('North (mas)')
                ax.set_title('Astrometric Trajectory (N/E)')
                ax.legend()
                ax.grid(True, alpha=0.3)
                ax.axis('equal')
                
                # Astrometry: Time series (bottom-right)
                ax = axes[1, 1]
                ax.plot(time, true_N_mas, 'b-', label='True N', linewidth=2, alpha=0.7)
                ax.plot(time, meas_N_mas, 'r.', label='Meas N', markersize=2, alpha=0.5)
                ax.plot(time, true_E_mas, 'g-', label='True E', linewidth=2, alpha=0.7)
                ax.plot(time, meas_E_mas, 'm.', label='Meas E', markersize=2, alpha=0.5)
                ax.set_xlabel('Time (days)')
                ax.set_ylabel('Centroid Shift (mas)')
                ax.set_title('Astrometry vs Time')
                ax.legend(fontsize=8)
                ax.grid(True, alpha=0.3)
            else:
                # Photometry only
                fig, ax = plt.subplots(1, 1, figsize=(10, 6))
                fig.suptitle(f'Smoke Test: {lc_file.stem}', fontsize=14)
                ax.errorbar(time, flux, yerr=flux_err, fmt='o', markersize=2, 
                           alpha=0.5, color='C0', label='Measured', zorder=1)
                if true_flux is not None:
                    ax.plot(time, true_flux, '-', linewidth=1.5, color='red', 
                           label='True', zorder=2, alpha=0.8)
                ax.axhline(1.0, color='k', linestyle='--', linewidth=1.5, 
                          label='Baseline', zorder=3)
                ax.set_xlabel('Time (days)')
                ax.set_ylabel('Relative Flux')
                ax.set_title('Light Curve')
                ax.legend()
                ax.grid(True, alpha=0.3)
            
            # Save plot
            plot_file = output_dir / f"{lc_file.stem}_plot.png"
            plt.tight_layout()
            plt.savefig(plot_file, dpi=150, bbox_inches='tight')
            plt.close()
            print(f"  Generated plot: {plot_file.name}")
            
        except Exception as e:
            print(f"  Warning: Could not plot {lc_file.name}: {e}")
            continue


def main() -> int:
    args = parse_args()
    build_bin = args.build_bin.resolve()
    if not build_bin.is_dir():
        raise SmokeTestError(f"Build directory not found: {build_bin}")

    selected = CASES if not args.cases else tuple((name, prm) for name, prm in CASES if name in args.cases)
    if not selected:
        print("No cases selected", file=sys.stderr)
        return 1

    espl_table = REPO_ROOT / "src" / "ESPL.tbl"
    if not espl_table.is_file():
        raise SmokeTestError(f"Missing ESPL.tbl at {espl_table}; copy it before running the smoke test.")

    prepared_cases, prep_failures = prepare_cases(build_bin, selected)
    if prep_failures:
        print("Smoke test setup issues:")
        for message in prep_failures:
            print(f" - {message}")
    if not prepared_cases:
        return 1

    env = os.environ.copy()
    base_dir = REPO_ROOT.as_posix() + "/"
    env["GULLS_BASE_DIR"] = base_dir
    env.setdefault("GULLS_STARS_DIR", base_dir)

    output_roots = {case.output_root for case in prepared_cases}
    if not args.keep_output:
        for root in output_roots:
            if root.exists():
                shutil.rmtree(root)
            root.mkdir(parents=True, exist_ok=True)
    else:
        for root in output_roots:
            root.mkdir(parents=True, exist_ok=True)

    failures: List[str] = []

    for case in prepared_cases:
        if case.output_dir.exists() and not args.keep_output:
            shutil.rmtree(case.output_dir)
        case.output_dir.mkdir(parents=True, exist_ok=True)

        cmd = [str(case.exe_path), "-i", str(case.param_path), "-s", args.instance]
        if args.field is not None:
            cmd.extend(["-f", str(args.field)])

        print(f"\n=== Running {case.name} with {case.param_path.name} ===")
        try:
            result = run_command(cmd, env, args.exec_timeout)
        except subprocess.TimeoutExpired as exc:
            if exc.stdout:
                print(exc.stdout)
            print(f"{case.name} exceeded {args.exec_timeout} seconds and was terminated.")
            failures.append(f"{case.name} timed out after {args.exec_timeout}s")
            continue

        print(result.stdout)
        if result.returncode != 0:
            failures.append(f"{case.name} exited with {result.returncode}")
            continue

        try:
            verify_outputs(case.output_dir)
            plot_lightcurves(case.output_dir)
        except SmokeTestError as exc:
            failures.append(f"{case.name}: {exc}")

    if failures:
        print("\nSmoke test failed:")
        for message in failures:
            print(f" - {message}")
        return 1

    print("\nAll smoke test cases completed successfully.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except SmokeTestError as exc:
        print(f"Smoke test setup error: {exc}", file=sys.stderr)
        sys.exit(1)
