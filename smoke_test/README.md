# Gulls Smoke Test Suite

This directory contains a minimal end-to-end smoke test for the gulls executables: `gulls_std`, `gulls_croin`, `gullsFish`, and `gulls_general`. The smoke test verifies that the executables can run successfully with synthetic input data and produce expected outputs including lightcurves, reports, and visualizations.

## Overview

The smoke test uses synthetic catalogs and simplified observatory configurations to quickly validate that:
- All executables compile and run without errors
- Lightcurve generation completes within timeout limits
- Photometry and astrometry outputs are generated correctly
- Visualization plots can be created from the output data

## Directory Structure

```
smoke_test/
├── README.md                  # This file
├── run_smoke_test.py          # Main test runner script
├── parameterfiles/            # Parameter files for each executable
│   ├── smoke_std.prm
│   ├── smoke_croin.prm
│   ├── smoke_fish.prm
│   └── smoke_general.prm
├── assets/                    # Test input data
│   ├── lenses/                # Synthetic lens catalogs
│   ├── sources/               # Synthetic source catalogs
│   ├── starfields/            # Starfield definitions
│   ├── observatories/         # Observatory configurations
│   ├── rates/                 # Event rate files
│   ├── weather/               # Weather profiles
│   └── planets/               # Planet configurations
└── output/                    # Test outputs (generated at runtime)
    ├── std/
    ├── croin/
    ├── fish/
    └── general/
```

## Catalog Helper Tool

Use `dat_tool.py` to inspect and adjust the whitespace-delimited catalogs under `assets/`.

- Inspect column layout with optional sample values:
  ```bash
  python dat_tool.py show assets/lenses/smoke_lens_catalog.dat
  ```
- Scale or offset a numeric column (writes to a new file unless `--inplace` is provided):
  ```bash
  python dat_tool.py apply assets/lenses/smoke_lens_catalog.dat \
      --column Mass --operation multiply --value 1.5 \
      --output assets/lenses/smoke_lens_catalog_scaled.dat
  ```
- Use `--format` to control numeric formatting (default `"{value:.7e}"`) and `--preview`
  to choose how many modified rows are echoed.

## Requirements

- Built gulls executables in `bin/` (see main README for build instructions)
- Python 3.9+ with the following packages:
  - `numpy`
  - `pandas`
  - `matplotlib` (optional, for plotting)
  - `astropy`
  - `scipy`
  - `VBMicrolensing`

Additional packages for `--bagle-joint-fit-sanity`:
- `bagle`
- `joblib`
- `dynesty`
- `ultranest`
- `pymultinest` (optional at runtime; if MultiNest is unavailable, the test falls back to scipy least-squares)

### Setting Up Python Dependencies

#### Option 1: Use the main Gulls environment (Recommended)
If you've already set up the main Gulls environment, you can use it directly:
```bash
conda activate gulls  # Already includes all smoke test dependencies
```

#### Option 2: Dedicated smoke test environment
For isolated testing or if you have conflicting dependencies:
```bash
conda env create -f smoke_test/smoke.yml
conda activate smoke
```

#### Option 3: Using pip
Alternatively, install packages with pip:
```bash
pip install numpy pandas matplotlib astropy scipy VBMicrolensing

# Additional dependencies for BAGLE joint-fit sanity checks
pip install bagle joblib dynesty ultranest pymultinest
```

**Note**: The main `environment.yml` includes all smoke test dependencies, so you don't need a separate smoke environment unless you have specific dependency conflicts.

## Running the Smoke Test

### Basic Usage

Run all test cases:
```bash
python smoke_test/run_smoke_test.py
```

Run a specific executable:
```bash
python smoke_test/run_smoke_test.py --cases gulls_std
```

Run with custom timeout (default: 180 seconds):
```bash
python smoke_test/run_smoke_test.py --exec-timeout 120
```

### Command-Line Options

- `--build-bin PATH`: Directory containing executables (default: `build/bin`)
- `--keep-output`: Skip cleaning existing output directories before running
- `--cases CASE [CASE ...]`: Subset of runs to execute. Accepts either executable names (`gulls_std`, `gulls_croin`, `gullsFish`, `gulls_general`) or case labels (`std-single`, `std-binary`, `std-heavy`, `croin-single`, `croin-binary`, `croin-heavy`, `fish-single`, `fish-binary`, `fish-heavy`, `general-single`, `general-binary`, `general-1s1l`). Case-specific run names are appended automatically (for example, `smoke_std_std-heavy`) so the heavy scenarios do not overwrite the baseline outputs.
- `--instance ID`: Instance identifier passed via `-s` flag (default: `0`)
- `--field N`: Field index passed via `-f` flag (default: `0`; use `-1` for auto-select)
- `--exec-timeout SECONDS`: Timeout per executable (default: `180`; `<=0` disables)
- `--astrometry-sanity`: Run astrometric self-consistency checks on generated `.lc/.out` files
- `--astrometry-strict-documented-columns`: With `--astrometry-sanity`, fail if documented `lens_parallax_x_mas/y_mas` columns are missing
- `--astrometry-long-baseline-years`: Coverage required on each side of `tref` (years) for long-baseline heliocentric-PM checks (default: `1.0`)
- `--astrometry-long-baseline-exclusion-te`: Exclude `|t-t0| <= N*tE` for long-baseline heliocentric-PM fits (default: `5.0`)
- `--astrometry-long-baseline-direction-tol-deg`: Minimum angular tolerance for long-baseline heliocentric-PM direction checks (default: `15.0`)
- `--bagle-joint-fit-sanity`: Run an additional BAGLE combined photometry+astrometry fit sanity check
- `--bagle-forward-sanity`: Run an additional BAGLE forward-model sanity check (no fitting) on explicit 1s1l events
- `--bagle-event-id`: With `--bagle-joint-fit-sanity`, force a specific `EventID` (otherwise auto-select one with single-lens chi2 < threshold)
- `--bagle-chi2-max`: With `--bagle-joint-fit-sanity`, event-selection threshold on `ObsGroup_0_chi2` (default: `100.0`)
- `--bagle-forward-event-id`: With `--bagle-forward-sanity`, force a specific explicit 1s1l `EventID` (otherwise auto-select first explicit 1s1l event)
- `--bagle-forward-obs-location`: With `--bagle-forward-sanity`, observer location passed to BAGLE (default: `earth`)
- `--bagle-n-live-points`: With `--bagle-joint-fit-sanity`, BAGLE nested-sampling live points (default: `200`)

### Examples

```bash
# Test only gulls_std with 2-minute timeout
python smoke_test/run_smoke_test.py --cases gulls_std --exec-timeout 120

# Test all executables and keep previous outputs
python smoke_test/run_smoke_test.py --keep-output

# Test with custom build directory
python smoke_test/run_smoke_test.py --build-bin /path/to/custom/build/bin

# Run only general case plus astrometry sanity checks
python smoke_test/run_smoke_test.py --cases general-single --astrometry-sanity

# Run general case plus BAGLE joint-fit sanity check
python smoke_test/run_smoke_test.py \
  --cases general-single \
  --bagle-joint-fit-sanity \
  --bagle-chi2-max 100 \
  --bagle-n-live-points 250

# Run explicit 1s1l case plus BAGLE forward-model sanity check (no fitting)
python smoke_test/run_smoke_test.py \
  --cases general-1s1l \
  --bagle-forward-sanity \
  --bagle-forward-obs-location earth
```

## Standalone Astrometry Sanity Checks

You can validate pre-generated outputs without rerunning executables:

```bash
# Auto-discover run directories under smoke_test/output/
python smoke_test/run_astrometry_sanity.py

# Check one run directory with explicit parameter file (for BJD consistency checks)
python smoke_test/run_astrometry_sanity.py \
  smoke_test/output/general/smoke_general \
  --params smoke_test/parameterfiles/smoke_general.prm

# Require >=2 years on each side of tref for long-baseline checks
python smoke_test/run_astrometry_sanity.py \
  /path/to/run_dir \
  --long-baseline-years 2.0 \
  --long-baseline-exclusion-te 8.0 \
  --long-baseline-direction-tol-deg 10.0
```

The astrometry sanity checker validates:
- required astrometry column presence and finite values
- lightcurve `#Astrometry_Frame` consistency with canonical `.out` event metadata
- declared frame metadata (`#Astrometry_EventToEcl`, `#Astrometry_Transform`, `#Astrometry_BAGLE`) is present for explicit convention tracking
- astrometric position at epoch nearest `tref` is close to canonical event pointing from `.out`
- relative proper-motion magnitude near `tref` from lightcurve source/lens tracks is consistent with `.out` (`murel_ref`, `thetaE`, `tE_ref`)
- centroid `(E,N)` ↔ `RA/Dec` conversion consistency
- lens-parallax `RA/Dec` offsets against documented `lens_parallax_x_mas/y_mas` columns
- long-baseline (far from event) parallax-corrected heliocentric proper-motion magnitude **and direction** consistency with `.out` (`murel_helio`, `murel_helio_alpha`, `murel_helio_delta`)
- BJD consistency with `SIMULATION_ZERO_TIME + Simulation_time`
- noise sanity using normalized residual statistics `(observed - true) / sigma`

## BAGLE Joint-Fit Sanity Check

This check performs a full **combined photometric + astrometric** fit with BAGLE on one event that satisfies:
- `NSource == 1`
- `ObsGroup_0_chi2 < 100` (configurable)
- `ObsGroup_0_FiniteSourceflag == 0` (PSPL-like event)

It then:
- fits `PSPL_PhotAstrom_Par_Param1` using BAGLE nested sampling
- generates a diagnostic plot with best-fit model vs data
- uses event timing in MJD (explicit JD/BJD -> MJD conversion)
- compares fitted proper-motion vector (`muRel`) against `.out` `murel_helio_alpha/delta`
- compares fitted parallax vector (`piE`) against `.out` `piEE/piEN`
- fails with detailed diagnostics if amplitude or direction mismatches exceed thresholds
- applies documented sign conversion when comparing to `.out` vectors:
  - `.out` stores lens-source convention, BAGLE `muRel`/`piE` uses source-lens convention
  - comparison uses `(-murel_helio_alpha, -murel_helio_delta)` and `(-piEE, -piEN)`
- if PyMultiNest / MultiNest runtime is unavailable, the check automatically falls back to a scipy least-squares BAGLE fit

Run it standalone on pre-generated outputs:

```bash
python smoke_test/run_bagle_fit_sanity.py \
  smoke_test/output/general/smoke_general \
  --params smoke_test/parameterfiles/smoke_general.prm \
  --chi2-max 100 \
  --n-live-points 250
```

Typical dependencies:
- `bagle` (BAGLE_Microlensing)
- `joblib`
- `pymultinest` (+ MultiNest runtime), `dynesty`, `ultranest` (optional when scipy fallback is used)
- `matplotlib` (for the diagnostic plot)
- `scipy` (for fallback optimization path)

## Test Configuration

### Simulation Parameters

The smoke test uses simplified parameters for quick execution:
- **Simulation length**: 100 days
- **Lightcurve timeout**: 10 seconds (enforced on VBMicrolensing)
- **Single observatory**: "SmokeScope" with F184 filter
- **Exposure time**: 5 seconds
- **Astrometry**: Enabled with 0.1 mas systematic floor
- **6 synthetic sources and 6 synthetic lenses** per test

### Observatory Configuration

The test observatory (`smoke.observatory`) is configured with:
- **Filter**: F184 (index 6, ~2 μm)
- **Detector**: 10M electron full well capacity
- **Exposure**: 5 seconds per observation
- **Zero magnitude**: 20.0
- **Sequence**: 100 days of observations

## Expected Outputs

For each test case, the following outputs are generated:

### 1. Summary Reports (`.out`)
Located in `smoke_test/output/{std,croin,fish,general}/smoke_{std,croin,fish,general}/`

Contains run statistics, parameter values, and execution timing.

### 2. Lightcurve Files (`.lc`)
ASCII tables with columns including:
- `Simulation_time`: Observation epoch (days)
- `measured_relative_flux`: Photometric flux (relative to baseline)
- `measured_relative_flux_error`: Photometric error
- `true_relative_flux`: True simulated flux
- Astrometric centroids in **observer-centric ecliptic EN** milliarcseconds
- Sky positions in ICRS `RA/Dec` (degrees) generated from ecliptic EN via declared transform
- Explicit frame/convention headers:
  - `#Astrometry_Frame`
  - `#Astrometry_EventToEcl`
  - `#Astrometry_Transform`
  - `#Astrometry_BAGLE`
- Astrometric errors
- Parallax information
- Source and lens positions

### 3. Visualization Plots (`.png`)
If matplotlib is available, the script generates diagnostic plots showing:
- **Top-left**: Lightcurve (flux vs time) with measured points, true flux line, and baseline
- **Top-right**: Astrometric position in RA/Dec with error bars
- **Bottom-left**: Astrometric trajectory in North/East coordinates
- **Bottom-right**: Astrometric centroid shifts vs time

## Validation Criteria

The smoke test passes when:
1. All selected executables complete without errors (exit code 0)
2. At least one `.out` summary file is created and non-empty
3. At least one `.lc` lightcurve file is generated
4. No timeouts occur during execution

## Troubleshooting

### "Missing executable" error
- Ensure you've built the project: `cmake --build build`
- Check that executables exist in `build/bin/`

### "Missing ESPL.tbl" error
- Copy the VBMicrolensing ESPL table: `cp VBMicrolensing/ESPL.tbl src/`

### Timeout errors
- Increase timeout: `--exec-timeout 300`
- Check for infinite loops in lightcurve generation
- Verify VBMicrolensing `LC_TIMEOUT` is enforced (see `src/pllxLightcurveGenerator.cpp`)

### Plotting errors
- Install required packages: `pip install numpy pandas matplotlib`
- Run without plotting by uninstalling matplotlib temporarily

### "Negative error values" warning
- This was fixed by using pandas with named columns
- Verify column names match between code and output files

## Technical Details

### VBMicrolensing Timeout Enforcement
The smoke test relies on timeout enforcement in VBMicrolensing to prevent infinite hangs during lightcurve generation. The timeout is set via `LC_TIMEOUT=10.0` in parameter files and enforced in `src/pllxLightcurveGenerator.cpp`.

### Pandas-Based Column Access
The plotting function uses pandas DataFrames with named columns for robust data access. This prevents column index misalignment issues and makes the code more maintainable.

### Data Format
Lightcurve files use whitespace-delimited ASCII format with:
- Comment lines starting with `#`
- Header line with column names
- Data rows with numerical values

## Performance

Typical execution times (MacBook Pro, M1):
- **gulls_std**: ~10 seconds
- **gulls_croin**: ~10 seconds  
- **gullsFish**: ~10 seconds
- **gulls_general**: similar to gulls_std (varies)
- **Total** (all): varies with selected cases

## Contributing

When modifying the smoke test:
1. Keep simulation parameters minimal for fast execution
2. Ensure all executables are tested
3. Verify plots are generated correctly
4. Update this README if you add new features or change behavior
