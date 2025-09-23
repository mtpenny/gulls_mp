# GULLS Contributors Guide (AGENTS)

This document orients new contributors to the GULLS microlensing simulator: what it does, how it’s structured, how to build and run it, and the key conventions (especially reference frames) that matter when modifying or extending the code.

If anything here disagrees with the parameter file documentation or source code, treat the code and parameter docs as authoritative.

---

## What GULLS is

GULLS simulates large populations of microlensing events by pairing sources and lenses drawn from Galactic population synthesis catalogs, then “observes” them with a configurable observatory model (sky tiling, cadence, filters, detector, and weather). It can produce both light curves and synthetic images and supports single lenses, binary lenses (incl. planet injections), and specialized runs (e.g., free-floating planets).

Key references: Penny et al. 2013, 2019.

---

## Architecture overview

At a high level, one run is orchestrated by a parameter (.prm) file that wires together:

- Catalog inputs
  - Source catalogs (limited to faint-end H ≤ ~25–27 depending on filter)
  - Lens catalogs (no mag limit; mass matters)
  - Starfield catalogs (non-lensing background)
- Observatory model
  - Observatory definition (.observatory): site/platform, filters, instrument
  - Field centers (.centres), observing sequence (.sequence)
  - Detector model (.detector), throughput (.throughput), and weather (.weather)
- Event generation
  - Event geometry and kinematics (u0, t0, tE, rE, theta_E, pi_E, rhos, alpha, murel, vt, gamma)
  - Optional planet parameters (Mp, a, inc, e, phase, q, s, period)
- Rendering/outputs
  - Lightcurves, images, diagnostics, and event catalogs with weights
  - Optional “reduce” step to HDF5 for downstream analysis

Core moving pieces (by concern):

- Executables (selected via Makefile target)
  - Examples: gulls, gullsSingle, gullsFFP, gullsFish, gullsFoM, gullsHZ, gullsBinaryStar, …
  - Each configures a different simulation mode and detection/diagnostic policy
- Physics / lensing backend
  - Uses VBMicrolensing (VBB) for accurate lensing calculations
  - ESPL.tbl lookup required by some single-lens calculations
- Scene/cadence/photometry
  - Observing cadence from .sequence; weather gates visibility on 0.25 day steps
  - Detector + throughput model govern SNR, PSF, saturation, noise, and zero-points
  - Photometry modes: ideal, aperture/PSF (configurable); image rendering optional
- IO and configuration
  - Plain-text parameter files and data catalogs
  - CFITSIO for image output and FITS handling

Supporting scripts (not C++):
- Launch utilities and reduction tools (e.g., scripts/reduce_gulls.py) to convert raw run outputs to HDF5 and apply detection cuts post hoc.

---

## Reference frames and conventions (important)

Understand these before changing kinematics, parallax, or photometry:

Coordinates and positions:
- Galactic (l, b) in degrees: primary frame for event rates, kinematics, and catalog selection.
- Equatorial (RA, Dec): used for sky localization and some outputs; conversion typically originates from the catalogs.
- Field centers are provided per observatory; events are placed within each field’s footprint.

Distances and scales:
- Distances: Ds (source), Dl (lens) in kpc (from catalogs).
- Einstein radius: rE (projected physical scale; see outputs), theta_E (angular, mas).
- Source radius in Einstein units: rhos (dimensionless).

Velocities and proper motions:
- Proper motions provided in Galactic components:
  - Source: smu_l, smu_b (mas/yr)
  - Lens: lmu_l, lmu_b (mas/yr)
- Relative proper motion murel (mas/yr); transverse speed vt (km/s) derived consistently with Ds, Dl.
- Event trajectory angle alpha (radians) defined in the lens plane (check the executable’s convention before changing).

Time and ephemerides:
- Times in Julian Days (JD): t0 is peak time, tE is Einstein timescale (days).
- SIMULATION_ZERO_TIME and SIMULATION_LENGTH govern absolute scheduling.
- Weather is sampled every 0.25 day; cadence is set by the .sequence file.
- Parallax:
  - pi_E (dimensionless microlens parallax) optionally enabled via PARALLAX=1.
  - Platform/observer geometry is governed by SPACE and ORBIT in .observatory and (implicitly) by the cadence.

Photometric system:
- Magnitudes are Vega unless specified.
- Zeropoint: DETECTOR ZEROMAG and ZEROFLUX define the instrument flux scale.
- Background in mag/arcsec^2 (BACKGROUND or SKY_BACKGROUND).
- Filters are enumerated integers that map to columns in the star catalogs (NFILTERS must match).

Detector/PSF:
- PIXELSCALE (arcsec/pixel), PSFFWHM (arcsec), PSFFILE kernel with subpixel placement (SUBPIX).
- Full-well, gain, read noise, dark current, and systematics affect SNR and saturation.
- PRETTY_PICS toggles image generation for diagnostics.

Units summary (typical; verify per output schema):
- l, b, RA, Dec: degrees
- Ds, Dl: kpc
- proper motions: mas/yr
- vt: km/s
- t0, tE: days (t0 in JD; tE duration)
- theta_E: mas
- rhos, u0, pi_E, q, s: dimensionless
- magnitudes: Vega
- fluxes: instrument-dependent counts/sec

---

## Data flow

1) Select catalogs and observatory
- Source/lens/starfield lists and files
- Observatory (fields, sequence), detector, throughput, weather
- Optional planets catalog(s)

2) Configure a .prm file
- Names the executable, input directories/lists, run timing, and output controls
- Enables diagnostics (PRETTY_PICS, OUTPUT_* flags)

3) Build the executable
- Makefile compiles a chosen target (e.g., gullsSingle)

4) Run
- bin/<executable>.x -i <paramfile.prm> -s <instance> [-f <field>] [-d …]
- Outputs to OUTPUT_DIR; FINAL_DIR for post-run organization

5) Reduce (optional but recommended)
- scripts/reduce_gulls.py to produce HDF5, normalize weights, and apply detection cuts downstream

---

## Inputs (by file type)

- .observatory: platform/location, filters, detector/throughput, pointing, cadence, weather
- .sequence: per-visit timing, stacks, interleaving, repeats
- .detector: instrument noise, PSF, pixel scale, saturation, zeropoints
- .throughput: passband + system response
- .weather: observing toggles per quarter-day
- Star catalogs:
  - sources/: typically constrained faint-end; require sufficient density (~1e5 per field)
  - lenses/: no magnitude limit; suggest ~1e4 per solid angle ~1e-4 deg^2
  - starfields/: non-lensing background split by brightness bins (H ≤15, 15–20, 20–25, >25)
- planets/: optional planet parameters for binary-lens simulations
- rates/: lensing rate model (e.g., bH.pmcorrected.rates)

---

## Outputs

- Event catalogs with “|”-delimited sections:
  - Field/coordinates, source properties, lens properties, event parameters, planet parameters, weights/flags, magnitudes per filter, blending fs, chi-square diagnostics
- Lightcurves:
  - .det.lc for detected events; .all.lc for all generated (if enabled)
- Images:
  - Baseline and peak per filter (if OUTPUT_IMAGES=1)
- Reduced HDF5 (via reduction script):
  - *_out*.hdf5: all events
  - *_det*.hdf5: detected subset (or apply cuts later)

Use final weight (w) for statistical analyses.

---

## Build and run (macOS quickstart)

Dependencies:
- GSL
- CFITSIO
- VBMicrolensing (VBB library)
- C++ toolchain (g++; some builds prefer GCC over Clang)

Install (Homebrew examples):
- brew install gsl cfitsio gcc
- Build VBMicrolensing and ensure lib path is available to linker
- Obtain and place licensed files:
  - headers/random.h, headers/zroots2.h
  - classes/random.cpp, classes/zroots2.cpp
- Copy ESPL.tbl from VBMicrolensing/VBMicrolensing/data to src/

Environment:
- export GULLS_BASE_DIR="/path/to/gulls/"
- Optional: export GULLS_STARS_DIR="/path/to/catalogs/"

Configure and build:
- ./configure.sh
- cd src
- make gullsSingle
  - Adjust Makefile BASEDIR, CFLAGS/CPPFLAGS include paths (VBB, CFITSIO)
  - Adjust LINKERFLAGS with -L (VBB, CFITSIO) and -lVBB -lcfitsio -lgsl -lgslcblas -lgfortran -lstdc++

Run (example):
- cd bin
- ./gullsSingle.x -i PATH/singlelens.prm -s 0 -f 83 -d

Diagnostics:
- For first runs: enable
  - PRETTY_PICS=1
  - OUTPUT_LC=1
  - OUTPUT_IMAGES=1
  - OUTPUT_ONALL=1

---

## Common pitfalls

- “White square” images
  - Use zscale in DS9
  - Too-bright stars or bad magnitudes; try larger PRETTY_PICS_DIMENSIONS and inspect catalogs
- “No valid stars”
  - NFILTERS mismatch vs. catalog filter columns
- Link errors
  - Missing VBB or CFITSIO include/lib paths; missing ESPL.tbl or missing random/zroots files
- Commenting in .prm
  - The parser is string-based; “commenting out” by prefixing with ‘#’ may not prevent overrides if the key still appears on the line. Remove or alter the key name to disable.

---

## Contributing

- Branching and PRs
  - Create feature branches, open PRs with a short design note if the change touches physics, reference frames, or file formats
- Code style
  - C++ (g++); prefer modern, readable code; follow existing patterns in src/
  - Guard platform-specific code; keep defaults portable
- Tests and validation
  - Add small synthetic parameter files and a short CI-like run (few events) where possible
  - Validate reference-frame changes with a known scene (fixed catalogs, deterministic RANDOM_SEED)
- Documentation
  - Update documentation/source/*.rst or top-level docs when changing formats or parameters
  - Note units and frames explicitly
- Performance
  - Avoid regressions in the tight loops (event generation, photometry); prefer pre-allocation and const references
- Backwards compatibility
  - Avoid breaking existing parameter files; add new keys with sensible defaults

---

## Where to look in the code (orientation)

- src/
  - Main executables wiring (e.g., gullsSingle, gullsFFP, gullsFish)
  - Event setup, reading .prm (e.g., readParamfile.cpp), planet handling (e.g., standardPlanet.cpp)
  - Photometry and imaging (e.g., photometry.cpp/.h), detector modeling
  - Catalog readers and IO (CFITSIO interactions)
- headers/ and classes/
  - Utility classes, math helpers, random/zroots (if provided)
- scripts/
  - reduce_gulls.py and launch helpers

Names may vary slightly across branches; search for the above to locate relevant files.

---

## Minimal checklist for a new feature

- Define the scope and which executable(s) it touches
- Confirm units and reference frames; write them down in code comments
- Extend .prm parsing if needed; provide defaults that don’t break older runs
- Add a tiny parameter file exercising the new path (runs in seconds)
- Update docs and AGENTS.md section(s)
- Run a smoke test:
  - PRETTY_PICS=1, few subruns, OUTPUT_* on
- Sanity-check outputs (units, distributions, flags)

---

## Useful commands

- macOS include/lib discovery (Homebrew x86 vs Apple Silicon):
  - Includes: /usr/local/include or /opt/homebrew/include
  - Libs: /usr/local/lib or /opt/homebrew/lib
- Example Make overrides:
  - Add to CFLAGS/CPPFLAGS: -I/path/to/VBMicrolensing -I/path/to/cfitsio/include
  - Add to LINKERFLAGS: -L/path/to/VBMicrolensing -L/path/to/cfitsio/lib -lVBB -lcfitsio

---

## Questions or issues

- Open an issue: https://github.com/gulls-microlensing/gulls/issues
- Include:
  - OS and compiler (gcc/clang versions)
  - The .prm file and any custom observatory/detector files
  - Build command and full error/output logs