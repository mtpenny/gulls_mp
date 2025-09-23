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

Understand these before changing kinematics, parallax, or photometry. The notes below reflect the current source code behavior.

Coordinates and positions:
- Galactic (l, b) in degrees: primary frame for event rates, kinematics, and catalog selection. Stored and output in degrees (`structures.h`, `buildEvent.cpp:377`).
- Equatorial (RA, Dec): used for sky localization and outputs. Internally in radians; written to outputs in degrees as `ra_deg`, `dec_deg` (`info.cpp:72`, `:214`).
- Conversions use `eq2gal`/`eq2horiz` in radians; in code the `'g'` flag to `eq2gal` converts Galactic → Equatorial (`timeSequencer.cpp:25`).

Distances and scales:
- Distances: Ds (source), Dl (lens) in kpc (from catalogs; see indices in `structures.h:404+`).
- Einstein radius rE in AU and thetaE in mas. Specifically: `rE = rEsun * sqrt(M⋅Ds⋅x(1−x))` with `rEsun=2.85412 AU` (M in Msun, Ds in kpc), and `thetaE = rE / Dl` (mas because AU/kpc = mas) (`headers/constants.h:26`, `buildEvent.cpp:512–516`).
- Source size in Einstein units: `rho` stored as `rs` (dimensionless) (`buildEvent.cpp:545`).

Velocities and proper motions:
- Catalog proper motions are Galactic components in mas/yr: `MUL`, `MUB` for each star (`structures.h:404+`).
- Event stores the heliocentric relative proper motion magnitude `murel` and components `murel_l`, `murel_b` in mas/yr (`buildEvent.cpp:526–540`).
- Transverse speed: `vt = 4.74047 × murel(mas/yr) × Dl(kpc)` km/s is implemented equivalently as `murel * Dl * AU / 1000 / SECINYR` (`buildEvent.cpp:536`).
- Trajectory angle `alpha` is in degrees (converted to radians where needed) (`buildEvent.cpp:500`, `fisher.cpp:46,109`).

Time and ephemerides:
- Epoch arrays (`jdtimes`) are in JD; however `t0` is stored and output in days relative to `SIMULATION_ZERO_TIME` (not an absolute JD) (`info.cpp:280`, `structures.h:117`).
- `tE_r` is the reference-frame Einstein timescale (days) and `tE_h` is heliocentric (days) (`info.cpp:291`, `classes/parallax.cpp:246,283`).
- Weather is sampled every 0.25 day; cadence comes from the `.sequence` file (`timeSequencer.cpp`).

Parallax:
- `piE` is dimensionless. The code tracks components `piEN`, `piEE` in the event’s reference frame using the ecliptic North/East basis; `piEll` and `piErp` are components // to and ⟂ to the Sun’s acceleration vector (`info.cpp:95–101`, `classes/parallax.cpp:312–334`).
- Parallax calculation and light-curve shifts are enabled when `pllxMultiplyer > 0` (commonly 1) and disabled when 0. This governs inclusion of `piEN/piEE` in fitting and shift application (`structures.h:143`, `fisher.cpp:70–100`).

Photometric system:
- Magnitudes are Vega unless specified.
- Zeropoint: detector `ZEROMAG`/`ZEROFLUX` define flux scale.
- Background in mag/arcsec^2 (`SKY_BACKGROUND` and zodiacal model).
- Filters are enumerated and mapped to catalog columns; `NFILTERS` must match.

Detector/PSF:
- `PIXELSCALE` (arcsec/pixel), `PSFFWHM` (arcsec), optional PSF kernel with subpixel placement.
- Full-well, gain, read noise, dark current, and systematics affect SNR and saturation.
- `PRETTY_PICS` toggles image generation for diagnostics.

Units summary (outputs; verify per schema):
- l, b: degrees; RA, Dec: degrees (internal radians)
- Ds, Dl: kpc; rE: AU; thetaE: mas; rho (`rs`): dimensionless
- Proper motions (murel, components): mas/yr; vt: km/s
- t0, tE: days (`t0` relative to `SIMULATION_ZERO_TIME`)
- u0, piE, q, s: dimensionless
- Magnitudes: Vega; Fluxes: instrument-dependent counts/sec

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

## Output Columns & Units (appendix)

Key columns written by the event catalog (see `gulls_mp/src/info.cpp`). Units and frames reflect the code’s current behavior.

- Position
  - `galactic_l`, `galactic_b` [deg]
  - `ra_deg`, `dec_deg` [deg] (internal RA/Dec are radians)

- Microlensing geometry
  - `u0lens1` [—]
  - `alpha` [deg] (converted to radians in computations)
  - `t0lens1` [days] relative to `SIMULATION_ZERO_TIME`
  - `tref` [days] reference epoch for parallax
  - `tcroin` [days], `ucroin` [—], `rcroin` [Einstein radii] (binary-lens parametrization; present when used)
  - `tE_ref`, `tE_helio` [days] (reference-frame vs heliocentric)
  - `rE` [AU]
  - `thetaE` [mas]
  - `rho` (aka `rs`) [—]

- Parallax
  - `piE` [—]
  - `piEN`, `piEE` [—] (components in ecliptic North/East of the reference frame)
  - `piEll`, `piErp` [—] (// and ⟂ to Sun acceleration in the reference frame)

- Proper motion (relative lens–source)
  - Heliocentric components: `murel_helio_alpha`, `murel_helio_delta` [mas/yr] (equatorial), `murel_helio_l`, `murel_helio_b` [mas/yr] (Galactic), `murel_helio_lambda`, `murel_helio_beta` [mas/yr] (ecliptic), `murel_helio` [mas/yr]
  - Reference-frame components: `murel_ref_alpha`, `murel_ref_delta`, `murel_ref_l`, `murel_ref_b`, `murel_ref_lambda`, `murel_ref_beta` [mas/yr], `murel_ref` [mas/yr]

- Velocities
  - `vtilde_helio`, `vtilde_ref`, `v_ref` [km/s] (projected lens, reference frame)
  - Component velocities: `vtilde_helio_N`, `vtilde_helio_E`, `vtilde_ref_N`, `vtilde_ref_E`, `v_ref_N`, `v_ref_E` [km/s]
  - `vt` [km/s] (lens transverse speed)

- Limb darkening
  - `LDgamma` [—]

- Planet parameters
  - `Planet_*` columns: names originate from the planet module; typical derived values include mass ratio `q` [—], separation `s` [Einstein radii], and period [days] when available.

- Weights and controls
  - `u0max` [—], `t0range` [days], `weight_scale` [—]
  - `raw_weight`, `weight` [—]

- Photometry and blending
  - Source and lens magnitudes per filter: `Source_<FILTER>`, `Lens_<FILTER>` [mag, Vega]
  - Blending per observatory: `Obs_<i>_fs` [—]
  - For multiple sources: `Source2_rho` [—], `Source2_s` [Einstein radii], `Source2_alpha` [deg], `Source2_inc` [deg], `Source2_phase` [—]; `Obs_<i>_fs2ofs1` [—]

- Diagnostics and groups
  - `NumObsGroups` [int], `ErrorFlag` [int]
  - Per group: `ObsGroup_<g>_flatlc` [int], `ObsGroup_<g>_flatchi2` [—], `ObsGroup_<g>_FiniteSourceflag` [int], `ObsGroup_<g>_chi2` [—] plus group-specific metrics
  - Optional error scaling: `scale_factor_300`, `scale_factor_n3sig3`, `scale_factor_n3sig6` [—]
  - `LCOutput` [0/1]

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
