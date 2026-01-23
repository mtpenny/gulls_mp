# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [3.1.0] - 2006-01-25 Astrometry Readdition



## [3.0.0] - 2026-01-15 Gulls General

> A new lightcurve generator has been developed to handle multi-lens and multi-source lightcurves with orbital motion - this is called general

### Limitations

- Can only handle up to binary stars as source or lens hosts, if higher multiples are provided then the behaviour is undefined (segfaults are likely)
- Can handle multiplanet systems or single-planet, multiple-moon systems in single or binary star systems. Can only put planets around one of the stars. No planets in source star systems
- Secondary sources can be dark

### Conventions
- A single star from the binary (source or lens) is chosen as the main lens (need not be the most massive or most bright). This is the host for wide-binary lenses. All relative quantities (e.g.,mass ratio) microlensing parameters (thetaE, rE, etc.) are relative to this object
- Orbits are defined by the vector addition of Keplerian ellipses, with objects internal to the currently considered object experiencing reflex orbital motion. It is the user's responsibility to ensure stability of the system in postprocessing
- Semimajor axis refers to the combined semimajor axis
- Input periods may be modified by the presence of additional masses

### Input formats

Binary stars are included via SynthPop's conventions for binary stars; either the primary or secondary star is drawn from the catalog, and its companion brought along with it. Physical parameters of the binary orbit are provided in the SynthPop catalog (period, eccentricity), others are free to be generated in gulls (currently eccentricity defaults to zero). Planets are supplied through planets files that list the following parameters:
```
Mass SemimajorAxis Eccentricity Inclination LongitudePerihelion LongitudeAscNode OrbitType
```
with angles in degrees. If inclination has a value >900, its inclination will be referenced to the binary star orbit's as `I = I_binary + (I-1000.0)`; if the star is not a binary, the inclination will be `I-1000`. OrbitType is an integer code, with a value of 1 or 2 indicating a planet, and 3 indicating a moon. Use multiple sets of these parameters for multiple planets/moons.


## [2.1.0] - 2025-12-27

### Added
- Basic astrometry support for microlensing simulations
- Six new output columns for centroid positions and uncertainties (true and observed)
- `ASTROMETRY_ON` parameter to enable/disable astrometric computations (default: 0)
- `ASTROMETRIC_SYS_FLOOR` parameter for systematic uncertainty floor in mas (default: 0.1)
- Comprehensive astrometry documentation (`documentation/source/astrometry.rst`)
- Centroid calculation using VBMicrolensing low-level functions (`BinaryMag2`, `MultiMag2`)
- Noise model implementation based on Gould & Yee (2014)
- Flux-weighted blending for centroids with lens and ambient stars
- Support for single and binary source configurations in astrometry

### Changed
- Variable naming: renamed `Asrc1`/`Asrc2` to `musrc1`/`musrc2` for consistency with magnification nomenclature
- Enhanced photometry module to compute astrometric uncertainties from photometric precision
- Updated lightcurve generators to compute and store true centroid positions

### Fixed
- Potential division by zero when θ_E (Einstein radius) is very small
- Added validation check with warning for events with θ_E < 1e-10 mas
- Spelling errors in comments: "oposite" → "opposite", "shif" → "shift", "abient" → "ambient"
- Buffer overflow in `snprintf` call (missing buffer size argument)
- Unbalanced parenthesis in documentation formula
- Trailing whitespace in parameter reading code
- Step numbering in astrometry documentation
- Pinned Sphinx version; `sphinx` and `sphinx-rtd-theme` became incompatable at version 7.

### Security
- Added validation to prevent division by zero in astrometric calculations
- Fixed buffer overflow vulnerability in string formatting


## [2.0.0] - 2025-10-20 Significant Quality of Life Updates

### Added
- Comprehensive documentation expansion with Sphinx/Read the Docs
- Input file format specifications and parameter reference
- Validation system for input catalogs and configuration files
- CI/CD pipeline with GitHub Actions
- Smoke test suite for automated testing
- CMake build system alongside traditional Makefile (includes `gulls_std`, `gulls_croin`, and `gullsFish`)
- Contributing guidelines and development workflow
- Example configurations and troubleshooting guides
- Missing docs build requirements in the `environment.yml`
- Smoke test output figures in Release Notes
- PSF generation utilities (`generateMoffat`, `precompute_psf`, `generateMoffatPSF`) to CMake build
- On-demand PSF file generation for smoke tests (no more 68MB files in git)
- Smart PSF caching system that reuses existing files when available

### Changed
- Improved error handling and user feedback
- Enhanced input validation with detailed error messages
- Streamlined installation process with better dependency management
- Updated documentation from legacy format to modern RST
- Simplified release workflow: patch > edit changelog > release
- PSF generation now uses proper subpixel sampling (9×9 = 81 variations)
- Smoke tests generate PSF files on-demand instead of requiring pre-committed files

### Fixed
- Buffer overflow issues in path handling
- Infinite loop bugs in random number generation
- Uninitialized memory issues in binary source calculations
- Off-by-one errors in catalog parsing
- Failure to build docs in the release workflow on GitHub
- PSF generation in CI environments (removed hardcoded local machine paths)
- PSF file size issues (now generates proper 68MB files with subpixel sampling)
- Simulation crashes due to missing or malformed PSF files

### Security
- Fixed potential buffer overflows in file path construction
- Improved input validation to prevent malformed data crashes


## [1.0.0] - 2013-2025 Historic Untracked Releases

### Added
- Core microlensing simulation framework
- Support for single and binary sources/lenses
- Realistic observing conditions and detector effects
- Multiple observatory and filter system support
- Photometric signal generation
- Lens orbital motion
- Support for low-level VBMicrolensing functions
- Detection statistics and survey planning tools
- Fisher uncertainty estimation
- Many other features, executables, and documentation evolutions

### Original Development
- Initial implementation by Matthew Penny and collaborators
- Published in Penny et al. (2013, 2014, 2019)
- Core algorithms for gravitational microlensing event simulation
