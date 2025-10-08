# VBMicrolensing astrometry integration plan

## Why this document exists
The VBM (VBMicrolensing) library already computes light curves *and* astrometric centroids, but the current GULLS astrometry rewrite has placeholder hooks (`pllxLightcurveGenerator.cpp`) that still need real numbers. This plan records the pieces of the upstream VBM stack we rely on, where they live in the source tree, and what we have to reproduce or call into from our own code.

Sources that were actually inspected are called out inline; everything is in the local clone at `/Users/malpas.1/Code/VBMicrolensing` unless noted otherwise.

---

## Pre-requisites before you can ask VBM for astrometry

| Step | What | Where to read | Notes |
| --- | --- | --- | --- |
| 1 | Instantiate `VBMicrolensing.VBMicrolensing()` | `VBMicrolensing/__init__.py` | Constructor wires default tables (`ESPL.tbl`, `SunEphemeris.txt`) via `SetESPLtablefile` and `SetSuntablefile`. |
| 2 | Point VBM at the event coordinates | `docs/python/Parallax.md`, §Target coordinates | Call `VBM.SetObjectCoordinates("RA DEC")` or `SetObjectCoordinates(coord_file, satellite_dir)`. Required before any parallax or astrometric call; the python bindings refuse to run otherwise. |
| 3 | (Optional) Load custom ESPL/Sun tables | `VBMicrolensing/lib/VBMicrolensingLibrary.cpp` (`LoadESPLTable`, `LoadSunTable`), docs linked above | Only needed if our date range falls outside the default 1990–2050 tables under `VBMicrolensing/data/`. |
| 4 | Choose parallax conventions | `docs/python/Parallax.md` | Defaults are HJD′ times and North/East parallax components. Override with `VBM.t_in_HJD = False` or `VBM.parallaxsystem = 0` if we want JD′ or $(\u03c0_{E,\parallel}, \pi_{E,\perp})`. |
| 5 | Enable/disable luminous companions | `docs/python/CentroidTrajectories.md` and `VBMicrolensing/lib/VBMicrolensingLibrary.cpp` (`BinaryAstroLightCurve*`) | `VBM.turn_off_secondary_lens = True` forces the companion dark. `VBM.lens_mass_luminosity_exponent` defaults to 4 for the flux scaling $F \propto M^4$. |

---

## Astrometry-enabled API surface
See `VBMicrolensing/lib/python_bindings.cpp` lines 1200–1440 for the authoritative parameter ordering and return layouts.

| Function | Parameter vector (python binding order) | Output arrays | Typical use |
| --- | --- | --- | --- |
| `PSPLAstroLightCurve` | `[u0, log_tE, t0, \pi_N, \pi_E, \mu_{S,N}, \mu_{S,E}, \pi_S, \theta_E]` | `[mag, c1_src, c2_src, c1_lens, c2_lens, y1, y2]` | Single lens, finite-source optional via `ESPL` tables. |
| `ESPLAstroLightCurve` | `[u0, log_tE, t0, log_rho, \pi_N, \pi_E, \mu_{S,N}, \mu_{S,E}, \pi_S, \theta_E]` | Same as PSPL | Single lens with limb darkening lookup. |
| `BinaryAstroLightCurve` | `[log_s, log_q, u0, \alpha, log_\rho, log_tE, t0, \pi_N, \pi_E, \mu_{S,N}, \mu_{S,E}, \pi_S, \theta_E]` | Same as PSPL | Static binary lens. Flux ratio applied via $q^{p}$ with `p = lens_mass_luminosity_exponent`. |
| `BinaryAstroLightCurveOrbital` | Binary parameters + `[w1, w2, w3]` | `[mag, c1_src, c2_src, c1_lens, c2_lens, y1, y2, separation]` | Circular orbital motion; angular velocity vector $(w1,w2,w3)$ in Einstein angles/day. |
| `BinaryAstroLightCurveKepler` | Orbital version + `[s_{z/s}, a_{r}]` | Same as orbital | Full Keplerian evolution. |
| `BinSourceAstroLightCurveXallarap` | See binding for full list (two-source parameters + orbital terms) | `[mag, c1_src, c2_src, c1_lens, c2_lens, y1_s1, y2_s1, y1_s2, y2_s2]` | Binary source with xallarap. |
| `BinSourceBinLensAstroLightCurve` | Binary lens + binary source + separation array | Adds `y*_s2` and `separation` arrays | Handles two luminous sources around a binary lens. |
| `TripleAstroLightCurve` | Triple-lens parameterization | `[mag, c1_src, c2_src, c1_lens, c2_lens, y1, y2]` | Triple systems with parallax. |
| `CombineCentroids` | Takes light-curve magnitudes, centroid arrays + blend ratio `g = F_L/F_S` | `[c1_combined, c2_combined]` | Flux-weighted centroid mix; see `VBMicrolensing/lib/VBMicrolensingLibrary.cpp` around line 4895. |

Reference walkthrough: `docs/python/CentroidTrajectories.md` demonstrates the PSPL case end-to-end and shows centroid blending via `CombineCentroids`.

---

## What PR #3 (mtpenny/gulls_mp) already does — and where it diverges

The draft astrometry pipeline in PR #3 is a useful reference point:

* ✅ **Centroid source** – Calls `BinaryAstroLightCurve` to get sky-frame centroids (mas) for each observatory/time grid.
* ✅ **Noise model** – Adds independent Gaussian noise with
  $\sigma_{\text{1D}} = \text{FWHM} / \text{SNR}$ per axis, combines with a configurable systematic floor (default 0.1 mas), and injects the result into the stored centroids.
* ✅ **Lens motion + RA/Dec export** – Converts Galactic lens proper motions (MUL/MUB) to equatorial, propagates from $t_0$ in years, and writes both offsets (mas) and absolute coordinates (deg) with wrap-around handling.
* ⚠️ **Frame conversions** – Parallax and source proper motions are still fed through the original “Step 1” approximations (piEN → piN, MUL/MUB → μRA/μDec) with a TODO warning; we must supply the full `coords::mulb2ad`/`coords::muecl2ad` conversions when we finish the port.
* ⚠️ **High-level dependency** – The entire centroid stream comes from the `*AstroLightCurve` helpers. Matt asked for a lower-level implementation that mirrors VBM’s internals, so we should treat that code as a validation scaffold only.

Action items for our branch:

1. Keep the **noise+systematics recipe** from PR #3 (FWHM/SNR ⊕ `ASTROMETRIC_SYS_FLOOR`) and the per-axis Gaussian perturbation.
2. Port the **coordinate transforms** for parallax and source μ into the low-level path so we can retire the Step 1 shortcuts.
3. Replace the high-level `BinaryAstroLightCurve` call with our own reproduction of `ComputeCentroids`, fed by the existing per-epoch `BinaryMag2` evaluations.

The notes below spell out how to achieve that step-by-step.

---

## Inside `VBMicrolensingLibrary.cpp`: what the centroids actually mean

Key functions (all around lines 4860–5180 of `VBMicrolensing/lib/VBMicrolensingLibrary.cpp`):

- **`ComputeCentroids`** transforms the instantaneous lens-plane centroid (`c1s`, `c2s`) into sky coordinates.
  - Inputs expect Einstein units along lens axes; the routine multiplies by `thetaE` to obtain mas.
  - It derives lens proper motions from the source motion plus parallax (`pai1`, `pai2`) and keeps everything in mas/day.
  - Rotation from the lens frame to North/East uses `PosAng = atan2(pai2, pai1) - alpha + dPosAng` with adjustments for orbital motion (`dPosAng`).
  - Lens position includes parallax offsets via the stored ephemeris vectors `Et` and `Ehel`.
- **`BinaryAstroLightCurve*`** wrappers populate `astrox1` / `astrox2` before calling `ComputeCentroids`. Those members carry the Einstein-frame centroid in the original x1/x2 system.
- **`CombineCentroids`** simply flux-weights the source and lens centroids: `c_tot = (mag * c_src + g * c_lens)/(mag + g)`.

The limb-darkening accumulator `annulus::LDastrox1` (header definition at `VBMicrolensing/lib/VBMicrolensingLibrary.h`, lines ~170–190) stores the centroid contribution per annulus while `ESPLMagDark` refines the integral.

These details matter if we re-implement the centroid math inside GULLS: we must match the same rotation, time conventions (HJD′ by default), and parallax corrections.

---

## Proposed workflow for GULLS

1. **Leverage the Python façade for validation**
   - Use the python bindings documented above to dump centroid trajectories for the exact lens/source configurations we simulate in C++.
   - Archive those arrays under `Astrometry/validation/` for regression diffing. (`centroid.ipynb` already exists here.)

2. **Mirror (and then _own_) the centroid math in `pllxLightcurveGenerator.cpp`**
  - Continue using the existing per-epoch binary magnification path (`BinaryMag2`, etc.) so we control the sampling.
  - After each magnification call, grab `VBM.astrox1/astrox2` (Einstein units) and run our local copy of the `ComputeCentroids` recipe:
    1. Multiply by `thetaE` to convert to mas.
    2. Rotate from lens axes into the sky frame with `PosAng = atan2(pai2, pai1) - alpha + dPosAng`.
    3. Apply source proper motion and parallax offsets using the preloaded `Et`, `Ehel`, and the equatorial parallax components.
    4. Add the lens centroid and optional flux-weighted lens contributions (`q^{p}` vs `turn_off_secondary_lens`).
  - This gives us c1/c2 without touching the `*AstroLightCurve` stack while still guaranteeing parity with VBM internals.

3. **Reuse VBM tables where possible**
   - The ESPL lookup cubes live in `VBMicrolensing/data/ESPL.tbl` and are loaded on demand. If we need native C++ access, lift `LoadESPLTable` (lines ~420–470) rather than reinventing the interpolation.
   - Sun ephemeris default range is 1990–2050 (see `VBMicrolensing/data/SunEphemeris.txt`). If future events fall outside, regenerate via NASA Horizons and update the data path before runtime.

4. **Handling luminous vs. dark companions**
   - Respect `turn_off_secondary_lens` when mixing lens centroids: the flux correction happens right after `ComputeCentroids` in each binary astrometry routine. Replicate the same `FR = q^{lens_mass_luminosity_exponent}` scaling.

5. **Testing hooks**
   - Write a thin C++ harness that calls `BinaryAstroLightCurve` via the existing shared library; compare against our hand-rolled centroid outputs to ensure sub-microarcsecond agreement before swapping implementations.

---

## Low-level centroid + noise checklist (implementation crib sheet)

1. **Build per-epoch geometry**
  - Use the same absolute JD array you pass into photometry.
  - Convert parallax components (ecliptic → equatorial) with `coords::muecl2ad` and log-scale inputs back to linear where needed.
  - Convert source proper motions (Galactic → equatorial) with `coords::mulb2ad`.

2. **Magnification + Einstein-plane centroid**
  - Call `BinaryMag2(s, q, y1, y2, rho)` (or the appropriate variant) so `VBM.astrox1/astrox2` are populated for that epoch.
  - Cache `VBM.astrometry = true` so the internal bookkeeping is active.

3. **Reimplement `ComputeCentroids` locally**
  - Coefficients to port: the parallax vectors (`Et`, `Ehel`), the rotation angle `PosAng`, and the lens proper motion terms `muL1`, `muL2` obtained from source μ, $\theta_E$, $\pi_E$, and the ephemerides.
  - Output arrays: source centroid N/E (`c1s`, `c2s`), lens centroid N/E (`c1l`, `c2l`), plus the combined centroid when blending.

4. **Lens proper motion + absolute coordinates**
  - Propagate the lens catalog proper motion from $t_0$ (years) and add to the centroid offsets before converting to RA/Dec.
  - Convert $\Delta\mathrm{Dec} = N / (3.6\times10^6)$, $\Delta\mathrm{RA} = E / (3.6\times10^6 \cos\delta_0)$ and wrap RA into [0, 360).

5. **Noise + systematics**
  - Photon term: `σ_{\text{photon}}[\text{mas}] = 1000 × \text{FWHM}[\text{arcsec}] / \text{SNR}`.
  - Floor: read `ASTROMETRIC_SYS_FLOOR` (mas) from the parameter file, add in quadrature.
  - Randomization: draw two independent `𝒩(0, σ)` deviates (Box–Muller is fine) and perturb the North/East offsets.
  - Store both perturbed centroids and their uncertainties in mas and in degrees.

6. **Diagnostics**
  - Keep PR #3’s VBM round-trip mode to cross-check the new implementation until we hit parity.
  - Log when the VBM numerical error estimate `≈ 50 ρ Tol θ_E` exceeds the nominal noise budget.

---

## Follow-up questions / open items

- **Parallax extrapolation warnings**: The python bindings print a warning if we ask for times outside the ephemeris range (`parallaxextrapolation > 0`). Decide whether GULLS should surface the same warning or fail hard.
- **Coordinate system sanity**: Double-check that our upstream pipeline always passes parallax in the North/East frame. If not, expose `parallaxsystem` as a configurable knob alongside other event metadata.
- **Data packaging**: consider vendoring `ESPL.tbl` and `SunEphemeris.txt` into the Astrometry repo so we are not implicitly relying on the installed VBM wheel.
- **Retire the high-level dependency**: Once the low-level path matches VBM to <1 μas, delete (or guard) the `*AstroLightCurve` shortcut so we do not regress.

---

## Pointers for future digging

- `VBMicrolensing/lib/VBMicrolensingLibrary.cpp`
  - Centroid machinery: lines ~4860–5200
  - Limb-darkening integration: lines ~430–700 (for ESPL)
- `VBMicrolensing/lib/python_bindings.cpp`
  - Parameter specs and return order: lines ~1200–1440
- `VBMicrolensing/docs/python/CentroidTrajectories.md`
  - Narrative walkthrough + plotting recipes
- `VBMicrolensing/docs/python/Parallax.md`
  - Coordinate setup, satellite ephemerides, parallax conventions
- `VBMicrolensing/data/`
  - Default ESPL table, Sun ephemeris, example satellite tables
- GULLS staging area
  - `gulls_mp/src/classes/VBMicrolensingLibrary.cpp`: local copy of upstream routines (search for `astrox1`)
  - `Astrometry/validation/centroid.ipynb`: existing notebook we can extend for regression tests

Document owner: Copilot (this file). Update as you discover more quirks or introduce wrapper utilities.

---

## Implementation Overview

The placeholder structure is already in place. We need to:

1. **Enable VBM astrometry** and extract centroids after each `BinaryMag2` call
2. **Transform centroids** from Einstein units in the lens frame to mas in the sky frame  
3. **Add noise model** to generate observed centroids with realistic uncertainties
4. **Implement coordinate transforms** for parallax and proper motion conversions

## Detailed Implementation Plan

### Step 1: Extract VBM centroids in pllxLightcurveGenerator.cpp

After the existing `BinaryMag2` calls (lines ~84 for single source, ~100 for binary source), we need to:

1. Enable astrometry before the first epoch: `Event->vbm->astrometry = true;`
2. After each `BinaryMag2` call, grab `Event->vbm->astrox1` and `Event->vbm->astrox2`
3. Transform from Einstein units to mas: multiply by `Event->thE`
4. Rotate from lens frame to sky frame using the parallax direction
5. Store in `Event->xctrue[idx]` and `Event->yctrue[idx]`

The rotation angle involves:
- `alpha` (trajectory angle in the lens plane)
- Parallax direction: `atan2(piE_E, piE_N)` where piE components come from the parallax structure
- The centroid shift is perpendicular to the source trajectory

### Step 2: Implement noise model in `photometry.cpp`

Currently lines 98-110 have placeholders. We need to:

1. Calculate photon noise per axis: `σ_photon[mas] = 1000 × FWHM[arcsec] / SNR`
2. Add systematic floor: `σ_total = sqrt(σ_photon² + σ_sys²)` where σ_sys ≈ 0.1 mas
3. Draw Gaussian noise using existing `gasdev()` function
4. Apply noise: `Event->xc[idx] = Event->xctrue[idx] + Event->xcerr[idx] * gasdev(...)`
5. Store the uncertainty values in the `*err` vectors

### Step 3: Add coordinate transformation functions

The code needs functions to convert:
- Galactic proper motion (μ_l, μ_b) → Equatorial (μ_α, μ_δ)
- Ecliptic parallax (π_EN, π_EE) → Equatorial (π_N, π_E)

These should go in astroFns.cpp/astroFns.h following the existing `eq2gal`, `gal2eclip` pattern.

### Key Challenges

1. **Frame conventions**: VBM centroids are in the lens barycenter frame. We need to transform to:
   - Sky North/East frame (for mas offsets)
   - Eventually RA/Dec (for absolute coordinates)

2. **Parallax integration**: The `parallax` class already computes `piEN`, `piEE` but we need to ensure it's being used correctly for the centroid rotation.

3. **Binary source handling**: Line 107-108 placeholder is for binary sources (xallarap). This adds complexity because we need to flux-weight the two source centroids.

4. **Lens luminosity**: If the lens is luminous, we need to add its centroid contribution weighted by the flux ratio (see `CombineCentroids` logic from VBM).

## Files to Edit

- `outputLightcurve.cpp`
- `photometry.cpp`
- `pllxLightcurveGenerator.cpp`
- `structures.h`
- `timeSequencer.cpp`