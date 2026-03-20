# PR: Astrometry output feature for the `general` executable

## Summary

This PR adds astrometric centroid computation and sky-coordinate output to the
`general` gulls executable, validated against the BAGLE microlensing fitting
framework.  A sign error in the parallax perpendicular-shift formula
(`parallax.cpp`) was identified and fixed as part of this work.

## What changed

### New feature: astrometric outputs (multiple commits)

The `general` executable now computes per-epoch flux-weighted centroids for
lensed source images, blends them with luminous lens light, converts the result
to absolute sky coordinates (ecliptic lambda/beta, then RA/Dec), adds
photon-noise--limited astrometric noise, and writes everything to the `.lc`
lightcurve file.

### Bug fix: `ushift` sign in `parallax.cpp` (commit `1876c52`)

**File:** `src/classes/parallax.cpp`, function `compute_tushifts`, line 152.

**Before:**

```cpp
ushift[i] = -piE * (-Nshift[i]*sn + Eshift[i]*cs);
```

**After:**

```cpp
ushift[i] =  piE * (-Nshift[i]*sn + Eshift[i]*cs);
```

The sign of the overall `piE` factor on `ushift` was changed from negative to
positive.  `tshift` was not changed.

### Bug fix: `u0` sign convention in BAGLE comparison (same commit)

**File:** `smoke_test/bagle_forward_sanity.py`, line 349.

**Before:**

```python
u0_for_bagle = float(row["u0lens1"])
```

**After:**

```python
u0_for_bagle = -float(row["u0lens1"])
```

---

## Why the `ushift` sign change is correct: mathematical proof

### The two rotations

The parallax perturbation flows through two successive 2D rotations:

1. **`compute_tushifts`** (parallax.cpp): rotates the ecliptic observer offset
   `(Nshift, Eshift)` into the `(tau, u)` microlensing frame using angle
   `phi_pi`.
2. **`omLightcurveGenerator`** (omLightcurveGenerator.cpp, line 721): rotates
   `(tau, u)` back to the event frame `(xs0, ys0)` using angle `alpha`, where
   `xs0` is ecliptic East and `ys0` is ecliptic North.

For astrometry to be correct, the round-trip `(N, E) → (tau, u) → (E, N)` must
recover the physical parallax shift without distortion.

### Angle definitions

From `provide_murel_h_ad` (parallax.cpp):

```
phi_pi = atan2(piEE, piEN) = atan2(mulam_r, mubet_r)
```

This is the angle of the pi_E vector measured from ecliptic North toward
ecliptic East.

From `omLightcurveGenerator.cpp` (line 165):

```
alpha = atan2(-mubet_r, -mulam_r)
```

This is the direction of source motion relative to the lens (= negative of the
lens-source relative proper motion), measured from ecliptic East toward ecliptic
North in the standard `atan2(y, x)` convention.

Since pi_E is parallel to the relative proper motion, and alpha points in the
opposite direction:

```
alpha = phi_pi + pi + (pi/2 - phi_pi - (pi/2 - theta)) ...
```

More precisely, using `theta = atan2(mubet_r, mulam_r)`:

```
phi_pi = pi/2 - theta
alpha  = theta + pi
```

Therefore:

```
cos(alpha) = -sin(phi_pi)
sin(alpha) = -cos(phi_pi)
```

### Perpendicular axis conventions differ

**In `compute_tushifts`**, the `(tau, u)` decomposition is:

```
tau_raw =  cos(phi_pi) * N + sin(phi_pi) * E  =  v_E * sin(phi_pi) + v_N * cos(phi_pi)
u_raw   = -sin(phi_pi) * N + cos(phi_pi) * E  =  v_E * cos(phi_pi) - v_N * sin(phi_pi)
```

The perpendicular axis `u_raw` points 90 degrees **clockwise** from the
pi_E direction.

**In the alpha rotation**, the decomposition is:

```
xs0 = tau * cos(alpha) - u * sin(alpha)    [ecliptic E]
ys0 = tau * sin(alpha) + u * cos(alpha)    [ecliptic N]
```

Here `u` is 90 degrees **counter-clockwise** from the tau direction.

Since pi_E and source motion are anti-parallel, "clockwise from pi_E" and
"counter-clockwise from source motion" point in **opposite** directions.
Therefore `u_raw = -u_alpha`, and the `u` component must be negated between the
two rotations.

### Verification: round-trip with the NEW sign

With the corrected formula:

```
tshift = -piE * tau_raw
ushift = +piE * u_raw     (note: opposite sign from tshift)
```

The event-frame parallax contribution is:

```
delta_E = tshift * cos(alpha) - ushift * sin(alpha)
delta_N = tshift * sin(alpha) + ushift * cos(alpha)
```

Substituting `cos(alpha) = -sin(phi_pi)`, `sin(alpha) = -cos(phi_pi)`,
expanding `tau_raw` and `u_raw`, and using `sin^2 + cos^2 = 1`:

```
delta_E = piE * Eshift
delta_N = piE * Nshift
```

This is a clean result with no residual dependence on `phi_pi`.  The old code
(with `-piE` on both) produces:

```
delta_E = -piE * [Eshift * cos(2*phi_pi) - Nshift * sin(2*phi_pi)]
delta_N = -piE * [Nshift * cos(2*phi_pi) + Eshift * sin(2*phi_pi)]
```

which spuriously rotates the parallax vector by `2*phi_pi`, corrupting
astrometry while leaving `|delta|` (and thus photometric magnification)
unaffected.

### Why photometry was never affected

Photometric magnification depends on `|u|^2 = tau^2 + u^2`.  Negating `u` does
not change `|u|`.  The sign of `u0` similarly cancels from `u0^2`.  Both bugs
are invisible to the lightcurve; only the astrometric centroid direction is
affected.

---

## Why the `u0` sign change is correct

GULLS defines `u0lens1` in the **lens-source (LS)** convention: the
perpendicular impact parameter measures the source position relative to the
lens.  BAGLE's `PSPL_PhotAstrom_Par_Param4_geoproj` model expects the
**source-lens (SL)** convention, which is negated.

Like `ushift`, `|u0|` cancels from photometric magnification, so the sign error
was invisible in lightcurves but placed the astrometric centroid on the wrong
side of the lens.

---

## How the BAGLE comparison works

### Model used

`PSPL_PhotAstrom_Par_Param4_geoproj` — a point-source point-lens model with
parallax, astrometry, and geocentric-projected parameterization.

### Parameter mapping: GULLS → BAGLE

| BAGLE parameter | Source | Conversion |
|---|---|---|
| `t0` (MJD) | `.out` `t0lens1` | `sim_zero_jd + t0lens1 - 2400000.5` |
| `u0` | `.out` `u0lens1` | **Negated** (LS → SL convention) |
| `tE` (days) | `.out` `tE_ref` | Direct (reference-frame timescale) |
| `thetaE` (mas) | `.out` `thetaE` | Direct |
| `piS` (mas) | `.out` `Source_Dist` | `1 / Source_Dist_kpc` |
| `piE_E` | `.out` `piE`, `murel_ref_alpha`, `murel_ref_delta` | `piE * murel_ref_alpha / |murel_ref|` |
| `piE_N` | (same) | `piE * murel_ref_delta / |murel_ref|` |
| `xS0` (arcsec) | `.lc` `true_RA_deg`, `true_Dec_deg` at reference epoch | `dRA * cos(dec) * 3600`, `dDec * 3600` |
| `muS_E` (mas/yr) | `.out` `mu_source_helio_alpha` | Direct (heliocentric, equatorial) |
| `muS_N` (mas/yr) | `.out` `mu_source_helio_delta` | Direct |
| `t0par` (MJD) | `.out` `tref` | Same as `t0` conversion |
| `raL`, `decL` (deg) | `.out` `ra_deg`, `dec_deg` | Direct |
| `mag_base` | `.lc` header `Obssrcmag` + `.out` `Obs_0_fs` | `source_mag + 2.5*log10(fs)` |
| `b_sff` | `.out` `Obs_0_fs` | Direct |

### Critical detail: piE direction projection

GULLS publishes `piE` (scalar), `piEE`, and `piEN` in **ecliptic** coordinates.
BAGLE expects `piE_E` and `piE_N` in **equatorial** (RA*cos(dec), Dec)
coordinates.  Rather than rotating the ecliptic `(piEE, piEN)` pair through
`muecl2ad`, the comparison code uses the already-published equatorial
reference-frame relative proper motion `(murel_ref_alpha, murel_ref_delta)` as
the direction vector, since pi_E is parallel to mu_rel:

```python
piE_E = piE_amp * murel_ref_alpha / |murel_ref|
piE_N = piE_amp * murel_ref_delta / |murel_ref|
```

This avoids a redundant coordinate rotation and uses the same published
quantities that downstream fitters will consume.

### What is compared

For each of 50 simulated 1-source-1-lens events:

1. **Photometry**: GULLS noiseless `true_relative_flux` vs
   `10^(-0.4*(BAGLE_mag - mag_base))`.
2. **Astrometry**: GULLS noiseless `(true_RA_deg, true_Dec_deg)` converted to
   tangent-plane `(dRA*cos(dec), dDec)` arcsec vs BAGLE
   `get_astrometry(t_mjd)`, with a constant translation matched at the
   reference epoch.

### Results (50 events, post-fix)

| Metric | Min | Median | Max |
|---|---|---|---|
| Astrometric RMS (mas) | 8.77e-05 | 1.19e-04 | 2.39e-04 |
| Photometric RMS (relative flux) | 6.71e-07 | 2.82e-05 | 1.83e-02 |

**Astrometric agreement: sub-0.001 mas for all 50 events.**  The residual
~0.1 microarcsecond floor is consistent with differences between GULLS's
Keplerian ephemeris and BAGLE's ERFA-based Earth position.

**Photometric outliers** (4 events with phot RMS > 0.001) all have small `|u0|`
(0.01--0.09) where finite-source effects dominate.  The PSPL forward model
cannot match finite-source magnification peaks; this is expected and is not an
astrometry bug.

---

## What the existing documentation got wrong

The previous `astrometry.rst` was written during development before the
coordinate pipeline was finalized.  The following claims were incorrect:

1. **"x,y orientation relative to celestial E,N is UNCERTAIN"** — the event
   frame IS ecliptic E,N by construction.  Alpha is computed from the
   reference-frame relative proper motion in ecliptic coordinates (line 165 of
   `omLightcurveGenerator.cpp`), and the `(tau, u) → (xs0, ys0)` rotation is
   designed to produce ecliptic `(E, N)`.

2. **"Does NOT account for lens proper motion"** — lens proper motion IS
   included (`pm_lam_mas`, `pm_beta_mas` in `photometry.cpp` line 230).

3. **Column names** — the documentation listed columns like `centroid_src_x_mas`,
   `RA_centroid_deg`, `RA_centroid_lpllx_deg` that do not exist in the actual
   output.  Actual column names are `blended_sources_only_x_thetaE`,
   `true_RA_deg`, etc.

4. **`#Astrometry_BAGLE` contract** — documented as `model_frame=lens_relative`
   with `blendless_columns=RA_centroid_src_only_deg,...`, but the actual header
   writes `model_frame=absolute` and `blendless_columns=none`.

5. **Units** — documented centroid columns as being in mas; actual event-frame
   columns are in theta_E.  Conversion to mas happens in `photometry.cpp`.

The rewritten `astrometry.rst` corrects all of these.

---

## Remaining known limitations

- **Ambient star blending**: ambient stars in the PSF aperture are not yet
  blended into the centroid (currently placeholder).
- **Commented-out overload**: the commented `compute_tushifts(vector<...>*)`
  overload at parallax.cpp line 187--188 still uses the old `-piE` convention.
  If ever re-enabled, it must be updated.
- **PSPL-only BAGLE comparison**: the forward sanity check uses PSPL, so
  finite-source events show photometric residuals.  An FSPL comparison would
  strengthen validation for high-magnification events.
- **Ephemeris differences**: GULLS uses Keplerian orbital elements for Earth's
  position; BAGLE uses ERFA (IAU SOFA derivative).  This produces ~0.1 µas
  astrometric floor that cannot be eliminated without ephemeris unification.

## Test plan

- [x] `smoke_test/run_bagle_forward_sanity.py` — 50 events, all astrometric
  RMS < 0.001 mas
- [x] Mathematical proof that the sign change produces a clean round-trip
  (no residual phi_pi dependence)
- [x] Photometric RMS unaffected (same-magnitude residuals as pre-fix)
- [ ] Re-run `smoke_test/run_astrometry_sanity.py` with full validation suite
- [ ] Re-run `smoke_test/run_bagle_fit_sanity.py` joint fit check
