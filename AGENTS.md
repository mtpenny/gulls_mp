Things we know from Synthpop:

Lens proper motion (mas yr-1)
Source proper motion (mas yr-1)
Lens magnitude
Source magnitude
Mass of the lens
Ra & dec? (this might be re-drawn in gulls because the source and lens have different coordinates and gulls do not actually simulate the galaxy, it forces things to align and calculates the probability of this happening)
Lens distance
Source distance

Things we know from Gulls:

s (calculated based on drawn orbital parameters and defined at time t_ref)
q (calculated from a drawn planet mass and synthpop star mass)
alpha (drawn randomly at time t_ref)
t0 (randomly drawn, usually forced to be in season, relative to the largest mass as the origin for the impact parameter)
Relative proper motion (lens-source from synthpop)
piE (magnitude based on distances and masses, direction from relative proper motion)
thetaE (calculated from synthpop parameters)
Multilens lens-plane positions (orbital dynamics are calculated hierarchically from most massive object to least)
l,b pointing in galactic coordinates is canonically true as presented in the out file
"Event frame" is ecliptic (internal simulator frame, not directly published sky EN)

---

Astrometry conventions (general executable only)
------------------------------------------------

- Internal `Event->x/y` centroid arrays are in an event frame, not directly sky EN.
  - Trajectory orientation in that frame is controlled by random `alpha` (deg).
  - Public astrometry outputs rotate event-frame centroids to observer-centric ecliptic EN:
    - `gamma_obs = atan2(murel_ref_beta, murel_ref_lambda) - alpha_deg`
    - `E = cos(gamma_obs)*x_evt - sin(gamma_obs)*y_evt`
    - `N = sin(gamma_obs)*x_evt + cos(gamma_obs)*y_evt`
- Public centroid columns in `.lc` (`x_centroid_mas`, `y_centroid_mas`, true/source/lens blended diagnostics) are written in ecliptic EN mas after this rotation.
- `RA_*` / `Dec_*` columns are produced from ecliptic EN via local tangent-plane transform using `coords::muecl2ad`.
  - Use `dRAcosDec` (not raw `dRA`) internally.
  - BAGLE conversion (documented in `.lc` header):
    - `x_E_arcsec = dRAcosDec_mas / 1000`
    - `y_N_arcsec = dDec_mas / 1000`
- `#Astrometry_BAGLE` is treated as a required contract for BAGLE checks:
  - `model_frame=lens_relative` means BAGLE model astrometry must use `get_astrometry - get_lens_astrometry`.
  - `blendless_columns=RA_centroid_src_only_deg,Dec_centroid_src_only_deg` are the preferred noiseless source-only columns for BAGLE fitting/validation.
  - `lens_columns=RA_lens_primary_deg,Dec_lens_primary_deg` provide primary-lens sky track in the same published RA/Dec frame.
- `lens_parallax_x_mas`, `lens_parallax_y_mas` are ecliptic EN lens-parallax terms:
  - `lens_parallax_E = -Eshift / D_L_kpc`
  - `lens_parallax_N = -Nshift / D_L_kpc`
- Time columns:
  - `.lc` `BJD` is JD-like; BAGLE expects MJD: `MJD = BJD - 2400000.5`.
  - If `BJD` is quantized/malformed, use `SIMULATION_ZERO_TIME + Simulation_time`, then JD->MJD.

Tests that document/enforce conventions
--------------------------------------

- `smoke_test/astrometry_sanity.py`:
  - validates declared frame transforms from header and RA/Dec reconstruction;
  - validates lpllx-vs-non-lpllx offsets against documented `lens_parallax_x/y`;
  - verifies non-constant observer orbit columns;
  - checks PM behavior near `tref` and long-baseline behavior (with parallax correction);
  - blend-aware rule: if blended-centroid long-baseline PM fails but source-only passes and `Obs_0_fs < 0.9`, warn/pass as likely blend-driven.
- `smoke_test/bagle_fit_sanity.py`:
  - runs BAGLE joint photometric+astrometric PSPL+parallax fit on selected single-source event (`ObsGroup_0_chi2 < 20` by default);
  - explicitly sets BAGLE `obsLocation` (default request `jwst`, automatic fallback to `earth` if initialization fails);
  - in `--fit-true-astrometry` mode, uses blendless source-only astrometry (`RA_centroid_src_only_deg/Dec_centroid_src_only_deg`) for fitting and strict RMS checks;
  - writes fit summary + diagnostic plot and records which noiseless astrometry columns were used;
  - compares GULLS primary-lens track (`RA_lens_primary_deg/Dec_lens_primary_deg`) against `best_model.get_lens_astrometry(t)` and records both raw RMS and XY-offset-removed RMS;
  - compares fitted PM and parallax vectors (amplitude + direction) against `.out` in the declared model frame convention;
  - falls back to scipy least-squares BAGLE fit if PyMultiNest runtime is unavailable.

Smoke observer notes
--------------------

- `smoke_test/assets/observatories/smoke.observatory` uses `SPACE=1` and `ORBIT=0` (Earth orbit path from `setupOrbit`).
- Weather is still loaded globally, but `smoke_test/assets/weather/smoke.weather` is all `1`s, so cadence is not weather-thinned.

Astrometric Deflection Sanity Check (Source-Dominated)
------------------------------------------------------
To verify astrometry is working under the hood, consider a source-dominated centroid track for a source at ~5 kpc ($\pi_S = 0.2$ mas). Instead of assuming a standard May-Sept season, we can mathematically verify the exact solar system locations appended to the `general-single` `.lc` file against the observed deflection:
- **The Geometry:** The Galactic bulge is at ecliptic longitude $\lambda \approx 266^\circ$. The parallactic shift is $\Delta\mathbf{s} = - \mathbf{r}_\oplus / D_S$. The observed displacement is anti-parallel to Earth's positional vector.
- **Simulation Timeline:** The `smoke_general` event simulates 200 days from $t=0$ (BJD 2458849, exactly Jan 1st, 2020) to mid-July. The event $t_0$ peaks around $t=75$ (mid-March).
- **Tracing the Deflection:** 
  - **Start (Jan 1, $t=0$):** Earth is at $\lambda_\oplus \approx 100^\circ$ ($x_\oplus \approx -0.16 \text{ AU}, y_\oplus \approx 0.97 \text{ AU}$). Projected against the $266^\circ$ line of sight, the transverse Eastward deflection ($\Delta E \propto x_\oplus \sin\lambda - y_\oplus \cos\lambda$) is a mild $\sim +0.05$ mas East.
  - **Peak Deflection (Late March, $t=80$):** Earth swings to $\lambda_\oplus \approx 170^\circ$ ($x \approx -0.97 \text{ AU}, y \approx 0.18 \text{ AU}$), maximizing its perpendicular distance to the bulge. The source shift arcs maximally to the **East** ($\sim +0.2$ mas).
  - **End (Mid-July, $t=200$):** Earth orbits around past opposition to $\lambda_\oplus \approx 296^\circ$ ($x \approx 0.45 \text{ AU}, y \approx -0.91 \text{ AU}$). The apparent source shift swings across 0 and deflects to the **West** ($\sim -0.1$ mas).
- **The Result:** The simulated blue source-dominated centroid traces an arc that smoothly swings $\sim +0.2$ mas *Eastward* around the March event peak before curving strongly back toward the *West* by July. The data physically match the exact ecliptic $x/y/z$ coordinates appended in the `parallax_shift_*` output columns, confirming the transformations and alignments are fully robust.
