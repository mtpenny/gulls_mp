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
  - in `--fit-true-astrometry` mode, uses blendless source-only astrometry (`RA_centroid_src_only_deg/Dec_centroid_src_only_deg`) for fitting and strict RMS checks;
  - writes fit summary + diagnostic plot and records which noiseless astrometry columns were used;
  - compares fitted PM and parallax vectors (amplitude + direction) against `.out` in the declared model frame convention;
  - falls back to scipy least-squares BAGLE fit if PyMultiNest runtime is unavailable.
