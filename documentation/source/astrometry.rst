Astrometry outputs and parameters
=================================

This page summarizes the astrometry options added to gulls, the new output columns, and the noise recipe used to generate observed astrometric positions.

Parameters
----------

- ``ASTROMETRY_ON`` (default: ``0``)
  - Enable astrometric computation and outputs when set to 1.
  - When 0, all sky-frame astrometry outputs are disabled and written as 0.0.

- ``ASTROMETRIC_SYS_FLOOR`` (units: mas, default: ``0.1``)
  - Per-axis systematic floor for astrometric uncertainty.
  - Combined in quadrature with the photon-limited term when producing per-epoch errors.

New output columns (NE sky frame)
---------------------------------

All angles below are in milliarcseconds (mas) unless stated otherwise.

- ``true_N_centroid_mas``: True (noise-free) North offset of the centroid in the sky frame.
- ``true_E_centroid_mas``: True (noise-free) East offset of the centroid in the sky frame.
- ``obs_N_centroid_mas``: Observed North centroid with noise applied.
- ``obs_E_centroid_mas``: Observed East centroid with noise applied.
- ``obs_N_centroid_err_mas``: 1-sigma per-axis uncertainty (North) for the observed centroid.
- ``obs_E_centroid_err_mas``: 1-sigma per-axis uncertainty (East) for the observed centroid.

Absolute RA/Dec centroid columns
--------------------------------

When astrometry is enabled, additional columns report the absolute sky position of the flux centroid in Equatorial coordinates per epoch, including lens proper motion drift:

- ``true_centroid_ra_deg``, ``true_centroid_dec_deg``: True (noise-free) absolute centroid coordinates in degrees.
- ``measured_centroid_ra_deg``, ``measured_centroid_dec_deg``: Observed centroid coordinates (noise applied) in degrees.
- ``measured_centroid_ra_error_deg``, ``measured_centroid_dec_error_deg``: 1-sigma uncertainties per axis in degrees.

Notes:
- The absolute centroid is computed as the catalog position at t0 plus the sum of two small-angle terms: (i) microlensing NE offsets and (ii) linear lens proper motion drift converted from Galactic to Equatorial components. For small angles, RA offset is scaled by cos(Dec).
- RA values are wrapped into [0, 360) degrees.

Notes on frames and units:
- ``x``/``y`` centroid columns in the lightcurve remain in the lens-frame coordinates (x1/x2) and are expressed in Einstein radii; the NE columns above are in the sky frame (mas). These are different coordinate systems by design.

Noise model and zeroing rules
-----------------------------

Per-epoch astrometric uncertainties and observed values are generated as follows:

1. Compute a photometric signal-to-noise ratio per epoch from the simulated photometry:
   - ``SNR = |Aobs| / max(Aerr, 1e-12)``
2. PSF width is taken from the instrument model as ``FWHM`` in arcsec; we convert to mas via ``FWHM_mas = 1000 * FWHM``.
3. Photon-limited per-axis precision (1D) is approximated as:
   - ``sigma_photon = FWHM_mas / SNR``
4. Total per-axis uncertainty combines the photon term with a systematic floor:
   - ``sigma_axis = sqrt(sigma_photon^2 + ASTROMETRIC_SYS_FLOOR^2)``
5. Observed NE centroids are produced by adding independent Gaussian noise to the true NE centroids:
   - ``obs_N = true_N + Normal(0, sigma_axis)``
   - ``obs_E = true_E + Normal(0, sigma_axis)``
   - ``obs_N_centroid_err_mas = obs_E_centroid_err_mas = sigma_axis``

Zeroing behavior (no sky orientation or disabled):
- If ``ASTROMETRY_ON = 0``, all NE outputs are set to 0.0.
- If no reliable sky orientation is available (e.g., neither parallax direction nor relative proper motion direction), the NE outputs are set to 0.0 for that event.
- In ideal photometry mode, the observed NE centroids equal the true NE centroids and their errors are set to the systematic floor (no random noise added).

Lens proper motion and absolute coordinates
-------------------------------------------

- NE outputs (true/obs) are centroid offsets relative to the lens-frame origin, rotated into the local sky North/East axes. Absolute RA/Dec columns incorporate the drift of the lens system through the sky due to its proper motion.
- Implementation details:
   - Convert Galactic proper motions (``μ_l, μ_b``) to Equatorial (``μ_RA*``, ``μ_Dec``) via ``coords.mulb2ad``.
   - Compute ``Δt = (epoch - t0) / 365.25`` years and form NE drifts: ``ΔE_PM = μ_RA* × Δt`` and ``ΔN_PM = μ_Dec × Δt`` (mas).
   - Total NE offsets = (microlensing NE) + (PM NE). Convert to degrees: ``ΔRA = ΔE / (cos(Dec) × 3600000)`` and ``ΔDec = ΔN / 3600000``. Add to base ``(RA,Dec)`` at ``t0``.
   - Uncertainties per axis are converted from mas to degrees using the same small-angle relations.
   - RA values are wrapped into [0, 360) degrees.
