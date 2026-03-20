Astrometry outputs and validation
=================================

This page documents the astrometry feature in the ``general`` gulls executable:
the coordinate pipeline, output columns, noise model, parameter conventions,
and validation against the BAGLE microlensing fitting framework.

.. contents:: Contents
   :local:
   :depth: 2


Parameters
----------

- ``ASTROMETRY_ON`` (default: ``0``)
  When 1, per-epoch astrometric centroids are computed, converted to sky
  coordinates, and written to the ``.lc`` file.  When 0, all astrometry
  columns are written as 0.0.

  **Warning**: enabling astrometry increases VBMicrolensing evaluation cost
  because centroids must be computed for each source image.

- ``ASTROMETRIC_SYS_FLOOR`` (units: mas, default: ``0.1``)
  Per-axis systematic floor added in quadrature with the photon-limited
  astrometric uncertainty.


Coordinate pipeline
-------------------

The astrometric centroid passes through five stages.  Each stage is stored
in the Event structure for diagnostic output.

Stage 1: VBM internal frame (theta_E)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

VBMicrolensing computes the flux-weighted centroid of the lensed source
images in its own internal coordinate system:

- *Single lens*: x1 along source-lens axis, x2 perpendicular.
- *Binary lens*: x1 along binary axis, x2 perpendicular, origin at center
  of mass.

Output columns: ``vbm_astrox1_source{i}_thE``, ``vbm_astrox2_source{i}_thE``.

Stage 2: event frame (theta_E)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``omLightcurveGenerator.cpp`` transforms VBM centroids to the event frame
by applying the same coordinate rotation used for source/lens positions.

The event frame is defined so that its axes align with **ecliptic East (E)
and ecliptic North (N)**, in theta_E units, with the origin at the primary
lens center of mass.  This alignment is achieved by computing the trajectory
angle alpha from the reference-frame ecliptic relative proper motion::

    alpha = atan2(-mubet_r, -mulam_r)

where ``mulam_r`` and ``mubet_r`` are the ecliptic (lambda, beta) components
of the reference-frame relative proper motion unit vector.  This is the
direction of source motion relative to the lens, measured from ecliptic East
toward ecliptic North (standard ``atan2(y, x)`` convention).

The source's (tau, u) coordinates are rotated by alpha to produce event-frame
(x, y)::

    xs0 = tau * cos(alpha) - u * sin(alpha)    [ecliptic E, theta_E]
    ys0 = tau * sin(alpha) + u * cos(alpha)    [ecliptic N, theta_E]

Multiple sources are flux-weighted using magnified fluxes::

    centroid = sum(flux_i * mu_i * centroid_i) / sum(flux_i * mu_i)

Output columns:
``blended_sources_only_x_thetaE``, ``blended_sources_only_y_thetaE``
(source-only centroid, no lens light).

Stage 3: luminous lens blending (theta_E)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In ``photometry.cpp``, the source-only centroid is flux-weighted with up to
two luminous lens positions::

    centroid = (centroid_src * f_src_total
                + pos_lens1 * fl1
                + pos_lens2 * fl2) / (f_src_total + fl1 + fl2)

where ``fl1``, ``fl2`` are computed from lens magnitudes relative to the
primary source zero point.

Output columns:
``blended_sources_lenses_x_thetaE``, ``blended_sources_lenses_y_thetaE``.

Stage 4: ecliptic tangent plane (degrees)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Still in ``photometry.cpp``, the blended centroid is converted to absolute
ecliptic coordinates:

1. Multiply by ``thE_mas`` to get mas offset from the lens center of mass.
2. Add the **lens reference-frame proper motion**: ``pm_lam_mas``,
   ``pm_beta_mas`` (ecliptic lambda and beta, scaled by dt from tref).
3. Add the **lens parallax**: ``-Eshift / D_L_kpc``, ``-Nshift / D_L_kpc``
   (the observer's transverse displacement, divided by lens distance, with
   a sign flip because the apparent lens motion is opposite to the observer's
   displacement).
4. Convert mas offset to a radians offset via the ecliptic tangent plane::

       dlambda_rad = dE_mas * mas_to_rad / cos(beta0)
       dbeta_rad   = dN_mas * mas_to_rad

5. Add to the reference ecliptic coordinates ``(lambda0, beta0)``.

Output columns:
``ecliptic_lambda_noiseless_deg``, ``ecliptic_beta_noiseless_deg``.

Stage 5: equatorial ICRS (degrees)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The absolute ecliptic coordinates are converted to equatorial (RA, Dec) by
``coords::ecl2ad``, a standard obliquity-based spherical transform.

Output columns (noiseless): ``true_RA_deg``, ``true_Dec_deg``.

Output columns (with noise): ``measured_RA_deg``, ``measured_Dec_deg``.


Noise model
-----------

Per-epoch astrometric uncertainties follow Gould & Yee (2014), computed
entirely in mas:

1. Fractional photometric error::

       sigma_phot = |Aerr / Aobs|

2. Photon-limited astrometric uncertainty::

       sigma_ast_psf_mas = FWHM_arcsec * 1000 * sigma_phot / sqrt(ln 256)

3. Total uncertainty with systematic floor::

       sigma_ast_mas = sqrt(sigma_ast_psf_mas^2 + ASTROMETRIC_SYS_FLOOR^2)

4. Independent Gaussian noise is added to RA and Dec::

       measured_RA_deg  = true_RA_deg  + sigma_ast_mas * N(0,1) * mas_to_deg / cos(dec0)
       measured_Dec_deg = true_Dec_deg + sigma_ast_mas * N(0,1) * mas_to_deg

Output columns: ``sigma_astrometric_mas``, ``measured_RA_error_deg``,
``measured_Dec_error_deg``.


Output column reference
-----------------------

All columns are per-epoch rows in the ``.lc`` file.

.. list-table:: Lightcurve astrometry columns
   :header-rows: 1
   :widths: 40 15 45

   * - Column
     - Units
     - Description
   * - ``blended_sources_only_x_thetaE``
     - theta_E
     - Flux-weighted source centroid, ecliptic E, no lens light
   * - ``blended_sources_only_y_thetaE``
     - theta_E
     - Flux-weighted source centroid, ecliptic N, no lens light
   * - ``blended_sources_lenses_x_thetaE``
     - theta_E
     - Blended centroid with luminous lenses, ecliptic E
   * - ``blended_sources_lenses_y_thetaE``
     - theta_E
     - Blended centroid with luminous lenses, ecliptic N
   * - ``ecliptic_lambda_noiseless_deg``
     - degrees
     - Noiseless blended centroid, ecliptic longitude
   * - ``ecliptic_beta_noiseless_deg``
     - degrees
     - Noiseless blended centroid, ecliptic latitude
   * - ``true_RA_deg``
     - degrees
     - Noiseless blended centroid, equatorial RA (ICRS)
   * - ``true_Dec_deg``
     - degrees
     - Noiseless blended centroid, equatorial Dec (ICRS)
   * - ``measured_RA_deg``
     - degrees
     - Observed centroid with noise, equatorial RA
   * - ``measured_Dec_deg``
     - degrees
     - Observed centroid with noise, equatorial Dec
   * - ``sigma_astrometric_mas``
     - mas
     - 1-sigma astrometric uncertainty per axis
   * - ``measured_RA_error_deg``
     - degrees
     - 1-sigma RA uncertainty (= sigma_ast_mas / cos(dec) in deg)
   * - ``measured_Dec_error_deg``
     - degrees
     - 1-sigma Dec uncertainty (= sigma_ast_mas in deg)
   * - ``fractional_total_source_flux``
     - dimensionless
     - Magnified total source flux / unmagnified primary source flux
   * - ``parallax_shift_t``
     - theta_E
     - Parallax shift in tau (along relative motion)
   * - ``parallax_shift_u``
     - theta_E
     - Parallax shift in u (perpendicular to relative motion)
   * - ``parallax_shift_x``
     - AU
     - Observer ecliptic-x position (3D, from orbital elements)
   * - ``parallax_shift_y``
     - AU
     - Observer ecliptic-y position (3D)
   * - ``parallax_shift_z``
     - AU
     - Observer ecliptic-z position (3D)
   * - ``vbm_astrox1_source{i}_thE``
     - theta_E
     - Raw VBM x1 centroid output for source i
   * - ``vbm_astrox2_source{i}_thE``
     - theta_E
     - Raw VBM x2 centroid output for source i
   * - ``source{i}_x_thE``
     - theta_E
     - Source i position, ecliptic E (event frame)
   * - ``source{i}_y_thE``
     - theta_E
     - Source i position, ecliptic N (event frame)
   * - ``source{i}_mu``
     - dimensionless
     - Magnification of source i
   * - ``lens{i}_x_thE``
     - theta_E
     - Lens i position, ecliptic E (event frame)
   * - ``lens{i}_y_thE``
     - theta_E
     - Lens i position, ecliptic N (event frame)


Lightcurve header metadata
--------------------------

The ``.lc`` header contains astrometry metadata needed for downstream analysis:

``#Astrometry_Frame``
    Reference position and scale::

        #Astrometry_Frame: RA_rad=X Dec_rad=Y RA_deg=X Dec_deg=Y
            lambda0_deg=X beta0_deg=Y thE_mas=Z tref=T
            frame=ecliptic_EN_barycenter

    The reference position is the catalog lens position.  ``thE_mas`` is the
    angular Einstein radius for unit conversion.  ``frame=ecliptic_EN_barycenter``
    declares that event-frame x/y axes are ecliptic East/North.

``#Astrometry_PM``
    Published proper motions for source, lens, and relative (lens minus source)
    in heliocentric and reference-frame variants, in equatorial and ecliptic
    components.  Units are mas/yr.  Equatorial components use the
    ``mu_alpha*cos(delta)`` convention (not raw ``mu_alpha``).

``#Astrometry_Contract``
    Machine-readable contract for downstream consumers::

        #Astrometry_Contract: version=v2 model_frame=absolute
            event_true_unit=thetaE noise_frame=ecliptic_tangent

    - ``model_frame=absolute``: sky astrometry columns are absolute RA/Dec,
      not lens-relative offsets.
    - ``event_true_unit=thetaE``: event-frame centroid columns are in theta_E.
    - ``noise_frame=ecliptic_tangent``: noise is applied in the ecliptic
      tangent plane.

``#Astrometry_Columns``
    Maps logical column roles to actual column names.

``#Astrometry_BAGLE``
    Contract for BAGLE-specific validation::

        #Astrometry_BAGLE: model_frame=absolute blendless_columns=none
            lens_columns=lens0_x_thE,lens0_y_thE lens_frame=event_xy_thetaE

    - ``model_frame=absolute``: BAGLE model astrometry should be compared
      directly to ``true_RA_deg``/``true_Dec_deg`` (not lens-relative).
    - ``blendless_columns=none``: no separate source-only RA/Dec columns are
      currently published (source-only centroids are available in theta_E only).
    - ``lens_columns=lens0_x_thE,lens0_y_thE``: primary lens track is
      available in event-frame theta_E.


The parallax sign fix
---------------------

A sign error in ``parallax.cpp::compute_tushifts`` (line 152) was corrected.
This section provides the mathematical proof that the fix is correct.

The problem
~~~~~~~~~~~

The parallax perturbation passes through two 2D rotations:

1. ``compute_tushifts``: ecliptic ``(Nshift, Eshift)`` in AU, rotated by
   ``phi_pi`` into ``(tshift, ushift)`` in theta_E.
2. ``omLightcurveGenerator``: ``(tau + tshift, u0 + ushift)`` rotated by
   ``alpha`` into event-frame ``(xs0, ys0)`` in theta_E.

The ``u``-perpendicular axis has **opposite orientation** in these two rotations:

- In ``compute_tushifts``, ``u_raw`` points 90 degrees **clockwise** from
  the pi_E direction (decomposition: ``u_raw = v_E * cos(phi_pi) - v_N * sin(phi_pi)``).
- In the alpha rotation, ``u`` points 90 degrees **counter-clockwise** from
  the tau direction (the standard 2D rotation convention).

Since pi_E and source motion (tau direction) are anti-parallel, "clockwise from
pi_E" and "counter-clockwise from source motion" point in opposite directions.
This means ``u_raw = -u_alpha``.

The fix
~~~~~~~

The old code used the same sign (``-piE``) for both tau and u::

    tshift[i] = -piE * ( Nshift[i]*cs + Eshift[i]*sn);   // correct
    ushift[i] = -piE * (-Nshift[i]*sn + Eshift[i]*cs);   // WRONG

The fix uses opposite signs::

    tshift[i] = -piE * ( Nshift[i]*cs + Eshift[i]*sn);   // unchanged
    ushift[i] =  piE * (-Nshift[i]*sn + Eshift[i]*cs);   // sign flipped

Round-trip proof
~~~~~~~~~~~~~~~~

With the corrected formula and the identities ``cos(alpha) = -sin(phi_pi)``,
``sin(alpha) = -cos(phi_pi)`` (derived from the angle definitions), the
event-frame parallax contribution is::

    delta_E = tshift * cos(alpha) - ushift * sin(alpha)
            = -piE * tau_raw * (-sin phi_pi) - piE * u_raw * (-(-cos phi_pi))

Expanding ``tau_raw`` and ``u_raw`` in terms of ``(Eshift, Nshift)`` and
applying ``sin^2 + cos^2 = 1``::

    delta_E = piE * Eshift
    delta_N = piE * Nshift

The result has **no residual dependence on phi_pi**, confirming the rotation
round-trip is clean.

With the old sign (``-piE`` on both), the result contains ``cos(2*phi_pi)``
and ``sin(2*phi_pi)`` terms that spuriously rotate the parallax vector —
corrupting astrometry while preserving ``|delta|`` (and hence photometry).

Why photometry was unaffected
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Photometric magnification depends on ``|u|^2 = tau^2 + u^2``, which is
invariant under sign flips of ``u``.  The sign of ``u0`` similarly cancels
from ``u0^2``.  Both the ``ushift`` and ``u0`` sign errors are invisible
to the lightcurve; only the astrometric centroid direction is affected.


Sign conventions
----------------

GULLS uses the **lens-source (LS)** convention for the relative proper motion
and impact parameter:

- ``murel = mu_lens - mu_source``
- ``u0lens1``: perpendicular source-lens separation at ``t0``, in LS convention
- ``piE``: parallel to ``murel``; ``piEN`` = ecliptic beta component,
  ``piEE`` = ecliptic lambda component
- ``phi_pi = atan2(piEE, piEN)``: angle of pi_E from ecliptic N toward E

BAGLE's ``PSPL_PhotAstrom_Par_Param4_geoproj`` uses the **source-lens (SL)**
convention for ``u0``.  When comparing GULLS to BAGLE:

- ``u0_bagle = -u0lens1``
- ``piE`` direction is the same (both codes use LS for pi_E)
- ``muS`` (heliocentric equatorial) passes through directly
- ``piE_E/piE_N`` must be projected from ecliptic to equatorial:
  the comparison code uses the published equatorial relative proper motion
  direction ``(murel_ref_alpha, murel_ref_delta)`` as the pi_E direction
  vector, since pi_E is parallel to mu_rel.


Validation against BAGLE
-------------------------

Forward model comparison
~~~~~~~~~~~~~~~~~~~~~~~~

``smoke_test/bagle_forward_sanity.py`` constructs a BAGLE
``PSPL_PhotAstrom_Par_Param4_geoproj`` model from the GULLS ``.out``
parameters for each 1-source-1-lens event, evaluates it at all GULLS
epochs, and compares:

- **Photometry**: GULLS noiseless ``true_relative_flux`` vs BAGLE model
  flux, converted to the same relative-flux convention.
- **Astrometry**: GULLS noiseless ``(true_RA_deg, true_Dec_deg)`` converted
  to tangent-plane arcsec vs BAGLE ``get_astrometry(t_mjd)``, with a
  constant translation matched at the reference epoch.

Results (50 simulated events):

.. list-table::
   :header-rows: 1
   :widths: 30 20 20 20

   * - Metric
     - Min
     - Median
     - Max
   * - Astrometric RMS (mas)
     - 8.77e-05
     - 1.19e-04
     - 2.39e-04
   * - Photometric RMS (rel flux)
     - 6.71e-07
     - 2.82e-05
     - 1.83e-02

All 50 events achieve astrometric RMS below 0.001 mas.  The ~0.1
microarcsecond residual floor is consistent with differences between GULLS's
Keplerian Earth ephemeris and BAGLE's ERFA-based ephemeris.

Photometric outliers (4 events with RMS > 0.001) all have small ``|u0|``
(0.01--0.09) where finite-source effects dominate.  The PSPL forward model
is expected to fail for these events; this is not an astrometry error.

The parameter mapping (GULLS ``.out`` columns to BAGLE constructor arguments)
is recorded in each event's JSON summary file under ``truth_mapping``, along
with both convention-explicit names and numerical values.

Astrometry sanity checks
~~~~~~~~~~~~~~~~~~~~~~~~~

``smoke_test/astrometry_sanity.py`` performs structural validation:

- Required column presence and header contract parsing.
- Observer orbit non-constancy (``parallax_shift_x/y/z``).
- Reference-frame proper motion recovery from finite-differenced source/lens
  positions.
- Long-baseline heliocentric PM comparison (far-field epochs, with parallax
  correction) against published ``murel_helio_alpha/delta``.
- Astropy cross-check of proper motion coordinate rotations.

BAGLE joint fit
~~~~~~~~~~~~~~~

``smoke_test/bagle_fit_sanity.py`` performs a full joint photometric +
astrometric PSPL + parallax fit (scipy least-squares fallback if PyMultiNest
is unavailable) on selected well-behaved events, and compares recovered
parameters against GULLS truth.


Affected source files
---------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - File
     - Role
   * - ``src/structures.h``
     - Event structure: centroid vectors, PM structs, sky coordinate arrays
   * - ``src/readParamfile.cpp``
     - Parse ``ASTROMETRY_ON``, ``ASTROMETRIC_SYS_FLOOR``
   * - ``src/omLightcurveGenerator.cpp``
     - Alpha definition, (tau,u) to event-frame rotation, VBM centroid
       transform, multi-source flux-weighted blending
   * - ``src/photometry.cpp``
     - Lens blending, event-frame to ecliptic to RA/Dec pipeline, noise
   * - ``src/outputLightcurve.cpp``
     - Header metadata, column definitions, per-epoch output
   * - ``src/classes/parallax.cpp``
     - ``compute_tushifts``: parallax (tau, u) shifts (sign fix here)
   * - ``src/classes/coords.cpp``
     - ``ecl2ad``, ``muecl2ad``, ``muad2ecl``: coordinate transforms
   * - ``src/info.cpp``
     - ``.out`` file columns including PM bundles and pi_E components


Known limitations
-----------------

- **Ambient star positions**: ambient stars in the PSF aperture are blended
  into photometry but not yet into the astrometric centroid.  This will
  slightly bias the centroid for crowded fields.

- **Source-only sky columns**: no separate ``RA/Dec`` columns for source-only
  (blendless) centroids are published.  Source-only centroids are available
  only in theta_E (``blended_sources_only_x/y_thetaE``).

- **Multi-lens (N >= 3) VBM transforms**: assumed VBM output is already in
  the event frame.  This has not been independently verified for N >= 3
  configurations.

- **Commented-out overload**: ``parallax.cpp`` contains a commented-out
  ``compute_tushifts(vector<...>*)`` overload (lines 187--188) that still
  uses the old ``-piE`` convention for both tau and u.  If re-enabled, it
  must be updated.

- **Ephemeris mismatch**: GULLS uses Keplerian orbital elements for the
  observer position; BAGLE uses ERFA (IAU SOFA derivative).  This produces
  a ~0.1 microarcsecond astrometric floor in comparisons.
