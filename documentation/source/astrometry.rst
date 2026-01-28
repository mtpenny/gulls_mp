Astrometry outputs and parameters
=================================

This page summarizes the astrometry options added to gulls, the new output columns, and the noise recipe used to generate observed astrometric positions.

Parameters
----------

- ``ASTROMETRY_ON`` (default: ``0``)
  - Enable astrometric computation and outputs when set to 1.
  - **Warning**: Enabling astrometry significantly slows VBMicrolensing computations.
  - When 0, all astrometry outputs are written as 0.0.

- ``ASTROMETRIC_SYS_FLOOR`` (units: mas, default: ``0.1``)
  - Per-axis systematic floor for astrometric uncertainty.
  - Combined in quadrature with the photon-limited term when producing per-epoch errors.

Coordinate Systems
------------------

Three coordinate systems are relevant to astrometry:

1. **VBM Internal Frame** (units: θ_E)
   - The coordinate system used internally by VBMicrolensing library.
   - Varies by lens configuration:
     - *Single lens*: x1 along source-lens axis, x2 = 0 by axial symmetry
     - *Binary lens*: x1 along binary axis, x2 perpendicular, origin at center of mass
     - *N-lens (N≥3)*: Same as input source coordinates
   - Raw VBM outputs are stored for debugging coordinate transforms.

2. **Lens-Centered Coordinates** (units: θ_E, converted to mas for output)
   - Simulation coordinate system used throughout gulls.
   - **Origin**: Lens center of mass (stationary in this frame).
   - **Units**: Angular Einstein radius (θ_E).
   - **x-axis**: Orientation relative to celestial East is **UNCERTAIN** - depends on α convention.
   - **y-axis**: Orientation relative to celestial North is **UNCERTAIN**.
   - Orientation defined by trajectory angle α, which may be inconsistently defined.
   - All centroid columns (except raw VBM) are in this frame.
   - **Validate x,y orientation with plots before trusting RA/Dec output!**

3. **ICRS Equatorial Coordinates** (RA/Dec in degrees)
   - Absolute positions in the International Celestial Reference System.
   - Base RA/Dec is written to lightcurve header as ``#Astrometry_Frame``.
   - **Two versions output** for validation:
     
     a. Without lens parallax (``_deg`` columns): assumes catalog RA/Dec is lens position
     b. With lens parallax (``_lpllx_deg`` columns): attempts to correct for observer position
   
   - Conversion assumes x~East, y~North (may be wrong!)::
   
         RA_deg = RA_base_deg + x_mas / (3600 × 1000 × cos(Dec_base))
         Dec_deg = Dec_base_deg + y_mas / (3600 × 1000)
   
   - **Does NOT account for lens proper motion** - base RA/Dec is catalog position at t_ref.

Output Columns
--------------

All centroid columns are in **mas** (milliarcseconds) unless noted otherwise.

Raw VBM Output (for debugging)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These columns store the raw VBMicrolensing library output without transformation:

- ``vbm_astrox1_raw_thE`` (θ_E): Raw VBM x1 centroid output.
- ``vbm_astrox2_raw_thE`` (θ_E): Raw VBM x2 centroid output.

**Use case**: Debugging coordinate transforms between VBM internal coordinates and lens-centered coordinates.

Flux-Weighted Blending Columns
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The photocenter is computed in stages, with each stage stored for diagnostic purposes:

1. **Source-only centroid** (lensed sources, no baseline blending):

   - ``centroid_src_x_mas``: Flux-weighted centroid of lensed source images (x).
   - ``centroid_src_y_mas``: Flux-weighted centroid of lensed source images (y).

   This is after:
   
   - VBM centroid transformed from VBM internal coordinates to lens-centered coordinates
   - Flux-weighted blend across multiple sources (if present)

2. **Source + lens centroid** (after blending with luminous lenses):

   - ``centroid_src_lens_x_mas``: Centroid after blending with lens 1 and lens 2 (if luminous).
   - ``centroid_src_lens_y_mas``: Same for y-component.

3. **Final centroid** (after blending with ambient stars):

   - ``centroid_final_x_mas``: Final blended centroid (x).
   - ``centroid_final_y_mas``: Final blended centroid (y).

   This equals the true centroid (``true_x_centroid_mas``).

True and Measured Centroids
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lens-centered coordinates (mas):

- ``true_x_centroid_mas``: Final blended centroid (noise-free).
- ``true_y_centroid_mas``: Final blended centroid (noise-free).
- ``true_x_centroid_error_mas``: Always 0.0 (no uncertainty on true value).
- ``true_y_centroid_error_mas``: Always 0.0.

- ``x_centroid_mas``: Measured centroid with Gaussian noise added.
- ``y_centroid_mas``: Measured centroid with Gaussian noise added.
- ``x_centroid_error_mas``: 1-σ uncertainty on measured centroid.
- ``y_centroid_error_mas``: 1-σ uncertainty on measured centroid.

Equatorial coordinates (ICRS, degrees) - **WARNING: x,y → E,N mapping is uncertain!**

Without lens parallax:

- ``RA_centroid_deg``: Observed centroid RA (with noise), relative to catalog position.
- ``Dec_centroid_deg``: Observed centroid Dec (with noise).
- ``RA_centroid_true_deg``: True centroid RA (noise-free).
- ``Dec_centroid_true_deg``: True centroid Dec (noise-free).

With lens parallax attempt:

- ``RA_centroid_lpllx_deg``: Observed RA with lens parallax correction.
- ``Dec_centroid_lpllx_deg``: Observed Dec with lens parallax correction.
- ``RA_true_lpllx_deg``: True RA with lens parallax correction.
- ``Dec_true_lpllx_deg``: True Dec with lens parallax correction.

Validation data:

- ``lens_dist_kpc``: Lens distance in kpc (for computing your own parallax corrections).
- ``lens_parallax_x_mas``: Lens parallax shift in x (mas) using observer position and lens distance.
- ``lens_parallax_y_mas``: Lens parallax shift in y (mas) using observer position and lens distance.

Source and Lens Positions
~~~~~~~~~~~~~~~~~~~~~~~~~

Per-epoch positions of all sources and lenses are output in lens-centered coordinates (θ_E):

- ``sourceN_x_thE``, ``sourceN_y_thE``: Position of source N (N=0,1,...).
- ``sourceN_mu``: Magnification of source N.
- ``lensN_x_thE``, ``lensN_y_thE``: Position of lens N (N=0,1,...).

These positions are in lens-centered coordinates. **The x,y orientation relative to celestial E,N is uncertain!**
Origin is at lens center of mass (lens 0 is at (0,0) for single lens events).

**Validation tips**:

1. Convert source/lens positions to mas (multiply by θ_E), compare against centroids.
2. The centroid should shift toward the lens during high magnification.
3. Plot VBM raw outputs vs lens-centered centroids to verify coordinate transforms.
4. Compare RA/Dec output (both with and without lens parallax) against existing fitting codes.

Blending Model
--------------

Flux-weighted centroid blending is applied in stages:

1. **Multiple sources** (in ``omLightcurveGenerator.cpp``):
   
   For each source, VBM computes the flux-weighted image centroid. If multiple sources exist, they are blended using magnified fluxes as weights::

       centroid = Σ(flux_i × magnification_i × centroid_i) / Σ(flux_i × magnification_i)

2. **Luminous lenses** (in ``photometry.cpp``):
   
   The source centroid is blended with up to two luminous lens positions::

       centroid = (centroid_src × f_src + pos_lens1 × f_lens1 + pos_lens2 × f_lens2) / (f_src + f_lens1 + f_lens2)

   Flux fractions ``fl1``, ``fl2`` should be set in ``buildEvent.cpp`` (currently defaults to old behavior if not set).

3. **Ambient stars** (in ``photometry.cpp``):
   
   Blend with ambient stars in the aperture. Currently uses lens position as proxy (TODO: proper positions).

Flux Fractions
~~~~~~~~~~~~~~

The following flux fractions are tracked per observatory:

- ``fs``: Source(s) flux fraction (including companion sources)
- ``fl1``: Primary lens flux fraction
- ``fl2``: Secondary lens flux fraction (if luminous)
- ``famb``: Ambient stars flux fraction

These should satisfy: ``fs + fl1 + fl2 + famb = 1``

**Note**: Proper computation of ``fl1``, ``fl2``, ``famb`` from lens magnitudes requires implementation in ``buildEvent.cpp``. Currently, if these are not set, all non-source flux is assumed to be at the primary lens position (legacy behavior).

Noise Model
-----------

Per-epoch astrometric uncertainties are computed following Gould & Yee (2014).
All noise calculations are performed entirely in mas:

1. Fractional photometric error::

       σ_phot = Aerr / Aobs

2. PSF FWHM converted to mas::

       FWHM_mas = FWHM_arcsec × 1000

3. Photon-limited astrometric uncertainty (mas)::

       σ_astro_mas = FWHM_mas × σ_phot / √(ln 256)

4. Total uncertainty with systematic floor (mas)::

       σ_total_mas = √(σ_astro_mas² + ASTROMETRIC_SYS_FLOOR²)

5. Convert final blended centroid to mas::

       xctrue_mas = xctrue_thE × θ_E

6. Apply noise to final blended centroid (mas)::

       xc_mas = xctrue_mas + σ_total_mas × N(0,1)

7. Storage and output:

   - ``xc``, ``yc`` (observed): stored and output in **mas**
   - ``xcerr``, ``ycerr``: stored and output in **mas**
   - ``xctrue``, ``yctrue`` (true): stored in θ_E, output as ``× θ_E`` → mas

Affected Files
--------------

- ``src/structures.h``: Astrometry parameters and intermediate centroid vectors.
- ``src/readParamfile.cpp``: Parse ``ASTROMETRY_ON`` and ``ASTROMETRIC_SYS_FLOOR``.
- ``src/omLightcurveGenerator.cpp``: Compute VBM centroids, transform to lens-centered coordinates, blend multiple sources.
- ``src/photometry.cpp``: Blend with luminous lenses and ambient stars, add noise.
- ``src/outputLightcurve.cpp``: Write astrometry columns to output files.

Lightcurve Header
-----------------

The lightcurve file contains astrometry metadata in the header:

.. code-block:: text

   #Astrometry_Frame: RA_rad=X Dec_rad=Y RA_deg=X Dec_deg=Y thE_mas=Z t0=T
   #Astrometry_Frame: origin=catalog_position, xy_orientation=UNCERTAIN(validate!)

This defines the coordinate origin (catalog lens position) in both radians and degrees,
along with the angular Einstein radius (θ_E in mas) needed for unit conversions.

Known Limitations
-----------------

- Ambient star positions: Currently uses lens position as proxy.
- Lens flux fractions: Requires ``buildEvent.cpp`` implementation for proper values.
- Multi-lens (N≥3) coordinate transforms: Assumed VBM output is already in lens-centered coordinates.
