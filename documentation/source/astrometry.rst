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

Coordinate systems and column names
-----------------------------------

Lens-frame (VBM) centroid
~~~~~~~~~~~~~~~~~~~~~~~~~

- ``true_x_centroid`` (Einstein radii): blended centroid x1 in the VBM lens frame (x1 along the binary axis).
- ``true_y_centroid`` (Einstein radii): blended centroid x2 in the VBM lens frame (x2 perpendicular to the binary axis).
- ``true_x_centroid_err`` (Einstein radii): uncertainty on ``true_x_centroid``. Value is 0.0 for all epochs.
- ``true_y_centroid_err`` (Einstein radii): uncertainty on ``true_y_centroid``. Value is 0.0 for all epochs.

Measurement uncertainty/noise approximations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``x_centroid`` (Einstein radii): blended centroid x1 in the VBM lens frame with a random noise component.
- ``y_centroid`` (Einstein radii): blended centroid x2 in the VBM lens frame with a random noise component.
- ``x_centroid_err`` (Einstein radii): dispersion of the noise on ``x_centroid``.
- ``y_centroid_err`` (Einstein radii): dispersion of the noise on ``y_centroid``.

Sky NE (North/East) centroid — lens-centric offsets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Not yet implemented.

Absolute RA/Dec centroid columns
--------------------------------

Not yet implemented.

Notes on frames and units:
- ``true_x_centroid``/``true_y_centroid`` are in the lens frame (Einstein radii).

Noise model and zeroing rules
-----------------------------

Per-epoch astrometric uncertainties and observed values are generated as follows:

1. Compute a fractional photometric error ratio per epoch from the simulated photometry:
   - ``sigma_phot = Aerr / max(Aobs, 1e-12)``
2. PSF width is taken from the instrument model as ``FWHM`` in arcsec; we convert to mas via ``FWHM_mas = 1000 * FWHM``.
3. Convert to Einstein-radius units using the event's ``thetaE_mas``:
   - ``FWHM_er = FWHM_mas / thetaE_mas``
4. Photon-limited astrometric uncertainty per axis per Gould & Yee (2014):
   - ``sigma_astro = FWHM_er * sigma_phot / sqrt(ln(256))``
5. Total per-axis uncertainty combines the photon term with a systematic floor:
   - ``sigmaAstro = sqrt(sigma_astro^2 + (ASTROMETRIC_SYS_FLOOR/thetaE_mas)^2)``

Zeroing behavior (no sky orientation or disabled):
- If ``ASTROMETRY_ON = 0``, NE and absolute RA/Dec astrometric outputs are disabled (set to 0.0).

Affected files
--------------

- ``src/structures.h`` : added parameters to ``Paramfile`` structure.
- ``src/readParamfile.cpp`` : read new parameters from the parameter file.
- ``src/pllxLightcurveGeneratorMultiple.cpp`` : set true astrometric values during lightcurve generation.
- ``src/pllxLightcurveGenerator.cpp`` : set true astrometric values during lightcurve generation.
- ``src/photometry.cpp`` : compute observed astrometric values and uncertainties during photometry step.
- ``timeSeriesOutput.cpp`` : write new astrometric columns to output files.
