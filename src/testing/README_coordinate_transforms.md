# Coordinate transform tests

From the repository root, run:

```sh
python3 -m pip install -r src/testing/requirements_coordinates.txt
python3 src/testing/test_coordinate_transforms_astropy.py
```

The runner compiles `coordinate_probe.cpp` with the production C++ sources and sends test cases through that executable. It needs a C++17 compiler and GSL development headers/libraries (`g++ libgsl-dev` on Ubuntu; `brew install gsl` on macOS). It does not configure or build the Fortran simulator.

The Astropy check covers every public position and proper-motion conversion in `coords`: Galactic/equatorial and equatorial/ecliptic positions (pointer and vector overloads), all four proper-motion directions and their cosine-latitude components, `XYZ`/`Rlb`, and `UVW` projections. It also checks the angle-unit and `fold` helpers. The probe covers the other coordinate converters in `astroFns`: `eq2gal`, `eq2eclip`, `gal2eclip`, `eq2horiz`, `lst`, and `geod2geoc`, plus the production consumer `parallax::set_lb`. The orbital state calculators such as `suncoo` and `mooncool` calculate object positions rather than convert coordinate frames. The separate `orbitalElements` transform helpers in `ephem.cpp` are not part of this interface and are outside this test's scope.

The frame choices and tolerances are explicit:

- The retained `eq2gal` matrix is compared in FK5/J2000 coordinates. Its rounded published coefficients differ from Astropy's FK5/J2000 transform by up to 9.2 mas in the all-sky sample. Comparing the same numeric coordinates as ICRS instead shifts them by about 25 mas because FK5 and ICRS have different axis realizations.
- Equatorial/ecliptic positions use `FK5(equinox=J2000)` and `BarycentricMeanEcliptic(equinox=J2000)`. Astropy's barycentric ecliptic frame is ICRS-based, while `coords` applies the fixed J2000 obliquity to the published equatorial coordinates; the measured frame offset is at most 28 mas in this sample, with a 35 mas bound.
- `eq2horiz` uses mean RA/Dec and local mean sidereal time. The test checks its sidereal time against Astropy's IAU 2006 mean sidereal time, then uses Astropy's `HADec` to `AltAz` local rotation. It avoids an FK5-to-observed-AltAz transform, which would add aberration and other apparent-place corrections absent from this legacy helper.
- `uvw2rvmuadlb` is checked against Astropy Cartesian velocities with U toward the Galactic center and V/W along Galactic +Y/+Z, distance in kpc, and velocity in km/s. The comparison applies the existing `4.74057` conversion factor explicitly; it is not replaced by Astropy's more exact unit constant.
- Galactic and ecliptic proper-motion components are compared directly away from their poles. At `|b|=89.9 deg`, the FK5/ICRS position offset rotates the Galactic longitude tangent basis enough to shift components by `3.4e-5 mas/yr`; near-pole behavior is instead checked with production finite differences, inverse conversions, and vector-norm preservation.

The first Galactic direction is inferred from the supplied table RA/Dec and exercises the actual `eq2gal` and `parallax::set_lb` paths. It does not claim to reproduce the historical event's `.out` inputs or its sigma light curves.
