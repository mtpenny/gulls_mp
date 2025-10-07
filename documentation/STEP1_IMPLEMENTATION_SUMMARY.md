# Step 1: VBM Centroid Extraction - Implementation Summary

## Date: October 7, 2025

## What Was Implemented

Added low-level astrometry centroid extraction to `src/pllxLightcurveGenerator.cpp` without using the high-level `BinaryAstroLightCurve` API.

### Changes Made

1. **Enabled VBM astrometry mode** (line ~51)
   - Added `Event->vbm->astrometry = true;` before the epoch loop
   - This tells VBMicrolensing to compute and store centroids in `astrox1` and `astrox2` after each `BinaryMag2` call

2. **Single source centroid extraction** (lines ~141-165)
   - After the `BinaryMag2` call for the primary source
   - Extract `Event->vbm->astrox1` and `Event->vbm->astrox2` (in Einstein radii, lens frame)
   - Convert to mas: multiply by `Event->thE` (Einstein angle in mas)
   - Rotate from lens frame to sky frame using trajectory angle `alpha`
   - Store in `Event->xctrue[idx]` (North offset, mas) and `Event->yctrue[idx]` (East offset, mas)

3. **Binary source centroid extraction** (lines ~102-139)
   - Get centroid for first source after first `BinaryMag2` call
   - Get centroid for second source after second `BinaryMag2` call
   - Compute flux-weighted centroid: `(flux1 * c1 + flux2 * c2) / (flux1 + flux2)`
   - Convert to mas and rotate to sky frame
   - Store flux-weighted result in `Event->xctrue[idx]` and `Event->yctrue[idx]`

## Coordinate System Conventions

### VBM Lens Frame
- Origin: Binary lens center of mass
- x1-axis: Along the binary axis (from primary to secondary)
- x2-axis: Perpendicular to binary axis
- Units: Einstein radii (θ_E)
- Stored in: `Event->vbm->astrox1`, `Event->vbm->astrox2`

### GULLS Lens Frame (used for trajectory)
- Origin: Binary lens center of mass (shifted by VBM_origin for VBM calculations)
- x-axis: Along source trajectory at angle `alpha` from binary axis
- y-axis: Perpendicular to trajectory
- Units: Einstein radii

### Sky Frame (output)
- Origin: Lens position (includes parallax offsets if enabled)
- x-axis: North direction
- y-axis: East direction  
- Units: milliarcseconds (mas)
- Stored in: `Event->xctrue[idx]`, `Event->yctrue[idx]`

## Rotation Transform

We rotate from lens frame (x1 along the binary axis, x2 perpendicular) to sky North/East using the sky orientation defined by the microlensing parallax vector πE = (πEN, πEE):

PosAng = atan2(πEE, πEN) − alpha  (+ dPosAng for orbital rotation if present)

Then apply the rotation to the centroid in mas:

```
N = cx_lens * cos(PosAng) + cy_lens * sin(PosAng)
E = -cx_lens * sin(PosAng) + cy_lens * cos(PosAng)
```

Where `cx_lens`, `cy_lens` are the centroid in the lens frame (converted to mas), and N, E are the offsets in the sky frame.

## Binary Source Flux Weighting

For binary sources (xallarap case):
- `flux1 = amp` (magnification of source 1)
- `flux2 = fsofs1 * amp2` (magnification contribution of source 2)
- Weighted centroid: `(flux1 * c1 + flux2 * c2) / (flux1 + flux2)`

Note: The total observed magnification is `amp + fsofs1 * (amp2 - 1)`, accounting for the baseline flux.

## What's NOT Yet Implemented

1. **Centroid uncertainties** - The `xctrueerr` and `yctrueerr` vectors are not yet populated
2. **Observational noise** - The `xc`, `yc`, `xcerr`, `ycerr` vectors (with noise) are not yet populated
3. **Parallax effects on centroids** - Currently just rotates to sky frame; parallax offsets not yet applied to centroid positions
4. **Lens luminosity contribution** - If lens is luminous, its centroid should be flux-weighted with source centroid
5. **Coordinate transformations** - No conversion of parallax/proper motion between Galactic/Ecliptic/Equatorial frames yet

## Next Steps (in order)

1. **Step 2**: Implement noise model in `photometry.cpp`
   - Calculate photon noise: σ = FWHM/SNR  
   - Add systematic floor (0.1 mas)
   - Generate Gaussian perturbations
   - Populate `xc`, `yc`, `xcerr`, `ycerr` vectors

2. **Step 3**: Add coordinate transformation functions in `astroFns.cpp`
   - `mulb2ad`: Galactic (μ_l, μ_b) → Equatorial (μ_α, μ_δ)
   - `muecl2ad`: Ecliptic (π_EN, π_EE) → Equatorial (π_N, π_E)

3. **Step 4**: Add RA/Dec absolute positions in output
   - Convert mas offsets to degrees
   - Propagate lens proper motion from t0
   - Handle RA wrap-around at 0°/360°

## Testing/Validation Notes

- Cannot easily compile due to build system complexity
- Once built, test with simple single-lens event first
- Compare with VBM Python `BinaryAstroLightCurve` output for validation
- Check that centroids are ~0 at high magnification (source approaches center)
- Check that centroid traces roughly follow source trajectory at baseline

## References

- VBMicrolensing header: `src/headers/VBMicrolensingLibrary.h` line 158 (astrox1, astrox2 members)
- VBM documentation: `VBMicrolensing/docs/C++/CentroidTrajectories.md`
- Implementation plan: `documentation/VBM_ASTROMETRY_PLAN.md`
- Event structure: `src/structures.h` lines 374-381 (centroid vectors)
