> A new lightcurve generator has been developed to handle multi-lens and multi-source lightcurves with orbital motion - this is called general

## Limitations

- Can only handle up to binary stars as source or lens hosts, if higher multiples are provided then the behaviour is undefined (segfaults are likely)
- Can handle multiplanet systems or single-planet, multiple-moon systems in single or binary star systems. Can only put planets around one of the stars. No planets in source star systems
- Secondary sources can be dark

## Conventions
- A single star from the binary (source or lens) is chosen as the main lens (need not be the most massive or most bright). This is the host for wide-binary lenses. All relative quantities (e.g.,mass ratio) microlensing parameters (thetaE, rE, etc.) are relative to this object
- Orbits are defined by the vector addition of Keplerian ellipses, with objects internal to the currently considered object experiencing reflex orbital motion. It is the user's responsibility to ensure stability of the system in postprocessing
- Semimajor axis refers to the combined semimajor axis
- Input periods may be modified by the presence of additional masses

## Input formats

Binary stars are included via SynthPop's conventions for binary stars; either the primary or secondary star is drawn from the catalog, and its companion brought along with it. Physical parameters of the binary orbit are provided in the SynthPop catalog (period, eccentricity), others are free to be generated in gulls (currently eccentricity defaults to zero). Planets are supplied through planets files that list the following parameters:
```
Mass SemimajorAxis Eccentricity Inclination LongitudePerihelion LongitudeAscNode OrbitType
```
with angles in degrees. If inclination has a value >900, its inclination will be referenced to the binary star orbit's as `I = I_binary + (I-1000.0)`; if the star is not a binary, the inclination will be `I-1000`. OrbitType is an integer code, with a value of 1 or 2 indicating a planet, and 3 indicating a moon. Use multiple sets of these parameters for multiple planets/moons.
