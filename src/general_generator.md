## A new lightcurve generator has been developed to handle multi-lens and multi-source lightcurves with orbital motion - this is called general

# Limitations

- Can only handle up to binary stars as source or lens hosts, if higher multiples are provided then the behaviour is undefined (segfaults are likely)
- Can handle multiplanet systems or single-planet, multiple-moon systems in single or binary star systems. Can only put planets around one of the stars. No planets in source star systems
- Secondary sources can be dark

#Conventions
- A single star from the binary (source or lens) is chosen as the main (need not be the most massive or most bright). This is the host for wide-binary lenses. All quantities (e.g.,mass ratio) are relative to this object
- Orbits are defined by the vector addition of Keplerian ellipses, with objects internal to the currently considered object experiencing reflex orbital motion. It is the user's responsibility to ensure stability of the system in postprocessing
- Semimajor axis refers to the combined semimajor axis
- Input periods may be modified by the presence of additional masses

#Input formats
