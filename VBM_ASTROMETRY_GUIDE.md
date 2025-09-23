Here’s a compact, C++‑focused guide to VBMicrolensing’s astrometric APIs with just what your “add astrometry to GULLS” agent needs.

## Overview

* Provides astrometric centroid trajectories alongside magnification for PSPL, ESPL, binary, triple, binary‑source, and orbital/Xallarap cases.
* Outputs centroids in sky coordinates (North/East ≈ Dec/RA axes), plus lens/source trajectory coordinates in the lens-centered frame.

## Build & Setup

* Add `VBMicrolensing/lib/VBMicrolensingLibrary.cpp` to your target and include `VBMicrolensing/lib/VBMicrolensingLibrary.h` (if not already linked via your build). In GULLS, VBMicrolensing is already used for light curves, so this is likely in place.
* Compile with C++17. No external deps for core use.
* Place data files where your binary can load them:
  - `VBMicrolensing/data/ESPL.tbl` (for ESPL/extended source)
  - `VBMicrolensing/data/SunEphemeris.txt` (for parallax)
  - `Optional satellite ephemerides` (e.g., `VBMicrolensing/data/satellite1.txt` …)

## Initialization (parallax/coords)

* Set object sky coordinates (J2000) and load the Sun ephemeris before any “Parallax” or “Astro” calls:
  - `VBM.SetObjectCoordinates("RA Dec")` or `VBM.SetObjectCoordinates("coordfile.txt", "satdir")`
  - `VBM.LoadSunTable("SunEphemeris.txt")` (or set `VBM.parallaxephemeris=false` for internal ephemerides)
* Enable astrometry outputs: `VBM.astrometry = true`;
* For ESPL: `VBM.LoadESPLTable("ESPL.tbl")`;

## Astro APIs (array form; recommended)

* PSPL: `PSPLAstroLightCurve(pr, t, mags, c1s, c2s, c1l, c2l, y1, y2, np)` extending PSPLLightCurveParallax.
* ESPL: `ESPLAstroLightCurve(...)` extending ESPLLightCurveParallax.
* Binary: `BinaryAstroLightCurve(...)` extending BinaryLightCurveParallax.
* Orbital/Kepler: `BinaryAstroLightCurveOrbital(...)`, `BinaryAstroLightCurveKepler(...)`.
* Binary source / triple lens: `BinSourceAstroLightCurveXallarap(...)`, `BinSourceBinLensAstroLightCurve(...)`, `TripleAstroLightCurve(...)`.

### Outputs per point:

* `mags[i]`: magnification
* `c1s[i], c2s[i]`: source image centroid [mas] in N/E (≈ Dec/RA) axes
* `c1l[i], c2l[i]`: lens centroid [mas] in N/E axes
* `y1[i], y2[i]`: source position in the lens frame (Einstein radii)

## Parameter arrays and units (core cases)

* Times `t[]`, `t0`: days; use absolute JD (or the same absolute day system used in the Sun ephemeris) for consistent parallax geometry.
* PSPLAstroLightCurve pr[9]:
  - [0] u0 [—]
  - [1] log(tE) [days]
  - [2] t0 [days, JD]
  - [3] pi_N, [4] pi_E [—] parallax components (North/East)
  - [5] muS_Dec, [6] muS_RA [mas/yr] source heliocentric proper motions
  - [7] pi_S [mas] source parallax
  - [8] thetaE [mas] Einstein angle
* ESPLAstroLightCurve pr[10]:
  - PSPL entries with rho inserted: [0]=u0, [1]=log(tE), [2]=t0, [3]=log(rho), [4]=pi_N, [5]=pi_E, then [6..9]=muS_Dec,muS_RA,pi_S,thetaE
* BinaryAstroLightCurve pr[13] minimal (no orbital):
  - [0]=log(s), [1]=log(q), [2]=u0, [3]=alpha [rad], [4]=log(rho), [5]=log(tE),
  - [6]=t0, [7]=pi_N, [8]=pi_E, then [9..12]=muS_Dec,muS_RA,pi_S,thetaE
* All “Astro” variants add the 4 astrometry terms (muS_Dec, muS_RA, pi_S, thetaE) after their corresponding Parallax light‑curve parameter list.

### Notes:

* `c1 = North (Dec)` and `c2 = East (RA)` components; rotation from lens axes to sky is handled internally.
* Parallax components are North/East in the sky (equatorial basis).

## Blending (optional)

* Combine source/lens centroids with flux ratio g = FL/FS:
  - `VBM.CombineCentroids(mags, c1s, c2s, c1l, c2l, c1comb, c2comb, g, np);`

## Choosing the right API (important)

- Magnification‑only, single point (scalar return):
  - `PSPLMag(u)`, `ESPLMag2(u, rho)`, `BinaryMag2(s, q, y1, y2, rho)` …
  - Inputs are in “physical” form (e.g., `s`, `q`, `u`, `rho`), not logs; output is a single magnification.
  - With `VBM.astrometry=true`, single‑lens (`PSPLMag`, `ESPLMag2`) populate `VBM.astrox1/astrox2` in lens coordinates (Einstein radii). For binaries you still lack sky rotation, lens motion, and mas scaling.

- Light‑curve arrays (photometry):
  - `...LightCurve(...)` and `...LightCurveParallax(...)` return arrays and take the standard parameter vector with logs for some entries (e.g., `log s`, `log q`, `log tE`, `log rho`).

- Astro light‑curve arrays (recommended for sky centroids):
  - `...AstroLightCurve(...)` compute magnifications and sky centroids in mas, including parallax, proper motions, and rotation to N/E.
  - These expect the array‑style “standard” parameter vector (logs where specified) plus the extra astrometry terms. They are not drop‑in replacements for single‑point `...Mag*` calls.

Practical guidance for GULLS:
- Keep your existing photometry path (e.g., using `BinaryMag2` per epoch) if desired.
- For astrometry, add a parallel call to the corresponding `...AstroLightCurve(...)` over the same time grid to obtain `c1s/c2s` and `c1l/c2l` in mas, using the parameter builders provided above.
- Alternatively, switch photometry to the array `...LightCurveParallax(...)` family for consistency (same parameterization/time grid) and use the matching `...AstroLightCurve(...)` for centroids.


## Minimal PSPL astrometry example

* Shows required calls and arrays; plug into your pipeline.
* Include `VBMicrolensingLibrary.h`; link `VBMicrolensingLibrary.cpp`.
* Ensure `SunEphemeris.txt` and `ESPL.tbl` (if ESPL) are accessible.

### Example (abbreviated):

* `VBMicrolensing VBM;`
* Set coords: 
  `VBM.SetObjectCoordinates("17:51:40.2082 -29:53:26.502");`
* Load Sun ephemerides: 
  `VBM.LoadSunTable("SunEphemeris.txt");`
* Turn on: 
  `VBM.astrometry=true;`
* Fill `pr` (see PSPLAstro layout above), `t[]` (JD), allocate outputs.
* Call `VBM.PSPLAstroLightCurve(pr, t, mags, c1s, c2s, c1l, c2l, y1, y2, np);`

## Advanced controls (when needed)

* Extended source tables: `VBM.LoadESPLTable("ESPL.tbl");`
* Limb darkening: `VBM.SetLDprofile(VBMicrolensing::LDquadratic)` or custom via function pointer.
* Multiple lenses method: `VBM.SetMethod(VBMicrolensing::Method::Multipoly|Nopoly|Singlepoly)`.
* Binary lens photometric center: `VBM.lens_mass_luminosity_exponent` (default 4.0); force dark companion: `VBM.turn_off_secondary_lens = true;`
* Satellite view: provide sat tables and `VBM.satellite = 1|2`.

## Integrating with GULLS (quick mapping)

* Time base: VB expects absolute days consistent with the Sun ephemeris (typically JD). In GULLS, `t0` is relative to `SIMULATION_ZERO_TIME`. Use `t0_abs = simulation_zerotime + Event->t0` and pass times `t[i]` similarly offset to absolute JD.
* Parallax components: GULLS stores `piEN`/`piEE` in ecliptic N/E (reference‑frame). VB expects sky North/East (equatorial). Convert to equatorial NE before putting into `pr[pi_N]`, `pr[pi_E]`. GULLS’ `coords` helpers and `parallax` class can transform between bases (use the event’s RA/Dec).
* Proper motions: GULLS catalogs are Galactic `MUL`/`MUB`. Convert the source heliocentric proper motions to equatorial (Dec, RA) [mas/yr] for `muS_Dec`, `muS_RA`.
* Einstein angle: use `Event->thE` [mas] for `thetaE`.
* Einstein time: use the reference‑frame timescale (days) for `tE`; VB uses it in the light‑curve kinematics and in combining with `piE`.
* Extended source: pass `rho = Event->rs`.
* Binary: map `s`, `q`, `alpha(rad)`, `u0` directly; alpha in radians.

## Common pitfalls

* Forgetting `VBM.astrometry = true` (centroids remain unset).
* Passing relative days instead of JD (parallax geometry becomes inconsistent).
* Missing `SetObjectCoordinates`/`LoadSunTable` before any Parallax/Astro call.
* ESPL without loading `ESPL.tbl`.
* Mixing frames: pi components and proper motions must be equatorial N/E (not Galactic, not ecliptic).

## Drop‑in GULLS example: build pr arrays + convert frames

```cpp
#include "headers/coords.h"     // GULLS coord transforms
#include "src/constdefs.h"      // TO_RAD, TO_DEG (or define your own)
#include "VBMicrolensing/lib/VBMicrolensingLibrary.h"

// Convert (piEN, piEE) in ecliptic N/E to equatorial N/E at (ra,dec)
static inline void ecl_to_equ_pi(double ra, double dec,
                                 double piEN, double piEE, double piE,
                                 double* piN_equ, double* piE_equ)
{
  if (piE <= 0) { *piN_equ = 0; *piE_equ = 0; return; }
  coords c;
  double mulam = piEE / piE; // ecliptic East unit component
  double mubet = piEN / piE; // ecliptic North unit component
  double mua, mud;           // equatorial RA/Dec unit components
  c.muecl2ad(ra, dec, mulam, mubet, &mua, &mud);
  *piN_equ = piE * mud;      // Dec (North)
  *piE_equ = piE * mua;      // RA (East)
}

// Convert source proper motion from Galactic (mul,mub) to equatorial (RA,Dec)
static inline void gal_to_equ_mu(double l_deg, double b_deg, double mul, double mub,
                                 double* muRA, double* muDec)
{
  coords c;
  c.mulb2ad(l_deg*TO_RAD, b_deg*TO_RAD, mul, mub, muRA, muDec);
}

// Build PSPL Astro parameter vector (size 9)
static inline void build_pspl_pr(const filekeywords* P, const event* E,
                                 const slcat* Sources, double pr[9])
{
  // Times
  double t0_abs = P->simulation_zerotime + E->t0;        // absolute days (JD)
  double tE_ref = E->tE_r;                                // days

  // Parallax components (equatorial N/E)
  double piN=0, piE=0;
  ecl_to_equ_pi(E->ra, E->dec, E->piEN, E->piEE, E->piE, &piN, &piE);

  // Source proper motions (equatorial) from catalog Galactic components
  int sn = E->source;
  double muRA=0, muDec=0;
  gal_to_equ_mu(E->l, E->b, Sources->data[sn][Sources->MUL],
                Sources->data[sn][Sources->MUB], &muRA, &muDec);

  // Source parallax in mas from distance in kpc
  double piS = 1000.0 / Sources->data[sn][Sources->DIST];

  // Fill pr array (PSPLAstro)
  pr[0] = E->u0;            // u0
  pr[1] = log(tE_ref);      // log tE
  pr[2] = t0_abs;           // t0 (absolute)
  pr[3] = piN;              // pi_N (equatorial North)
  pr[4] = piE;              // pi_E (equatorial East)
  pr[5] = muDec;            // muS_Dec [mas/yr]
  pr[6] = muRA;             // muS_RA  [mas/yr]
  pr[7] = piS;              // source parallax [mas]
  pr[8] = E->thE;           // thetaE [mas]
}

// Build ESPL Astro parameter vector (size 10)
static inline void build_espl_pr(const filekeywords* P, const event* E,
                                 const slcat* Sources, double pr[10])
{
  double pr_pspl[9];
  build_pspl_pr(P,E,Sources,pr_pspl);
  pr[0] = pr_pspl[0];             // u0
  pr[1] = pr_pspl[1];             // log tE
  pr[2] = pr_pspl[2];             // t0
  pr[3] = log(E->rs);             // log rho
  pr[4] = pr_pspl[3];             // pi_N
  pr[5] = pr_pspl[4];             // pi_E
  pr[6] = pr_pspl[5];             // muS_Dec
  pr[7] = pr_pspl[6];             // muS_RA
  pr[8] = pr_pspl[7];             // pi_S
  pr[9] = pr_pspl[8];             // thetaE
}

// Build Binary Astro parameter vector (size 13); requires s,q,alpha
static inline bool build_binary_pr(const filekeywords* P, const event* E,
                                   const slcat* Sources, double pr[13])
{
  // Require s and q present in Event->params
  // SS = separation in Einstein radii, QQ = mass ratio
  // See src/columnCodes.h for indices
  #ifndef SS
  #define SS (NPLANETINPUT + 1)
  #define QQ (NPLANETINPUT + 0)
  #endif
  if (E->params.size() == 0) return false;

  double pr_pspl[9];
  build_pspl_pr(P,E,Sources,pr_pspl);

  pr[0] = log(std::max(1e-12, E->params[SS])); // log s
  pr[1] = log(std::max(1e-12, E->params[QQ])); // log q
  pr[2] = E->u0;                               // u0
  pr[3] = E->alpha * TO_RAD;                   // alpha (rad) from degrees
  pr[4] = log(std::max(1e-12, E->rs));         // log rho
  pr[5] = pr_pspl[1];                          // log tE
  pr[6] = pr_pspl[2];                          // t0
  pr[7] = pr_pspl[3];                          // pi_N
  pr[8] = pr_pspl[4];                          // pi_E
  pr[9] = pr_pspl[5];                          // muS_Dec
  pr[10]= pr_pspl[6];                          // muS_RA
  pr[11]= pr_pspl[7];                          // pi_S
  pr[12]= pr_pspl[8];                          // thetaE
  return true;
}

// Usage sketch
void run_pspl_astro(const filekeywords* P, const obsfilekeywords* W,
                    const event* E, const slcat* Sources)
{
  VBMicrolensing VBM;
  // Coordinates must be set once (string form or via file); ensure Sun ephemeris loaded
  // VBM.SetObjectCoordinates("RA Dec");
  // VBM.LoadSunTable("SunEphemeris.txt");
  VBM.astrometry = true;

  // Build parameter vector
  double pr[9];
  build_pspl_pr(P,E,Sources,pr);

  // Times: pass absolute JD times; here we reuse GULLS per‑obs epochs
  const std::vector<double>& t = E->jdtimes[E->obsidx[0]]; // example: first obs stream
  int np = (int)t.size();
  std::vector<double> mags(np), c1s(np), c2s(np), c1l(np), c2l(np), y1(np), y2(np);

  VBM.PSPLAstroLightCurve(pr, const_cast<double*>(t.data()), mags.data(),
                          c1s.data(), c2s.data(), c1l.data(), c2l.data(),
                          y1.data(), y2.data(), np);

  // c1s/c2s and c1l/c2l are astrometric centroids in mas
}
```

This snippet shows the two critical conversions (ecliptic→equatorial for parallax; Galactic→equatorial for source μ), builds the parameter arrays in VB’s expected order, and calls the PSPL Astro function with absolute JD times.
