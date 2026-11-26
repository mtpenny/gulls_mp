# Upgrading the vendored VBMicrolensing library

`src/classes/VBMicrolensingLibrary.cpp` and `src/headers/VBMicrolensingLibrary.h` are a
**customized** copy of [valboz/VBMicrolensing](https://github.com/valboz/VBMicrolensing)
(`VBMicrolensing/lib/`). Do not overwrite them with a fresh upstream copy: that would
silently drop the gulls changes listed below.

## Current state

| | |
|---|---|
| Upstream version | tag `v5.6.1` (commit `db4f952`) |
| Merge base for the next upgrade | upstream tag `v5.6.1` (pass the tag name; the script below normalizes it) |

## gulls customizations (what must survive an upgrade)

1. **`complex` → `VBcomplex` rename** everywhere (class, locals, signatures). This avoids
   a clash with gulls' own complex type. A whole-word replace reproduces it exactly:
   `perl -pe 's/\bcomplex\b/VBcomplex/g'`.
2. **Opt-in timeouts and error policy** (`VBMTimeoutError`, `TimeoutConfig`,
   `SetTimeouts`, `ErrorPolicy`, `LastError`, plus `TimeBudget`, `ScopedBudget`,
   `ShouldCheck` and `CheckTimeout` in an anonymous namespace). Public entry points are
   wrapped in `try { TimeBudget …; ScopedBudget …; … } catch (const VBMTimeoutError&)`,
   and long loops call `ShouldCheck/CheckTimeout`. Timeouts are disabled by default.
   New public functions that upstream adds need the same wrapper: in 5.6 these were
   `TripleAstroLightCurveOrbital` (Astrometry), `TripleLightCurveOrbital`
   (Magnification) and `t0_from_t0_par` (Parallax). Each wrapper opens with
   `ApiScope api_scope(this);` (declared before the `try`; upstream has nothing there).
   It counts nesting depth: only the outermost call clears `LastError`, and
   `HandleTimeoutError` rethrows whenever the call is nested, so the policy is applied
   once by the function the caller actually invoked. Without it, an inner wrapper
   swallowed the timeout under `ReturnNaN` and the next inner call cleared `LastError`,
   so the outer call could return a wrong value with no error set. Functions that catch
   but do not clear (`MultiMag`, `ComputeParallax`) also take an `ApiScope`. Both
   `BinaryMag0_shear` overloads have no upstream timeout handling, so the 8-argument one
   carries a gulls wrapper (Magnification budget, returns NaN with `*Images == nullptr`).
   Without it, the root solver was the outermost call, so under `ReturnNaN` it absorbed
   the timeout and the shear magnification was computed from incomplete roots. The private
   helper `du_par` still calls `CheckTimeout` after each `ComputeParallax` call. That is
   now redundant, but harmless.
   See `VBM_timeout_only_edits.md`. The wrappers re-indent function bodies, so a plain
   3-way merge produces many whitespace-only conflicts (see below).
3. **RAII ownership** in `BinaryMag`, `BinaryMag0`, `MultiMag`, `MultiMag0` and the
   `*Dark` functions, plus both `PlotCrit` overloads. These use `std::unique_ptr`
   (`images_owner`, `thetas_owner`, `stheta_owner`, `prov_owner`, `AnnulusChainPtr`;
   in `PlotCrit`, `CriticalCurves`, `samples` and `caustic`) instead of raw
   `new`/`delete`, so a timeout that exits early frees what was built so far. Upstream
   `PlotCrit` also leaks its per-step sample curve on every successful call; the
   `samples` owner fixes that too. `Prov` lists are
   drained with `while (Prov->first) { …; Prov->drop(scan1); }`. On upstream's own
   `return -1` paths, `MultiMag` still hands back the image list
   (`images_owner.release()`), which preserves upstream's contract with `MultiMagSafe`.
   `BinaryMag` does not: its `return -1` paths (no images found on the contour) leave
   `*Images == nullptr`. A timeout never reaches the callers' null checks (nested calls
   rethrow), so a null list there means VBM's own calculation failed, for example
   `BinaryMagSafe`'s last retry for a wide planet whose source covers the host. The
   callers (`BinaryMag2`, `BinaryMagDark`, `MultiMag2`, `MultiMagDark`) raise it with
   `ThrowMissingImages`, an `Unknown`-category `VBMTimeoutError`, so the outermost call
   sets `LastError` and returns NaN (Gulls maps `Unknown` to `LCGEN_VBM_ERR` and
   discards the event). Before this, they returned NaN with no error set and Gulls wrote
   all-NaN lightcurves. Restoring upstream's behaviour instead would return a finite
   magnification there, but upstream is inaccurate in that regime and its
   `BinaryMagSafe` leaves the centroid from its last, failed retry. Every `annulus` is allocated
   value-initialized (`new annulus()`), so the `AnnulusChainPtr` deleter never follows an
   uninitialized `next` when a timeout exits early. `BinaryMag2` and `MultiMag2` likewise
   raise the same error if `BinaryMag0`/`MultiMag0` hand back `Images == nullptr`. Upstream uses
   `new annulus;` and does not need the null checks, so re-check both after a merge.
4. **`max_num_images`**, a public counter updated in `MultiMag2`/`MultiMagDark` and used
   by gulls diagnostics.
5. **Diagnostics**: the extended `"Source too large!"` message in `ESPLMag`.
6. The unused `dist_mp` member is removed.

## Upgrade procedure (3-way merge)

Use a merge with the upstream base, not a hand diff. Put scratch files in a private
workspace, not in the repo.

```bash
W=$(mktemp -d); cd $W
git clone -q https://github.com/valboz/VBMicrolensing.git up
norm() { tr -d '\r' | perl -pe 's/\bcomplex\b/VBcomplex/g'; }
OLD=v5.6.1     # upstream version the current gulls copy is based on (see table above)
NEW=v5.7       # target upstream tag. Upstream's default branch is `main`; NEW=origin/dev previews unreleased changes
G=/path/to/gulls

mkdir mr && cd mr && git init -q
for f in cpp h; do git -C ../up show $OLD:VBMicrolensing/lib/VBMicrolensingLibrary.$f | norm > VBMicrolensingLibrary.$f; done
git add . && git commit -qm base && B=$(git rev-parse HEAD)
git checkout -qb theirs
for f in cpp h; do git -C ../up show $NEW:VBMicrolensing/lib/VBMicrolensingLibrary.$f | norm > VBMicrolensingLibrary.$f; done
git commit -qam theirs
git checkout -qb ours $B
tr -d '\r' < $G/src/classes/VBMicrolensingLibrary.cpp > VBMicrolensingLibrary.cpp
tr -d '\r' < $G/src/headers/VBMicrolensingLibrary.h  > VBMicrolensingLibrary.h
git commit -qam ours
git -c merge.conflictStyle=diff3 merge -Xdiff-algorithm=histogram -Xignore-space-change theirs
```

Then:

- Resolve conflicts. Take upstream's logic and re-apply the customizations above. For
  heavily rewritten functions (for example, `cmplx_roots_multigen` and `ComputeParallax`
  in 5.6), the cleanest approach is to take the upstream function whole and re-add the
  try/budget wrapper and loop checks.
- **Check for silently misplaced wrappers.** In 5.6, the auto-merge moved
  `TripleAstroLightCurve`'s `catch` block onto the new `TripleAstroLightCurveOrbital`.
  This prints every `HandleTimeoutError(err, "X")` that is not inside function `X`; it
  should print nothing once conflicts are resolved. Upstream renames show up here too:
  ```bash
  awk '/^[A-Za-z_].*VBMicrolensing::[~A-Za-z_0-9]+\(/{match($0,/VBMicrolensing::[~A-Za-z_0-9]+/); fn=substr($0,RSTART+16,RLENGTH-16)}
       /HandleTimeoutError\(err, "/{match($0,/"[^"]+"/); n=substr($0,RSTART+1,RLENGTH-2); if (n!=fn) print FILENAME ":" NR ": " n " is inside " fn}' VBMicrolensingLibrary.cpp
  ```
- **List upstream's new public functions** and add timeout wrappers to any that can run
  for a long time:
  ```bash
  pub() { git -C ../up show $1:VBMicrolensing/lib/VBMicrolensingLibrary.h | tr -d '\r' \
    | grep -oE '^\s+[A-Za-z_][A-Za-z_0-9<>\*& ]*\s\**([A-Za-z_][A-Za-z_0-9]*)\(' \
    | grep -oE '[A-Za-z_][A-Za-z_0-9]*\($' | tr -d '(' | sort -u; }
  comm -13 <(pub $OLD) <(pub $NEW)
  ```
- Audit with `diff -w <upstream NEW normalized> <merged>`. Every remaining difference
  should map to the customization list above.
- Convert back to CRLF before copying into gulls (`perl -pe 's/\r?\n/\r\n/'`). That
  keeps the line endings of the existing file and upstream.
- Update the version line at the top of both files and the table above.

### Dry run against upstream `dev` (2026-11-18)

Running the script above with `OLD=v5.6.1` and `NEW=origin/dev` (upstream `4b452c2`,
about 190 changed library lines) auto-merged everything except four small conflicts. All
four are upstream renames next to gulls timeout lines: `LoadSunTable` →
`LoadRefEphemeris` and `t_in_HJD` → `t_in_BJD` (three in `ComputeParallax`). Each
resolves as "keep the timeout lines, take the rename". The misplaced-wrapper check then
flags the `LoadSunTable` catch label inside `LoadRefEphemeris`, and the new-function
listing shows `LoadRefEphemeris`, `SatellitePosition` and `SetRefEphemerisFile`. With
`NEW=main` (no library changes since `v5.6.1`), the merge is a no-op and reproduces the
gulls files exactly.

## Validation

- Build: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build`
  (on macOS with Homebrew gfortran, `export SDKROOT=$(xcrun --show-sdk-path)` first).
- Run the smoke tests for the old and new builds and compare outputs. `true_relative_flux`
  should agree to about VBM tolerance (≲1e-6 relative). `measured_*` columns can differ
  in short windows: photometric noise uses a rejection-sampled `poisson()` on a shared
  RNG stream, so tiny changes in magnification shift some draws. Fisher errors on poorly
  constrained parameters can move by several percent because they use
  finite-difference derivatives.
- Exercise the timeout paths, because the smoke tests never hit them. Build a small
  program against the library with AddressSanitizer
  (`c++ -std=c++17 -O1 -g -fsanitize=address,undefined -I src/headers t.cpp src/classes/VBMicrolensingLibrary.cpp`)
  that calls `BinaryMag2`, `MultiMag2` and `ESPLMag2` (gulls' per-epoch calls) under
  `ErrorPolicy::ReturnNaN`. Sweep budgets from about 1e-7 s up, with `check_interval` 1–64
  and with `magnification`/`root_solver` enabled together and separately. Each call should
  either match the untimed result or set `HasLastError()` and return NaN, with
  `GetLastError().api` naming the called function and no sanitizer reports. Run the
  sweep a few times, because the lost-error case it guards against is timing-dependent.
  This is how the `annulus` and `Images == nullptr` crashes above were found.
