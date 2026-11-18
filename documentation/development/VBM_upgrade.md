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
   See `VBM_timeout_only_edits.md`. The wrappers re-indent function bodies, so a plain
   3-way merge produces many whitespace-only conflicts (see below).
3. **RAII ownership** in `BinaryMag`, `BinaryMag0`, `MultiMag`, `MultiMag0` and the
   `*Dark` functions. These use `std::unique_ptr` (`images_owner`, `thetas_owner`,
   `stheta_owner`, `AnnulusChainPtr`) instead of raw `new`/`delete`. `Prov` lists are
   drained with `while (Prov->first) { …; Prov->drop(scan1); }`. Callers check
   `Images == nullptr` and return NaN after a timeout. On upstream's own `return -1`
   paths, `MultiMag` still hands back the image list (`images_owner.release()`), which
   preserves upstream's contract with `MultiMagSafe`.
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
NEW=v5.7       # target upstream tag
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
  Check that every `HandleTimeoutError(err, "X")` is inside function `X`.
- Audit with `diff -w <upstream NEW normalized> <merged>`. Every remaining difference
  should map to the customization list above.
- Convert back to CRLF before copying into gulls (`perl -pe 's/\r?\n/\r\n/'`). That
  keeps the line endings of the existing file and upstream.
- Update the version line at the top of both files and the table above.

## Validation

- Build: `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build`
  (on macOS with Homebrew gfortran, `export SDKROOT=$(xcrun --show-sdk-path)` first).
- Run the smoke tests for the old and new builds and compare outputs. `true_relative_flux`
  should agree to about VBM tolerance (≲1e-6 relative). `measured_*` columns can differ
  in short windows: photometric noise uses a rejection-sampled `poisson()` on a shared
  RNG stream, so tiny changes in magnification shift some draws. Fisher errors on poorly
  constrained parameters can move by several percent because they use
  finite-difference derivatives.
