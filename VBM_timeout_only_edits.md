# Minimal Timeout Integration Plan

**Goal:** Add optional timeouts to the `VBMicrolensingLibrary.cpp` monolith with minimal edits. Timeouts must be **opt‑in** and **disabled by default**, preserving legacy behavior unless explicitly configured.

---

## Principles

- **Minimal invasive changes**: introduce small, localized additions.
- **No behavior change by default**: timeouts disabled unless configured.
- **No silent fallbacks**: when a timeout triggers, throw a clear exception (or return a documented error code if exceptions are not acceptable in this repo).
- **Single-file edits preferred**: keep changes inside `VBMicrolensingLibrary.h/.cpp` unless absolutely necessary.

---

## Step 0 — Decide Timeout Error Mechanism

Choose one of the following (pick **one** based on the repo’s exception policy):

1. **Exceptions (preferred)**
   - Add a simple exception type `VBMTimeoutError : public std::runtime_error`.
   - Throw it on timeout.

2. **Error code return (fallback)**
   - Return `-1` or a sentinel value and set a new public flag (e.g., `VBM.timeout_triggered`).
   - Document the new flag and error value clearly.

If exceptions are already used in the monolith, use option 1.

---

## Step 1 — Add Minimal Timeout Config to `VBMicrolensingLibrary.h`

Add a small struct and setters. Keep defaults at **0 (disabled)**.

**Add to class `VBMicrolensing` public section:**

```cpp
struct TimeoutConfig {
  double root_solver_seconds = 0.0;
  double magnification_seconds = 0.0;
  double parallax_seconds = 0.0;
  double critical_curves_seconds = 0.0;
  double astrometry_seconds = 0.0;
  int check_interval = 0; // 0 = default interval
};

void SetTimeouts(const TimeoutConfig& cfg) { timeout_config_ = cfg; }
TimeoutConfig GetTimeouts() const { return timeout_config_; }
```

**Add to private section:**

```cpp
TimeoutConfig timeout_config_;
```

---

## Step 2 — Add Minimal Timeout Utility Helpers in `VBMicrolensingLibrary.cpp`

Add a tiny helper class and per-thread budget pointer (thread_local). Keep it **self‑contained**.

```cpp
namespace {
constexpr int kDefaultTimeoutCheckInterval = 256;

class TimeBudget {
 public:
  TimeBudget() : enabled_(false), budget_(0.0), start_(std::chrono::steady_clock::now()) {}
  static TimeBudget Disabled() { return TimeBudget(); }
  static TimeBudget FromSeconds(double seconds) {
    if (seconds <= 0.0) return Disabled();
    return TimeBudget(seconds);
  }
  bool enabled() const { return enabled_; }
  bool expired() const {
    if (!enabled_) return false;
    const auto elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(
      std::chrono::steady_clock::now() - start_).count();
    return elapsed > budget_;
  }

 private:
  explicit TimeBudget(double seconds)
      : enabled_(true), budget_(seconds), start_(std::chrono::steady_clock::now()) {}
  bool enabled_;
  double budget_;
  std::chrono::steady_clock::time_point start_;
};

thread_local const TimeBudget* g_budget = nullptr;
thread_local int g_check_interval = 0;

struct ScopedBudget {
  const TimeBudget* prev_budget;
  int prev_interval;
  ScopedBudget(const TimeBudget* b, int interval)
      : prev_budget(g_budget), prev_interval(g_check_interval) {
    g_budget = b;
    g_check_interval = interval;
  }
  ~ScopedBudget() { g_budget = prev_budget; g_check_interval = prev_interval; }
};

inline bool ShouldCheck(int iter, int fallback_interval) {
  int interval = g_check_interval > 0 ? g_check_interval : fallback_interval;
  if (interval <= 0) return true;
  return (iter % interval) == 0;
}

inline void CheckTimeout(const char* where) {
  if (!g_budget || !g_budget->enabled()) return;
  if (!g_budget->expired()) return;
  throw std::runtime_error(std::string("VBMTimeoutError: ") + (where ? where : ""));
}
} // namespace
```

**Note:** include `<chrono>` and `<stdexcept>` at top of `.cpp`.

---

## Step 3 — Wire Budgets at Public Entry Points (Minimal Insertions)

Add a scoped budget at the **top** of public, heavy-entry functions. Do **not** change signatures.

### Magnification entry points
Add this at top of each public magnification API (e.g., `PSPLMag`, `ESPLMag*`, `BinaryMag*`, `MultiMag*`, `BinaryLightCurve*`):

```cpp
TimeBudget budget = TimeBudget::FromSeconds(timeout_config_.magnification_seconds);
ScopedBudget scoped(&budget, timeout_config_.check_interval);
```

### Root solver
Wrap high‑level root solver entry points (e.g., `cmplx_roots_gen`, `cmplx_roots_multigen`) with:

```cpp
TimeBudget budget = TimeBudget::FromSeconds(timeout_config_.root_solver_seconds);
ScopedBudget scoped(&budget, timeout_config_.check_interval);
```

### Critical curves
At the top of `PlotCrit()` and `PlotCrit(a,q)`:

```cpp
TimeBudget budget = TimeBudget::FromSeconds(timeout_config_.critical_curves_seconds);
ScopedBudget scoped(&budget, timeout_config_.check_interval);
```

### Parallax / Astrometry
At the top of key public entry points (`SetObjectCoordinates`, `LoadSunTable`, and any astrometry light‑curve/centroid API):

```cpp
TimeBudget budget = TimeBudget::FromSeconds(timeout_config_.parallax_seconds);
ScopedBudget scoped(&budget, timeout_config_.check_interval);
```

### Loop checks
In long loops inside those functions (or their helpers), add:

```cpp
if (ShouldCheck(iter, kDefaultTimeoutCheckInterval)) {
  CheckTimeout("FunctionName");
}
```

Use an existing loop counter where possible to avoid new variables.

---

## Step 4 — Keep Defaults Non‑Breaking

- All timeout values default to **0.0** => timeouts disabled.
- The only new behavior occurs when the user calls `SetTimeouts()` with non‑zero values.

---

## Step 5 — Minimal Documentation (inline comments + README note)

Add a short block comment near the new `TimeoutConfig` explaining:

- `0.0` = disabled
- `check_interval` controls how often loops check timeouts
- Exceptions thrown on timeout

If README exists, add a **3‑line** note pointing users to the new timeout config.

---

## Step 6 — Manual Verification (No Tests Required)

Since no harness is available, do a manual run:

1. Build the library.
2. Create a small sample program:
   - Set `critical_curves_seconds = 0.01`.
   - Call `PlotCrit()` with large `NPcrit`.
   - Confirm timeout exception is raised.
3. Repeat with all timeouts at `0.0` to ensure no behavior change.

---

## Risks / Notes

- **Thread safety**: `thread_local` budgets keep configuration per‑thread; `SetTimeouts` is per‑instance (good). If this repo is heavily multithreaded, warn users to set timeouts **before** entering threaded work.
- **Exception policy**: if exceptions are not allowed, replace throws with error flags and document them.
- **Performance**: `check_interval` avoids overhead; default 256.

---

## Summary of Files to Edit

- `VBMicrolensingLibrary.h`
  - Add `TimeoutConfig`, `SetTimeouts`, `GetTimeouts`, `timeout_config_`.
- `VBMicrolensingLibrary.cpp`
  - Add `TimeBudget`, `ScopedBudget`, `CheckTimeout` helpers.
  - Add scoped budgets in entry points.
  - Add loop timeout checks in heavy loops.
- Optional: `README` for minimal note.

---

## Minimal Diff Strategy

- Keep new code in **one helper block** near top of `.cpp`.
- Only add **1–2 lines** per public entry point to enable the scoped budget.
- Insert loop checks only in **clearly long loops** (annuli loops, critical curve loops, root‑solver iterations).

This yields a low‑risk, low‑diff timeout integration with no default behavior changes.

# Minimal Timeout Integration Implementation

## 1. Scope Actually Implemented

This implementation ended up including four related pieces of work:

1. VBM timeout support inside `VBMicrolensingLibrary` (opt-in, default off).
2. Parameter-file wiring for all VBM timeout knobs in gulls.
3. Explicit timeout classification + timeout-specific process exit codes.
4. Two adjacent production-safety fixes discovered while validating:
   - image filename construction bug producing literal `%s/%d` paths.
   - hard guard against accidental use of `random.cpp` fallback/stub in science runs.

The sections below document each code change and the validation evidence.

---

## 2. VBM Timeout Core Implementation

### 2.1 `VBMicrolensingLibrary.h` API and error type

`/Users/malpas.1/Code/gulls_general/src/headers/VBMicrolensingLibrary.h`

- Added `TimeoutConfig` to `VBMicrolensing` public API:
  - `root_solver_seconds`
  - `magnification_seconds`
  - `parallax_seconds`
  - `critical_curves_seconds`
  - `astrometry_seconds`
  - `check_interval`
- Added setters/getters:
  - `SetTimeouts(const TimeoutConfig&)`
  - `GetTimeouts() const`
- Added private storage:
  - `timeout_config_`

Locations:
- `TimeoutConfig`: line ~198
- Set/Get: lines ~207-208
- private storage: line ~359

Also upgraded timeout exception metadata:

- `class VBMTimeoutError` now includes:
  - `enum class TimeoutCategory { Unknown, RootSolver, Magnification, Parallax, CriticalCurves, Astrometry }`
  - `category()` accessor
  - `where()` accessor
  - static `CategoryName(...)`
  - structured message format: `VBMTimeoutError[Category]: function_name`

Locations:
- `VBMTimeoutError` block: lines ~69-115

Rationale:
- This avoids brittle message parsing and gives deterministic machine-readable category data.

### 2.2 `VBMicrolensingLibrary.cpp` timeout machinery

`/Users/malpas.1/Code/gulls_general/src/classes/VBMicrolensingLibrary.cpp`

Added/used:

- `TimeBudget` with `enabled()/expired()` and `FromSeconds(...)`.
- `thread_local` state:
  - current budget pointer (`g_budget`)
  - check interval (`g_check_interval`)
  - timeout category (`g_timeout_category`)
- `ScopedBudget` now carries and restores all three.
- Timeout check helper `CheckTimeout(...)` now throws:
  - `VBMTimeoutError(g_timeout_category, where)`

Key lines:
- helper block starts ~46
- category thread-local and scoped propagation: lines ~81-107
- timeout throw site: line ~141

### 2.3 Entry-point category wiring

All existing timeout-scoped entry points were updated to pass explicit category into `ScopedBudget(...)`.

Categories now attached at entry:

- `Magnification` (many `BinaryMag*`, `MultiMag*`, etc.) examples:
  - lines ~592, ~605, ~658, ... (multiple locations)
- `Astrometry`:
  - lines ~5091, ~5128, ~5165, ...
- `Parallax`:
  - lines ~6581, ~6710, ~6876, ~6928
- `CriticalCurves`:
  - lines ~7268, ~7363
- `RootSolver`:
  - lines ~8080, ~8194

Result:
- timeout throws now include both subsystem category and function source.

---

## 3. Parameter Wiring in gulls

### 3.1 New `filekeywords` fields

`/Users/malpas.1/Code/gulls_general/src/structures.h`

Added:

- `vbm_timeout_root_solver`
- `vbm_timeout_magnification`
- `vbm_timeout_parallax`
- `vbm_timeout_critical_curves`
- `vbm_timeout_astrometry`
- `vbm_timeout_check_interval`
- `allow_random_stub` (for random fallback protection)

Location:
- lines ~171-178

### 3.2 Parameter parser defaults and assignments

`/Users/malpas.1/Code/gulls_general/src/readParamfile.cpp`

New defaults:

- `VBM_TIMEOUT_ROOT_SOLVER=0.0`
- `VBM_TIMEOUT_MAGNIFICATION=0.0`
- `VBM_TIMEOUT_PARALLAX=0.0`
- `VBM_TIMEOUT_CRITICAL_CURVES=0.0`
- `VBM_TIMEOUT_ASTROMETRY=0.0`
- `VBM_TIMEOUT_CHECK_INTERVAL=0`
- `ALLOW_RANDOM_STUB=0`

Location:
- lines ~66-73

New assignments:

- parse all six VBM timeout fields into `Paramfile`
- parse `ALLOW_RANDOM_STUB`

Location:
- lines ~240-247

### 3.3 Applying timeout config to VBM instance

`/Users/malpas.1/Code/gulls_general/src/gulls.cpp`

`VBM.SetTimeouts(...)` now uses dedicated VBM parameters (not `LC_TIMEOUT`):

- `timeout_cfg.root_solver_seconds = Paramfile.vbm_timeout_root_solver`
- `timeout_cfg.magnification_seconds = Paramfile.vbm_timeout_magnification`
- `timeout_cfg.parallax_seconds = Paramfile.vbm_timeout_parallax`
- `timeout_cfg.critical_curves_seconds = Paramfile.vbm_timeout_critical_curves`
- `timeout_cfg.astrometry_seconds = Paramfile.vbm_timeout_astrometry`
- `timeout_cfg.check_interval = Paramfile.vbm_timeout_check_interval`

Location:
- lines ~232-239

This leaves `LC_TIMEOUT` in place for the existing outer lightcurve-generator watchdog behavior.

---

## 4. Runtime Timeout Behavior in Science Runs (Current Semantics)

### 4.1 Top-level catch and exit code mapping

`/Users/malpas.1/Code/gulls_general/src/gulls.cpp`

Added timeout exit code mapping:

- `40` = `Unknown`
- `41` = `RootSolver`
- `42` = `Magnification`
- `43` = `Parallax`
- `44` = `CriticalCurves`
- `45` = `Astrometry`

Mapping function:
- lines ~62-78

Top-level timeout catch in `main`:
- logs:
  - `FATAL: VBMTimeoutError[...]`
  - `Timeout category: ...`
  - `Timeout source: ...`
- writes same to logfile if open
- returns mapped exit code

Location:
- catch block ~507-523

### 4.2 Exactly what happens when VBM timeout triggers

Current implemented behavior is **process-level fail-fast**, not per-event skip:

1. VBM call exceeds configured budget.
2. `VBMTimeoutError` is thrown.
3. Error propagates up to `main` catch in `gulls.cpp`.
4. Process exits immediately with timeout-specific exit code (`41`-`45`).
5. Remaining events in the subrun are **not processed**.

Important distinction:

- This is different from the existing `LC_TIMEOUT` watchdog path (`LCGEN_TIMEOUT_ERR`), which can mark an event as timed out and continue to next event.
- VBM timeout currently does **not** downgrade to event-local `lcerror`; it aborts the run.

If your desired science behavior is "skip only the bad event and continue":

- you would need to catch `VBMTimeoutError` inside per-event lightcurve generation logic, convert it to event-level timeout state, and continue looping.
- that behavior is not part of this implementation.

---

## 5. Smoke Parameter Files Updated

All runnable smoke `.prm` files under:

`/Users/malpas.1/Code/gulls_general/smoke_test/parameterfiles/`

were normalized to include:

- `ALLOW_RANDOM_STUB=1` (for intentional CI/smoke usage of fallback random backend)
- full `VBM_TIMEOUT_*` block, with each timeout value mirrored from that file’s `LC_TIMEOUT`
- `VBM_TIMEOUT_CHECK_INTERVAL=256`

Files updated:

- `smoke_croin.prm`
- `smoke_croin_binary.prm`
- `smoke_croin_heavy.prm`
- `smoke_fish.prm`
- `smoke_fish_binary.prm`
- `smoke_fish_heavy.prm`
- `smoke_general.prm`
- `smoke_general_binary.prm`
- `smoke_std.prm`
- `smoke_std_binary.prm`
- `smoke_std_heavy.prm`
- `smoke_std_houston_seed1.prm`
- `smoke_std_houston_seed2.prm`
- `smoke_std_houston_seed3.prm`

The placeholder `smoke.prm` was intentionally not modified (non-runnable doc placeholder).

---

## 6. Additional Miscellaneous Fixes (During Validation)

### 6.1 FITS filename construction bug (`%s/%d` literals)

`/Users/malpas.1/Code/gulls_general/src/outputLightcurve.cpp`

Issue:

- Field-selected branch used a comma expression:
  - `tmp1 = "%s%s_%d_%d_%d", ...`
- This resulted in malformed literal `%s/%d` filenames in repo root.
- There was also a second path concatenation issue:
  - names were assembled as `basefname + tmp1 + ...` where `tmp1` was being mutated in loop.

Fix:

- Corrected `tmp1` construction to plain C++ string concatenation.
- Introduced `obs_prefix = basefname + "." + obsidx + "_"`.
- Built FITS names as `obs_prefix + imtype + extension + ".fits"`.

Locations:

- corrected branch: lines ~540-543
- loop prefix and filename assembly: lines ~550, ~556, ~576

### 6.2 Random stub guard to protect science runs

Files:

- `/Users/malpas.1/Code/gulls_general/src/headers/random.h`
- `/Users/malpas.1/Code/gulls_general/src/classes/random.cpp`
- `/Users/malpas.1/Code/gulls_general/src/gulls.cpp`

Added backend metadata API:

- `bool gulls_random_is_stub();`
- `const char* gulls_random_backend_name();`

Stub implementation (`src/classes/random.cpp`) returns:

- `gulls_random_is_stub() == true`
- backend name: `"gsl_fallback_stub"`

`gulls.cpp` guard:

- after parsing params, if stub backend and `ALLOW_RANDOM_STUB != 1`, abort immediately with clear fatal message.

Current safety outcome:

- accidental science runs with fallback stub now fail hard unless explicitly overridden in `.prm`.

---

## 7. Validation Performed (Detailed)

### 7.1 Build validation

Repeated full builds succeeded after each major change set:

```bash
cmake --build build -j8
```

No compile/link failures introduced by timeout, exit code, parser, filename, or random-guard changes.

### 7.2 Smoke validation for general executable

Primary runtime command used throughout:

```bash
GULLS_BASE_DIR=/Users/malpas.1/Code/gulls_general/ \
GULLS_STARS_DIR=/Users/malpas.1/Code/gulls_general/ \
./bin/gulls_general.x -i smoke_test/parameterfiles/smoke_general.prm -s 0 -f 0
```

Observed successful full completion multiple times (program timings printed, clean end).

### 7.3 Timeout trigger sweep (one knob at a time)

Method:

- For each timeout parameter in `smoke_general.prm`:
  - set that one to `1e-9`
  - set others to `0.0`
  - set `VBM_TIMEOUT_CHECK_INTERVAL=1`
  - run `gulls_general`

Results (saved in `/tmp/vbm_timeout_sweep_results.tsv`):

- `VBM_TIMEOUT_ROOT_SOLVER`: triggered (`exit 134` before top-level timeout catch existed)
- `VBM_TIMEOUT_MAGNIFICATION`: triggered (`exit 134` before top-level timeout catch existed)
- `VBM_TIMEOUT_PARALLAX`: not triggered in this smoke scenario
- `VBM_TIMEOUT_CRITICAL_CURVES`: not triggered in this smoke scenario
- `VBM_TIMEOUT_ASTROMETRY`: not triggered in this smoke scenario

Interpretation:

- root/magnification timeout machinery was proven active.
- parallax/critical/astrometry code paths are not heavily exercised by `smoke_general` configuration.

### 7.4 Post-categorization exit code tests

After adding structured categories + timeout exit code mapping:

- forced root timeout: `exit=41`
  - log shows:
    - `FATAL: VBMTimeoutError[RootSolver]: cmplx_roots_gen`
    - `Timeout category: RootSolver`
    - `Timeout source: cmplx_roots_gen`
- forced magnification timeout: `exit=42`
  - log shows:
    - `FATAL: VBMTimeoutError[Magnification]: cmplx_roots_gen`
    - `Timeout category: Magnification`
    - `Timeout source: cmplx_roots_gen`

These checks confirm deterministic category propagation and exit code mapping.

### 7.5 Random stub guard validation

Fail case:

- temporary `.prm` with `ALLOW_RANDOM_STUB=0`
- observed:
  - `exit_code=1`
  - fatal guard messages:
    - backend is CI/testing stub
    - set `ALLOW_RANDOM_STUB=1` to override
    - aborting to protect science runs

Pass case:

- standard smoke `.prm` with `ALLOW_RANDOM_STUB=1`
- observed:
  - normal full run completion
  - program timings and execution-end stamp present

### 7.6 FITS filename bug validation

After filename fix:

- no new root-level weird files matching `%s%s_*`
- no FITS path/create errors in smoke run log
- generated FITS names are valid and located under expected output directories (e.g. `smoke_test/output/general/smoke_general/...`).

---

## 8. Behavior Summary for Review / PR Discussion

### What is guaranteed now

- VBM timeout controls are configurable from `.prm`.
- Timeout checks are category-aware and produce explicit diagnostics.
- Timeout-triggered process exits are category-specific (`40`-`45`).
- Stub RNG backend cannot be used accidentally in science runs unless explicitly overridden.

### What is **not** implemented

- Per-event continuation after VBM timeout.

Current behavior on VBM timeout is **fail run immediately**.

If the desired behavior is "keep run alive and skip only timed-out events", that requires a follow-up design change in event-level lightcurve generation error handling.

---

## 9. Notes for Maintainers

- `ALLOW_RANDOM_STUB` default is intentionally `0` (safe-by-default for production/science).
- Smoke files set `ALLOW_RANDOM_STUB=1` intentionally to preserve CI/smoke behavior with fallback RNG.
- For external wrappers/orchestration, timeout cause can be inferred from process exit code (`41`-`45`) and logfile content.
