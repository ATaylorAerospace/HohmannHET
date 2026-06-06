# 🚀 Optimizer efficiency + CI reliability

Performance improvements to the optimization solver plus the CI fixes that get
all jobs green, applied identically across Python, MATLAB, and C++ to preserve
the `1e-6` cross-language parity guarantee.

## ⚡ Solver efficiency

- **Golden-section caching** (`optimization.py` / `optimization.hpp` / `Optimization.m`): the two interior probe points and their objective values are cached between iterations, so each step makes **one** new objective evaluation instead of two: roughly **2x fewer** `exp()`-bearing evaluations. Interval-narrowing decisions are unchanged, so results are numerically identical.
- **C++ template solver**: replaced `std::function<double(double)>` with a template type parameter so the callable inlines (no type-erasure / heap allocation on the hot path); dropped the `<functional>` include. Compiles clean under `-Wall -Wextra -Wpedantic`.

## 🔧 Parity fix

- Python `optimize_isp` now treats a bare-float `delta_v` as **m/s**, matching the MATLAB and C++ APIs (previously km/s in Python only). Fixes a latent unit mismatch and aligns all three signatures.

## 🧪 CI reliability

- Added **`python/README.md`** so the Hatchling editable install resolves its `readme` field (the build previously failed with `Readme file does not exist`, which cascaded and cancelled the 3.10/3.12 matrix).
- **pip cache** on the Python jobs; **prebuilt GoogleTest** via `libgtest-dev`/`libgmock-dev` with a `find_package(GTest)` fallback to `FetchContent`, so GoogleTest is no longer downloaded/recompiled each run.
- Removed the fragile `cpp/build` cache (a known CMake stale-cache footgun) now that GoogleTest is prebuilt.

## ✅ Verification (exact CI steps)

- **Python:** `pip install -e ".[dev]"` succeeds; `pytest tests/` -> **86 passed**.
- **C++:** apt GoogleTest -> `find_package` -> configure/build -> `ctest` -> **3/3 suites pass**.
- **Cross-language parity:** Python vs C++ optimizer agree to `|isp| = 4e-10`, `|mp| = 2e-7` (well under `1e-6`).

## 📋 Contribution rules honored

Tri-language parity maintained; no `Claude`/`AI`/`LLM`/`anthropic` mentions; no em dashes in README; `Author: A Taylor` headers intact; commits authored by A Taylor.

## Notes

- The **Node.js 20 lines** in the run log are deprecation **warnings** (June 2026 runtime change), not failures: `checkout@v4`/`setup-python@v5` are current.
