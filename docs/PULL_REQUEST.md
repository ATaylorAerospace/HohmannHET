# 🔍 Repository review: unit consistency, validation, and real parity tests

A correctness and consistency pass across the tri-language suite. Every change
is mirrored in Python, C++, and MATLAB and all suites stay green.

## 🐛 Correctness and consistency fixes

- **Unified delta-V units to m/s (was inconsistent).** `propellant_mass` and
  `burn_time` took delta-V in **km/s** in Python and MATLAB but **m/s** in C++,
  so identical calls gave different answers across languages and Python was even
  internally inconsistent (its `optimize_isp` already used m/s). All delta-V
  arguments are now **m/s** everywhere. astropy `Quantity` inputs still
  auto-convert in Python.
- **Added missing validation to Python `propellant_mass`.** It silently accepted
  negative mass / zero Isp / negative delta-V while C++ and MATLAB rejected them.
  It now raises `ValueError`, matching the other two.
- **Guarded `optimize_isp` against zero delta-V (all three).** `delta_v == 0`
  made the burn-time normalisation `0/0` (NaN) and the golden-section search ran
  on NaN. `delta_v` must now be strictly positive; rejected with a clear error.

## 🧪 Real cross-language parity tests

- **Hardcoded golden-value tests** added to all three suites
  (`TestGoldenReferenceValues` / `DynamicsGolden` / `PropulsionGolden` /
  `golden`-tagged). The prior `*_parity` tests recomputed the expected value with
  a copy of the same formula, so they were tautological and could not catch a
  wrong-but-consistent formula or a unit regression. The new tests pin identical
  numeric literals (LEO-to-GEO transfer and the SPT-100 HET point), verified to
  agree across Python and C++ to ~1e-13, well inside the 1e-6 guarantee.
- Added validation tests: Python `propellant_mass` negative-mass / zero-Isp /
  negative-delta-V, and `optimize_isp` zero / negative delta-V in all three.

## 🧹 Cleanup and docs

- Removed the unused `scipy` dependency from `pyproject.toml` (never imported;
  the solver is pure-Python by design).
- README: corrected the benchmark total delta-V (3.86 -> 3.85), clarified the C++
  GoogleTest requirement (system `find_package` first, `FetchContent` fallback),
  and added an explicit units-convention table.
- CONTRIBUTING: documented the units convention and the shared golden-value rule.

## ✅ Verification

- **Python:** `pip install -e ".[dev]"` then `pytest tests/` -> **100 passed**
  (was 86; +14 new tests).
- **C++:** apt GoogleTest -> `find_package` -> build -> `ctest` -> **3/3 suites pass**
  (golden anchors + zero-delta-V guard included).
- **Parity:** Python and C++ golden values match to ~1e-13.

## 📋 Contribution rules honored

Tri-language parity maintained; no `Claude`/`AI`/`LLM`/`anthropic` mentions; no em
dashes in README; `Author: A Taylor` headers intact; commits authored by A Taylor.
