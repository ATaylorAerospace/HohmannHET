# 🚀 hohmann-het (Python)

Python implementation of **HohmannHET**, a tri-language library for low-thrust
orbital transfers combining Keplerian Hohmann mechanics, high fidelity Hall
Effect Thruster (HET) propulsion models, and mission optimization solvers.

Physical quantities use `astropy.units`. Results match the C++20 and MATLAB
ports to within a floating point tolerance of `1e-6`.

**Author:** A Taylor | **Reference:** Vallado / Curtis

## 📦 Modules

- **dynamics**: Keplerian Hohmann transfer, circular velocity, orbital period
- **propulsion**: HET operating point, Tsiolkovsky propellant budget, burn time
- **optimization**: Golden-section solver for minimum-propellant / minimum-TOF design

## 🚀 Install

```bash
pip install -e ".[dev]"
```

## 🧪 Test

```bash
pytest tests/ -v
```

## 🔧 Usage

```python
import astropy.units as u
from hohmann_het import compute_hohmann, compute_het_state, optimize_isp

transfer = compute_hohmann(400.0 * u.km, 35786.0 * u.km)
print(f"Total dv : {transfer.total_dv:.4f}")
print(f"TOF      : {transfer.tof_hours:.3f}")

het = compute_het_state(300.0, 1350.0, 0.50)
print(f"Isp      : {het.isp:.1f}")

result = optimize_isp(1000.0 * u.kg, transfer.total_dv, 5000.0 * u.W, 0.55)
print(f"Optimal Isp : {result.optimal_isp:.1f}")
print(f"Propellant  : {result.propellant_mass:.2f}")
```

See the repository root `README.md` for the full tri-language documentation,
C++ build instructions, and the MATLAB usage guide.
