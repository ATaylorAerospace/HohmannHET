# Author: A Taylor | Purpose: pytest suite for optimization module | Ref: Vallado/Curtis
"""pytest test suite for hohmann_het.optimization.

Validates golden-section search correctness, optimization physics,
and the min_propellant_transfer convenience wrapper.
"""
import math
import pytest
import astropy.units as u

from hohmann_het.optimization import (
    golden_section_minimize,
    optimize_isp,
    min_propellant_transfer,
)

TOL = 1e-6

M0   = 1000.0   # kg
DV   = 3852.6   # m/s  (LEO-GEO Hohmann)
P    = 5000.0   # W
ETA  = 0.55


class TestGoldenSectionMinimize:
    def test_quadratic_minimum(self):
        """Find minimum of f(x) = (x - 3)^2 on [0, 10]."""
        result = golden_section_minimize(lambda x: (x - 3.0)**2, 0.0, 10.0)
        assert abs(result - 3.0) < 1e-6

    def test_offset_quadratic(self):
        result = golden_section_minimize(lambda x: (x - 1.5)**2, -5.0, 5.0)
        assert abs(result - 1.5) < 1e-6

    def test_at_boundary_low(self):
        """Monotone decreasing - minimum at left boundary."""
        result = golden_section_minimize(lambda x: x, 1.0, 10.0)
        assert result < 2.0  # converges towards lower bound

    def test_tighter_tolerance(self):
        result = golden_section_minimize(
            lambda x: (x - 2.718281828)**2,
            0.0, 5.0, tol=1e-12,
        )
        assert abs(result - 2.718281828) < 1e-9


class TestOptimizeIsp:
    def test_returns_result(self):
        r = optimize_isp(M0, DV, P, ETA)
        assert r is not None

    def test_optimal_isp_units(self):
        r = optimize_isp(M0, DV, P, ETA)
        assert r.optimal_isp.unit == u.s

    def test_propellant_mass_units(self):
        r = optimize_isp(M0, DV, P, ETA)
        assert r.propellant_mass.unit == u.kg

    def test_burn_time_units(self):
        r = optimize_isp(M0, DV, P, ETA)
        assert r.burn_time.unit == u.s

    def test_isp_within_bounds(self):
        r = optimize_isp(M0, DV, P, ETA, isp_min=500.0, isp_max=5000.0)
        assert 500.0 <= r.optimal_isp.value <= 5000.0

    def test_propellant_mass_positive(self):
        r = optimize_isp(M0, DV, P, ETA)
        assert r.propellant_mass.value > 0.0

    def test_propellant_mass_less_than_initial(self):
        r = optimize_isp(M0, DV, P, ETA)
        assert r.propellant_mass.value < M0

    def test_burn_time_positive(self):
        r = optimize_isp(M0, DV, P, ETA)
        assert r.burn_time.value > 0.0

    def test_objective_value_finite(self):
        r = optimize_isp(M0, DV, P, ETA)
        assert math.isfinite(r.objective_value)

    def test_higher_lambda_lowers_optimal_isp(self):
        """Heavier burn-time penalty pushes optimizer toward lower Isp."""
        r_low  = optimize_isp(M0, DV, P, ETA, lambda_weight=1e-6)
        r_high = optimize_isp(M0, DV, P, ETA, lambda_weight=1e-1)
        assert r_high.optimal_isp.value <= r_low.optimal_isp.value

    def test_accepts_quantity_inputs(self):
        r = optimize_isp(
            M0 * u.kg,
            DV * (u.m / u.s),
            P * u.W,
            ETA,
        )
        assert r.optimal_isp.value > 0.0

    def test_invalid_negative_mass(self):
        with pytest.raises(ValueError):
            optimize_isp(-100.0, DV, P, ETA)

    def test_invalid_efficiency(self):
        with pytest.raises(ValueError):
            optimize_isp(M0, DV, P, 1.5)

    def test_invalid_isp_bounds(self):
        with pytest.raises(ValueError):
            optimize_isp(M0, DV, P, ETA, isp_min=3000.0, isp_max=1000.0)


class TestMinPropellantTransfer:
    def test_end_to_end_leo_geo(self):
        """Convenience wrapper produces positive propellant mass for LEO-GEO."""
        r = min_propellant_transfer(
            h1=400.0, h2=35786.0,
            m_initial=M0, discharge_power=P,
            anode_efficiency=ETA,
        )
        assert r.propellant_mass.value > 0.0
        assert r.optimal_isp.value > 0.0

    def test_consistency_with_optimize_isp(self):
        """Wrapper result matches direct optimize_isp call."""
        from hohmann_het.dynamics import compute_hohmann
        transfer = compute_hohmann(400.0, 35786.0)

        r_direct = optimize_isp(M0, transfer.total_dv, P, ETA)
        r_wrap   = min_propellant_transfer(400.0, 35786.0, M0, P, ETA)

        assert abs(r_direct.optimal_isp.value - r_wrap.optimal_isp.value) < TOL
        assert abs(r_direct.propellant_mass.value - r_wrap.propellant_mass.value) < TOL
