# Author: A Taylor | Purpose: pytest suite for dynamics module | Ref: Vallado/Curtis
"""pytest test suite for hohmann_het.dynamics.

Validates Keplerian dynamics against known LEO-to-GEO benchmark values
(Vallado, Table 6-1) and verifies cross-language parity constants.
Precision requirement: all core values within 1e-6 of reference.
"""
import math
import pytest
import astropy.units as u

from hohmann_het.dynamics import (
    compute_hohmann,
    circular_velocity,
    orbital_period,
    MU,
    R_EARTH,
)

# ---------------------------------------------------------------------------
# Reference constants (same across all three language implementations)
# ---------------------------------------------------------------------------
MU_REF     = 398600.4418   # km^3/s^2
R_EARTH_REF = 6378.137     # km
TOL        = 1e-6          # cross-language parity tolerance

# LEO-to-GEO benchmark
H1 = 400.0    # km
H2 = 35786.0  # km
R1_REF = R_EARTH_REF + H1   # 6778.137 km
R2_REF = R_EARTH_REF + H2   # 42164.137 km
A_REF  = (R1_REF + R2_REF) / 2.0  # 24471.137 km


# ---------------------------------------------------------------------------
# Helper: compute reference values analytically (same formula as source)
# ---------------------------------------------------------------------------
def _ref_hohmann(h1, h2):
    mu = MU_REF
    r1 = R_EARTH_REF + h1
    r2 = R_EARTH_REF + h2
    a  = (r1 + r2) / 2.0
    mua = mu / a
    vc1 = math.sqrt(mu / r1)
    vc2 = math.sqrt(mu / r2)
    vt1 = math.sqrt(2.0 * mu / r1 - mua)
    vt2 = math.sqrt(2.0 * mu / r2 - mua)
    dv1 = vt1 - vc1
    dv2 = vc2 - vt2
    tof = math.pi * a * math.sqrt(a / mu)
    return dv1, dv2, dv1 + dv2, tof


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
class TestConstants:
    def test_mu_value(self):
        assert abs(MU.value - MU_REF) < TOL

    def test_r_earth_value(self):
        assert abs(R_EARTH.value - R_EARTH_REF) < TOL

    def test_mu_units(self):
        assert MU.unit == u.km**3 / u.s**2

    def test_r_earth_units(self):
        assert R_EARTH.unit == u.km


# ---------------------------------------------------------------------------
# compute_hohmann - valid inputs
# ---------------------------------------------------------------------------
class TestComputeHohmann:
    def test_radii_units(self):
        r = compute_hohmann(H1, H2)
        assert r.r1.unit == u.km
        assert r.r2.unit == u.km
        assert r.a_transfer.unit == u.km

    def test_dv_units(self):
        r = compute_hohmann(H1, H2)
        assert r.dv_departure.unit == u.km / u.s
        assert r.dv_arrival.unit == u.km / u.s

    def test_tof_units(self):
        r = compute_hohmann(H1, H2)
        assert r.tof.unit == u.s

    def test_radii_values(self):
        r = compute_hohmann(H1, H2)
        assert abs(r.r1.value - R1_REF) < TOL
        assert abs(r.r2.value - R2_REF) < TOL

    def test_semi_major_axis(self):
        r = compute_hohmann(H1, H2)
        assert abs(r.a_transfer.value - A_REF) < TOL

    def test_dv_departure_parity(self):
        """Delta-V departure matches reference formula to 1e-6 km/s."""
        ref_dv1, _, _, _ = _ref_hohmann(H1, H2)
        r = compute_hohmann(H1, H2)
        assert abs(r.dv_departure.value - ref_dv1) < TOL

    def test_dv_arrival_parity(self):
        """Delta-V arrival matches reference formula to 1e-6 km/s."""
        _, ref_dv2, _, _ = _ref_hohmann(H1, H2)
        r = compute_hohmann(H1, H2)
        assert abs(r.dv_arrival.value - ref_dv2) < TOL

    def test_total_dv_parity(self):
        """Total delta-V matches reference to 1e-6 km/s."""
        _, _, ref_total, _ = _ref_hohmann(H1, H2)
        r = compute_hohmann(H1, H2)
        assert abs(r.total_dv.value - ref_total) < TOL

    def test_tof_parity(self):
        """TOF matches reference formula to 1e-6 s."""
        _, _, _, ref_tof = _ref_hohmann(H1, H2)
        r = compute_hohmann(H1, H2)
        assert abs(r.tof.value - ref_tof) < TOL

    def test_known_leo_geo_dv_departure(self):
        """LEO-to-GEO departure burn approx 2.4 km/s (Vallado Table 6-1)."""
        r = compute_hohmann(H1, H2)
        assert 2.3 < r.dv_departure.value < 2.5

    def test_known_leo_geo_dv_arrival(self):
        """LEO-to-GEO arrival burn approx 1.5 km/s."""
        r = compute_hohmann(H1, H2)
        assert 1.3 < r.dv_arrival.value < 1.6

    def test_known_leo_geo_total_dv(self):
        """LEO-to-GEO total dv approx 3.9 km/s."""
        r = compute_hohmann(H1, H2)
        assert 3.7 < r.total_dv.value < 4.1

    def test_known_leo_geo_tof_hours(self):
        """LEO-to-GEO transfer time approx 5.3 hours."""
        r = compute_hohmann(H1, H2)
        assert 4.5 < r.tof_hours.value < 6.0

    def test_accepts_quantity_input(self):
        r = compute_hohmann(H1 * u.km, H2 * u.km)
        assert abs(r.r1.value - R1_REF) < TOL

    def test_same_altitude_zero_dv(self):
        """Same-altitude transfer requires zero delta-V."""
        r = compute_hohmann(H1, H1)
        assert abs(r.dv_departure.value) < TOL
        assert abs(r.dv_arrival.value) < TOL
        assert abs(r.total_dv.value) < TOL

    def test_reverse_transfer_positive_dv(self):
        """Reverse (GEO-to-LEO) transfer yields positive total delta-V."""
        r = compute_hohmann(H2, H1)
        assert r.total_dv.value > 0.0

    def test_reverse_symmetry_total_dv(self):
        """Forward and reverse transfers have equal total delta-V."""
        fwd = compute_hohmann(H1, H2)
        rev = compute_hohmann(H2, H1)
        assert abs(fwd.total_dv.value - rev.total_dv.value) < 1e-6

    def test_reverse_symmetry_tof(self):
        """Forward and reverse transfers have equal TOF."""
        fwd = compute_hohmann(H1, H2)
        rev = compute_hohmann(H2, H1)
        assert abs(fwd.tof.value - rev.tof.value) < TOL

    def test_semi_major_axis_between_radii(self):
        """Transfer semi-major axis lies between r1 and r2 for ascending transfer."""
        r = compute_hohmann(H1, H2)
        assert r.r1.value < r.a_transfer.value < r.r2.value

    def test_total_dv_equals_sum(self):
        """total_dv property equals |dv1| + |dv2|."""
        r = compute_hohmann(H1, H2)
        expected = abs(r.dv_departure.value) + abs(r.dv_arrival.value)
        assert abs(r.total_dv.value - expected) < TOL

    def test_tof_hours_conversion(self):
        r = compute_hohmann(H1, H2)
        assert abs(r.tof_hours.value - r.tof.value / 3600.0) < TOL


# ---------------------------------------------------------------------------
# compute_hohmann - invalid inputs
# ---------------------------------------------------------------------------
class TestComputeHohmannErrors:
    def test_negative_h1(self):
        with pytest.raises(ValueError):
            compute_hohmann(-100.0, H2)

    def test_negative_h2(self):
        with pytest.raises(ValueError):
            compute_hohmann(H1, -100.0)

    def test_zero_h1(self):
        with pytest.raises(ValueError):
            compute_hohmann(0.0, H2)

    def test_zero_h2(self):
        with pytest.raises(ValueError):
            compute_hohmann(H1, 0.0)


# ---------------------------------------------------------------------------
# circular_velocity
# ---------------------------------------------------------------------------
class TestCircularVelocity:
    def test_leo_value(self):
        vc = circular_velocity(R1_REF)
        assert abs(vc.value - math.sqrt(MU_REF / R1_REF)) < TOL

    def test_units(self):
        vc = circular_velocity(R1_REF)
        assert vc.unit == u.km / u.s

    def test_negative_radius(self):
        with pytest.raises(ValueError):
            circular_velocity(-1.0)


# ---------------------------------------------------------------------------
# orbital_period
# ---------------------------------------------------------------------------
class TestOrbitalPeriod:
    def test_geo_period_approx_24h(self):
        """GEO orbital period is approximately 24 hours."""
        T = orbital_period(R2_REF)
        assert abs(T.to(u.hour).value - 24.0) < 0.1

    def test_units(self):
        T = orbital_period(R1_REF)
        assert T.unit == u.s

    def test_formula_consistency(self):
        a = A_REF
        T = orbital_period(a)
        expected = 2.0 * math.pi * a * math.sqrt(a / MU_REF)
        assert abs(T.value - expected) < TOL
