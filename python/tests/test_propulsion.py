# Author: A Taylor | Purpose: pytest suite for propulsion module | Ref: Vallado/Curtis
"""pytest test suite for hohmann_het.propulsion.

Validates HET models against SPT-100 class operating points and
verifies cross-language parity for Tsiolkovsky propellant calculations.
"""
import math
import pytest
import astropy.units as u

from hohmann_het.propulsion import (
    compute_het_state,
    propellant_mass,
    burn_time,
    G0,
    E_CHARGE,
    M_XE,
)

TOL = 1e-6

# SPT-100 reference operating point
V_D_REF  = 300.0   # V
P_D_REF  = 1350.0  # W
ETA_REF  = 0.50    # anode efficiency

# Reference exhaust velocity (from physics)
VE_REF = math.sqrt(2.0 * ETA_REF * E_CHARGE * V_D_REF / M_XE)
ISP_REF = VE_REF / G0.value


class TestConstants:
    def test_g0_value(self):
        assert abs(G0.value - 9.80665) < TOL

    def test_e_charge(self):
        assert abs(E_CHARGE - 1.602176634e-19) < 1e-28

    def test_m_xe_order(self):
        assert 2.17e-25 < M_XE < 2.19e-25


class TestComputeHETState:
    def test_isp_units(self):
        s = compute_het_state(V_D_REF, P_D_REF, ETA_REF)
        assert s.isp.unit == u.s

    def test_exhaust_velocity_units(self):
        s = compute_het_state(V_D_REF, P_D_REF, ETA_REF)
        assert s.exhaust_velocity.unit == u.m / u.s

    def test_thrust_units(self):
        s = compute_het_state(V_D_REF, P_D_REF, ETA_REF)
        assert s.thrust.unit == u.N

    def test_mass_flow_units(self):
        s = compute_het_state(V_D_REF, P_D_REF, ETA_REF)
        assert s.mass_flow.unit == u.kg / u.s

    def test_isp_parity(self):
        """Isp matches beam-voltage formula to 1e-6 s."""
        s = compute_het_state(V_D_REF, P_D_REF, ETA_REF)
        assert abs(s.isp.value - ISP_REF) < TOL

    def test_exhaust_velocity_parity(self):
        """ve matches sqrt formula to 1e-6 m/s."""
        s = compute_het_state(V_D_REF, P_D_REF, ETA_REF)
        assert abs(s.exhaust_velocity.value - VE_REF) < TOL

    def test_isp_reasonable_range(self):
        """Isp for SPT-100 class should be 1000-2500 s."""
        s = compute_het_state(V_D_REF, P_D_REF, ETA_REF)
        assert 1000.0 < s.isp.value < 2500.0

    def test_thrust_power_relation(self):
        """Thrust = 2 * eta * P / ve."""
        s = compute_het_state(V_D_REF, P_D_REF, ETA_REF)
        expected_F = 2.0 * ETA_REF * P_D_REF / VE_REF
        assert abs(s.thrust.value - expected_F) < TOL

    def test_thrust_mass_flow_ratio(self):
        """ve = F / mdot."""
        s = compute_het_state(V_D_REF, P_D_REF, ETA_REF)
        ve_check = s.thrust.value / s.mass_flow.value
        assert abs(ve_check - s.exhaust_velocity.value) < TOL

    def test_higher_voltage_higher_isp(self):
        s1 = compute_het_state(300.0, P_D_REF, ETA_REF)
        s2 = compute_het_state(600.0, P_D_REF, ETA_REF)
        assert s2.isp.value > s1.isp.value

    def test_accepts_quantity_inputs(self):
        s = compute_het_state(V_D_REF * u.V, P_D_REF * u.W, ETA_REF)
        assert abs(s.isp.value - ISP_REF) < TOL

    def test_stored_efficiency(self):
        s = compute_het_state(V_D_REF, P_D_REF, ETA_REF)
        assert s.anode_efficiency == ETA_REF

    def test_stored_discharge_power(self):
        s = compute_het_state(V_D_REF, P_D_REF, ETA_REF)
        assert abs(s.discharge_power.value - P_D_REF) < TOL


class TestComputeHETStateErrors:
    def test_negative_voltage(self):
        with pytest.raises(ValueError):
            compute_het_state(-300.0, P_D_REF, ETA_REF)

    def test_zero_voltage(self):
        with pytest.raises(ValueError):
            compute_het_state(0.0, P_D_REF, ETA_REF)

    def test_negative_power(self):
        with pytest.raises(ValueError):
            compute_het_state(V_D_REF, -1350.0, ETA_REF)

    def test_efficiency_zero(self):
        with pytest.raises(ValueError):
            compute_het_state(V_D_REF, P_D_REF, 0.0)

    def test_efficiency_one(self):
        with pytest.raises(ValueError):
            compute_het_state(V_D_REF, P_D_REF, 1.0)

    def test_efficiency_above_one(self):
        with pytest.raises(ValueError):
            compute_het_state(V_D_REF, P_D_REF, 1.5)


class TestPropellantMass:
    # Reference: m0=1000 kg, dv=3.8526 km/s (LEO-GEO), Isp=ISP_REF
    M0 = 1000.0   # kg
    DV = 3.8526   # km/s

    def _mp_ref(self, isp):
        ve = isp * G0.value
        dv_ms = self.DV * 1000.0
        return self.M0 * (1.0 - math.exp(-dv_ms / ve))

    def test_units(self):
        mp = propellant_mass(self.M0, self.DV, ISP_REF)
        assert mp.unit == u.kg

    def test_parity_reference(self):
        """Propellant mass matches Tsiolkovsky formula to 1e-6 kg."""
        mp = propellant_mass(self.M0, self.DV, ISP_REF)
        assert abs(mp.value - self._mp_ref(ISP_REF)) < TOL

    def test_higher_isp_lower_prop_mass(self):
        mp1 = propellant_mass(self.M0, self.DV, 1500.0)
        mp2 = propellant_mass(self.M0, self.DV, 3000.0)
        assert mp2.value < mp1.value

    def test_zero_dv_zero_prop(self):
        mp = propellant_mass(self.M0, 0.0, ISP_REF)
        assert abs(mp.value) < TOL

    def test_positive_result(self):
        mp = propellant_mass(self.M0, self.DV, ISP_REF)
        assert mp.value > 0.0

    def test_less_than_initial_mass(self):
        mp = propellant_mass(self.M0, self.DV, ISP_REF)
        assert mp.value < self.M0


class TestBurnTime:
    def test_positive_result(self):
        s = compute_het_state(V_D_REF, P_D_REF, ETA_REF)
        tb = burn_time(
            s.thrust,
            s.mass_flow,
            3.8526 * (u.km / u.s),
            1000.0 * u.kg,
        )
        assert tb.value > 0.0

    def test_units(self):
        s = compute_het_state(V_D_REF, P_D_REF, ETA_REF)
        tb = burn_time(s.thrust, s.mass_flow, 3.8526, 1000.0)
        assert tb.unit == u.s

    def test_larger_dv_longer_burn(self):
        s = compute_het_state(V_D_REF, P_D_REF, ETA_REF)
        tb1 = burn_time(s.thrust, s.mass_flow, 1.0, 1000.0)
        tb2 = burn_time(s.thrust, s.mass_flow, 3.0, 1000.0)
        assert tb2.value > tb1.value
