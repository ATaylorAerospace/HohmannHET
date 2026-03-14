// Author: A Taylor | Purpose: GoogleTest suite for propulsion module | Ref: Vallado/Curtis
/**
 * @file test_propulsion.cpp
 * @brief GoogleTest suite for hohmann_het::propulsion.
 *
 * Validates HET models against SPT-100 class operating points and
 * verifies cross-language parity for Tsiolkovsky propellant calculations.
 * Precision requirement: all core values within 1e-6 of reference.
 *
 * @author A Taylor
 */

#include <gtest/gtest.h>
#include <cmath>
#include <stdexcept>

#include "hohmann_het/propulsion.hpp"

using namespace hohmann_het;
using namespace hohmann_het::constants;
using namespace hohmann_het::prop_constants;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------
static constexpr double TOL    = 1e-6;
static constexpr double V_D    = 300.0;    // V  (SPT-100 reference)
static constexpr double P_D    = 1350.0;   // W
static constexpr double ETA    = 0.50;     // anode efficiency
static constexpr double VE_REF =           // reference exhaust velocity [m/s]
    // computed at compile time would require constexpr sqrt - use inline
    0.0;  // filled in tests via std::sqrt

// ---------------------------------------------------------------------------
// Fixture
// ---------------------------------------------------------------------------
class PropulsionTest : public ::testing::Test {
protected:
    const double ve_ref  = std::sqrt(2.0 * ETA * E_CHARGE * V_D / M_XE);
    const double isp_ref = ve_ref / G0;
    const HETState state = compute_het_state(V_D, P_D, ETA);
};


// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------
TEST(PropulsionConstants, G0Value) {
    EXPECT_DOUBLE_EQ(G0, 9.80665);
}

TEST(PropulsionConstants, EChargeOrder) {
    EXPECT_GT(E_CHARGE, 1.60e-19);
    EXPECT_LT(E_CHARGE, 1.61e-19);
}

TEST(PropulsionConstants, MXeOrder) {
    EXPECT_GT(M_XE, 2.17e-25);
    EXPECT_LT(M_XE, 2.19e-25);
}


// ---------------------------------------------------------------------------
// compute_het_state
// ---------------------------------------------------------------------------
TEST_F(PropulsionTest, IspParity) {
    EXPECT_NEAR(state.isp, isp_ref, TOL);
}

TEST_F(PropulsionTest, ExhaustVelocityParity) {
    EXPECT_NEAR(state.exhaust_velocity, ve_ref, TOL);
}

TEST_F(PropulsionTest, IspReasonableRange) {
    // SPT-100 class: 1000-2500 s
    EXPECT_GT(state.isp, 1000.0);
    EXPECT_LT(state.isp, 2500.0);
}

TEST_F(PropulsionTest, ThrustPowerRelation) {
    const double expected_F = 2.0 * ETA * P_D / ve_ref;
    EXPECT_NEAR(state.thrust, expected_F, TOL);
}

TEST_F(PropulsionTest, ThrustMassFlowRatio) {
    // ve = F / mdot
    const double ve_check = state.thrust / state.mass_flow;
    EXPECT_NEAR(ve_check, state.exhaust_velocity, TOL);
}

TEST_F(PropulsionTest, HigherVoltageHigherIsp) {
    const auto s1 = compute_het_state(300.0, P_D, ETA);
    const auto s2 = compute_het_state(600.0, P_D, ETA);
    EXPECT_GT(s2.isp, s1.isp);
}

TEST_F(PropulsionTest, StoredEfficiency) {
    EXPECT_DOUBLE_EQ(state.anode_efficiency, ETA);
}

TEST_F(PropulsionTest, StoredDischargePower) {
    EXPECT_NEAR(state.discharge_power, P_D, TOL);
}


// ---------------------------------------------------------------------------
// compute_het_state - invalid inputs
// ---------------------------------------------------------------------------
TEST(PropulsionValidation, NegativeVoltageThrows) {
    EXPECT_THROW(compute_het_state(-300.0, P_D, ETA), std::invalid_argument);
}

TEST(PropulsionValidation, ZeroVoltageThrows) {
    EXPECT_THROW(compute_het_state(0.0, P_D, ETA), std::invalid_argument);
}

TEST(PropulsionValidation, NegativePowerThrows) {
    EXPECT_THROW(compute_het_state(V_D, -1350.0, ETA), std::invalid_argument);
}

TEST(PropulsionValidation, EfficiencyZeroThrows) {
    EXPECT_THROW(compute_het_state(V_D, P_D, 0.0), std::invalid_argument);
}

TEST(PropulsionValidation, EfficiencyOneThrows) {
    EXPECT_THROW(compute_het_state(V_D, P_D, 1.0), std::invalid_argument);
}

TEST(PropulsionValidation, EfficiencyAboveOneThrows) {
    EXPECT_THROW(compute_het_state(V_D, P_D, 1.5), std::invalid_argument);
}


// ---------------------------------------------------------------------------
// propellant_mass (Tsiolkovsky)
// ---------------------------------------------------------------------------
class TsiolkovskyTest : public ::testing::Test {
protected:
    static constexpr double M0     = 1000.0;   // kg
    static constexpr double DV_MS  = 3852.6;   // m/s  (LEO-GEO)
    const double isp_ref = std::sqrt(2.0 * ETA * E_CHARGE * V_D / M_XE) / G0;

    double mp_ref(double isp) const {
        const double ve = isp * G0;
        return M0 * (1.0 - std::exp(-DV_MS / ve));
    }
};

TEST_F(TsiolkovskyTest, PropellantMassParity) {
    const double mp = propellant_mass(M0, DV_MS, isp_ref);
    EXPECT_NEAR(mp, mp_ref(isp_ref), TOL);
}

TEST_F(TsiolkovskyTest, ZeroDvZeroProp) {
    EXPECT_NEAR(propellant_mass(M0, 0.0, isp_ref), 0.0, TOL);
}

TEST_F(TsiolkovskyTest, HigherIspLowerPropMass) {
    const double mp1 = propellant_mass(M0, DV_MS, 1500.0);
    const double mp2 = propellant_mass(M0, DV_MS, 3000.0);
    EXPECT_LT(mp2, mp1);
}

TEST_F(TsiolkovskyTest, PositiveResult) {
    EXPECT_GT(propellant_mass(M0, DV_MS, isp_ref), 0.0);
}

TEST_F(TsiolkovskyTest, LessThanInitialMass) {
    EXPECT_LT(propellant_mass(M0, DV_MS, isp_ref), M0);
}

TEST_F(TsiolkovskyTest, NegativeMassThrows) {
    EXPECT_THROW(propellant_mass(-100.0, DV_MS, isp_ref), std::invalid_argument);
}

TEST_F(TsiolkovskyTest, NegativeDvThrows) {
    EXPECT_THROW(propellant_mass(M0, -1.0, isp_ref), std::invalid_argument);
}

TEST_F(TsiolkovskyTest, ZeroIspThrows) {
    EXPECT_THROW(propellant_mass(M0, DV_MS, 0.0), std::invalid_argument);
}


// ---------------------------------------------------------------------------
// burn_time
// ---------------------------------------------------------------------------
TEST(BurnTimeTest, PositiveResult) {
    const auto s = compute_het_state(V_D, P_D, ETA);
    EXPECT_GT(burn_time(s.thrust, s.mass_flow, 3852.6, 1000.0), 0.0);
}

TEST(BurnTimeTest, LargerDvLongerBurn) {
    const auto s  = compute_het_state(V_D, P_D, ETA);
    const double tb1 = burn_time(s.thrust, s.mass_flow, 1000.0, 1000.0);
    const double tb2 = burn_time(s.thrust, s.mass_flow, 3000.0, 1000.0);
    EXPECT_GT(tb2, tb1);
}

TEST(BurnTimeTest, ZeroDvZeroBurnTime) {
    const auto s = compute_het_state(V_D, P_D, ETA);
    EXPECT_NEAR(burn_time(s.thrust, s.mass_flow, 0.0, 1000.0), 0.0, TOL);
}

TEST(BurnTimeTest, InvalidThrustThrows) {
    EXPECT_THROW(burn_time(-1.0, 1e-4, 3000.0, 1000.0), std::invalid_argument);
}
