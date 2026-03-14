// Author: A Taylor | Purpose: GoogleTest suite for dynamics module | Ref: Vallado/Curtis
/**
 * @file test_dynamics.cpp
 * @brief GoogleTest suite for hohmann_het::dynamics.
 *
 * Validates Keplerian dynamics against known LEO-to-GEO benchmark values
 * (Vallado, Table 6-1) and verifies cross-language parity constants.
 * Precision requirement: all core values within 1e-6 of reference.
 *
 * @author A Taylor
 */

#include <gtest/gtest.h>
#include <cmath>
#include <stdexcept>

#include "hohmann_het/dynamics.hpp"

using namespace hohmann_het;
using namespace hohmann_het::constants;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------
static constexpr double TOL         = 1e-6;
static constexpr double LEO_ALT     = 400.0;
static constexpr double GEO_ALT     = 35786.0;
static constexpr double R1_REF      = R_EARTH + LEO_ALT;    // 6778.137
static constexpr double R2_REF      = R_EARTH + GEO_ALT;    // 42164.137
static constexpr double A_REF       = (R1_REF + R2_REF) / 2.0;

// ---------------------------------------------------------------------------
// Reference implementation (same formula as source, for parity checks)
// ---------------------------------------------------------------------------
struct RefResult { double dv1, dv2, total, tof; };

static RefResult ref_hohmann(double h1, double h2)
{
    const double r1  = R_EARTH + h1;
    const double r2  = R_EARTH + h2;
    const double a   = (r1 + r2) * 0.5;
    const double mua = MU / a;
    const double vc1 = std::sqrt(MU / r1);
    const double vc2 = std::sqrt(MU / r2);
    const double vt1 = std::sqrt(2.0 * MU / r1 - mua);
    const double vt2 = std::sqrt(2.0 * MU / r2 - mua);
    return { vt1 - vc1, vc2 - vt2, (vt1 - vc1) + (vc2 - vt2),
             PI * a * std::sqrt(a / MU) };
}

// ---------------------------------------------------------------------------
// Fixture
// ---------------------------------------------------------------------------
class DynamicsTest : public ::testing::Test {
protected:
    const HohmannTransfer leo_geo = compute_hohmann(LEO_ALT, GEO_ALT);
};


// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------
TEST(DynamicsConstants, MuValue) {
    EXPECT_DOUBLE_EQ(MU, 398600.4418);
}

TEST(DynamicsConstants, REarthValue) {
    EXPECT_DOUBLE_EQ(R_EARTH, 6378.137);
}

TEST(DynamicsConstants, G0Value) {
    EXPECT_DOUBLE_EQ(G0, 9.80665);
}


// ---------------------------------------------------------------------------
// compute_hohmann - valid inputs
// ---------------------------------------------------------------------------
TEST_F(DynamicsTest, RadiiValues) {
    EXPECT_NEAR(leo_geo.r1, R1_REF, TOL);
    EXPECT_NEAR(leo_geo.r2, R2_REF, TOL);
}

TEST_F(DynamicsTest, SemiMajorAxis) {
    EXPECT_NEAR(leo_geo.a_transfer, A_REF, TOL);
}

TEST_F(DynamicsTest, DvDepartureParity) {
    // Cross-language parity: matches reference formula to 1e-6 km/s
    const auto ref = ref_hohmann(LEO_ALT, GEO_ALT);
    EXPECT_NEAR(leo_geo.dv_departure, ref.dv1, TOL);
}

TEST_F(DynamicsTest, DvArrivalParity) {
    const auto ref = ref_hohmann(LEO_ALT, GEO_ALT);
    EXPECT_NEAR(leo_geo.dv_arrival, ref.dv2, TOL);
}

TEST_F(DynamicsTest, TotalDvParity) {
    const auto ref = ref_hohmann(LEO_ALT, GEO_ALT);
    EXPECT_NEAR(leo_geo.total_dv(), ref.total, TOL);
}

TEST_F(DynamicsTest, TofParity) {
    const auto ref = ref_hohmann(LEO_ALT, GEO_ALT);
    EXPECT_NEAR(leo_geo.tof, ref.tof, TOL);
}

TEST_F(DynamicsTest, KnownLeoGeoDvDeparture) {
    // LEO-to-GEO departure burn ~2.4 km/s (Vallado Table 6-1)
    EXPECT_GT(leo_geo.dv_departure, 2.3);
    EXPECT_LT(leo_geo.dv_departure, 2.5);
}

TEST_F(DynamicsTest, KnownLeoGeoDvArrival) {
    EXPECT_GT(leo_geo.dv_arrival, 1.3);
    EXPECT_LT(leo_geo.dv_arrival, 1.6);
}

TEST_F(DynamicsTest, KnownLeoGeoTotalDv) {
    // LEO-to-GEO total dv ~3.9 km/s
    EXPECT_GT(leo_geo.total_dv(), 3.7);
    EXPECT_LT(leo_geo.total_dv(), 4.1);
}

TEST_F(DynamicsTest, KnownLeoGeoTofHours) {
    EXPECT_GT(leo_geo.tof_hours(), 4.5);
    EXPECT_LT(leo_geo.tof_hours(), 6.0);
}

TEST_F(DynamicsTest, SameAltitudeZeroDv) {
    const auto t = compute_hohmann(LEO_ALT, LEO_ALT);
    EXPECT_NEAR(t.dv_departure, 0.0, TOL);
    EXPECT_NEAR(t.dv_arrival,   0.0, TOL);
    EXPECT_NEAR(t.total_dv(),   0.0, TOL);
}

TEST_F(DynamicsTest, ReverseTransferPositiveDv) {
    const auto t = compute_hohmann(GEO_ALT, LEO_ALT);
    EXPECT_GT(t.total_dv(), 0.0);
}

TEST_F(DynamicsTest, ReverseSymmetryTotalDv) {
    const auto fwd = compute_hohmann(LEO_ALT, GEO_ALT);
    const auto rev = compute_hohmann(GEO_ALT, LEO_ALT);
    EXPECT_NEAR(fwd.total_dv(), rev.total_dv(), TOL);
}

TEST_F(DynamicsTest, ReverseSymmetryTof) {
    const auto fwd = compute_hohmann(LEO_ALT, GEO_ALT);
    const auto rev = compute_hohmann(GEO_ALT, LEO_ALT);
    EXPECT_NEAR(fwd.tof, rev.tof, TOL);
}

TEST_F(DynamicsTest, TotalDvEqualsSum) {
    EXPECT_NEAR(leo_geo.total_dv(),
                leo_geo.dv_departure + leo_geo.dv_arrival, TOL);
}

TEST_F(DynamicsTest, TofHoursConversion) {
    EXPECT_NEAR(leo_geo.tof_hours(), leo_geo.tof / 3600.0, TOL);
}

TEST_F(DynamicsTest, SemiMajorAxisBetweenRadii) {
    EXPECT_GT(leo_geo.a_transfer, leo_geo.r1);
    EXPECT_LT(leo_geo.a_transfer, leo_geo.r2);
}

TEST_F(DynamicsTest, SmallAltitudeDifference) {
    const auto t = compute_hohmann(400.0, 401.0);
    EXPECT_GT(t.total_dv(), 0.0);
    EXPECT_LT(t.total_dv(), 0.1);
}

TEST_F(DynamicsTest, LargeAltitudeDifference) {
    EXPECT_NO_THROW(compute_hohmann(400.0, 100000.0));
}


// ---------------------------------------------------------------------------
// compute_hohmann - invalid inputs
// ---------------------------------------------------------------------------
TEST(DynamicsValidation, NegativeH1Throws) {
    EXPECT_THROW(compute_hohmann(-100.0, GEO_ALT), std::invalid_argument);
}

TEST(DynamicsValidation, NegativeH2Throws) {
    EXPECT_THROW(compute_hohmann(LEO_ALT, -100.0), std::invalid_argument);
}

TEST(DynamicsValidation, ZeroH1Throws) {
    EXPECT_THROW(compute_hohmann(0.0, GEO_ALT), std::invalid_argument);
}

TEST(DynamicsValidation, ZeroH2Throws) {
    EXPECT_THROW(compute_hohmann(LEO_ALT, 0.0), std::invalid_argument);
}

TEST(DynamicsValidation, BothNegativeThrows) {
    EXPECT_THROW(compute_hohmann(-100.0, -200.0), std::invalid_argument);
}


// ---------------------------------------------------------------------------
// circular_velocity
// ---------------------------------------------------------------------------
TEST(DynamicsCircularVelocity, LeoValue) {
    const double vc = circular_velocity(R1_REF);
    EXPECT_NEAR(vc, std::sqrt(MU / R1_REF), TOL);
}

TEST(DynamicsCircularVelocity, NegativeRadiusThrows) {
    EXPECT_THROW(circular_velocity(-1.0), std::invalid_argument);
}

TEST(DynamicsCircularVelocity, ZeroRadiusThrows) {
    EXPECT_THROW(circular_velocity(0.0), std::invalid_argument);
}


// ---------------------------------------------------------------------------
// orbital_period
// ---------------------------------------------------------------------------
TEST(DynamicsOrbitalPeriod, GeoPeriodApprox24h) {
    const double T = orbital_period(R2_REF);
    EXPECT_NEAR(T / 3600.0, 24.0, 0.1);
}

TEST(DynamicsOrbitalPeriod, FormulaConsistency) {
    const double T = orbital_period(A_REF);
    const double expected = 2.0 * PI * A_REF * std::sqrt(A_REF / MU);
    EXPECT_NEAR(T, expected, TOL);
}

TEST(DynamicsOrbitalPeriod, NegativeSmaThrows) {
    EXPECT_THROW(orbital_period(-1.0), std::invalid_argument);
}


// ---------------------------------------------------------------------------
// Parameterized transfer scenarios
// ---------------------------------------------------------------------------
struct TransferScenario { double h1, h2, dv_lo, dv_hi; };

class DynamicsParameterizedTest
    : public ::testing::TestWithParam<TransferScenario> {};

TEST_P(DynamicsParameterizedTest, DeltaVInExpectedRange) {
    const auto [h1, h2, dv_lo, dv_hi] = GetParam();
    const auto t = compute_hohmann(h1, h2);
    EXPECT_GE(t.total_dv(), dv_lo);
    EXPECT_LE(t.total_dv(), dv_hi);
    EXPECT_GT(t.tof, 0.0);
}

INSTANTIATE_TEST_SUITE_P(
    ScenarioSuite,
    DynamicsParameterizedTest,
    ::testing::Values(
        TransferScenario{200.0,   400.0,   0.0, 1.0},
        TransferScenario{400.0,   35786.0, 3.0, 5.0},
        TransferScenario{35786.0, 400.0,   3.0, 5.0},
        TransferScenario{1000.0,  2000.0,  0.0, 2.0}
    )
);
