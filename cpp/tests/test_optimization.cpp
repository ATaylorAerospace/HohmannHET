// Author: A Taylor | Purpose: GoogleTest suite for optimization module | Ref: Vallado/Curtis
/**
 * @file test_optimization.cpp
 * @brief GoogleTest suite for hohmann_het::optimization.
 *
 * Validates golden-section search correctness and mission optimization
 * physics behavior.
 *
 * @author A Taylor
 */

#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "hohmann_het/optimization.hpp"

using namespace hohmann_het;

static constexpr double TOL  = 1e-6;
static constexpr double M0   = 1000.0;   // kg
static constexpr double DV   = 3852.6;   // m/s
static constexpr double P    = 5000.0;   // W
static constexpr double ETA  = 0.55;


// ---------------------------------------------------------------------------
// golden_section_minimize
// ---------------------------------------------------------------------------
TEST(GoldenSectionTest, QuadraticMinimum) {
    // Minimum of f(x) = (x-3)^2 on [0, 10]
    const double x_min = golden_section_minimize(
        [](double x) { return (x - 3.0) * (x - 3.0); }, 0.0, 10.0);
    EXPECT_NEAR(x_min, 3.0, 1e-6);
}

TEST(GoldenSectionTest, OffsetQuadratic) {
    const double x_min = golden_section_minimize(
        [](double x) { return (x - 1.5) * (x - 1.5); }, -5.0, 5.0);
    EXPECT_NEAR(x_min, 1.5, 1e-6);
}

TEST(GoldenSectionTest, TighterTolerance) {
    constexpr double target = 2.718281828;
    const double x_min = golden_section_minimize(
        [](double x) { return (x - 2.718281828) * (x - 2.718281828); },
        0.0, 5.0, 1e-12);
    EXPECT_NEAR(x_min, target, 1e-9);
}

TEST(GoldenSectionTest, NarrowInterval) {
    const double x_min = golden_section_minimize(
        [](double x) { return (x - 7.3) * (x - 7.3); }, 7.0, 8.0);
    EXPECT_NEAR(x_min, 7.3, 1e-6);
}

TEST(GoldenSectionTest, InvalidBoundsThrows) {
    EXPECT_THROW(golden_section_minimize([](double x) { return x * x; }, 5.0, 1.0),
                 std::invalid_argument);
    EXPECT_THROW(golden_section_minimize([](double x) { return x * x; }, 2.0, 2.0),
                 std::invalid_argument);
}

TEST(GoldenSectionTest, NonPositiveTolThrows) {
    EXPECT_THROW(golden_section_minimize([](double x) { return x * x; }, 0.0, 5.0, 0.0),
                 std::invalid_argument);
    EXPECT_THROW(golden_section_minimize([](double x) { return x * x; }, 0.0, 5.0, -1e-6),
                 std::invalid_argument);
}

TEST(GoldenSectionTest, NonFiniteThrows) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    EXPECT_THROW(golden_section_minimize([](double x) { return x * x; }, nan, 5.0),
                 std::invalid_argument);
    EXPECT_THROW(golden_section_minimize([](double x) { return x * x; }, 0.0, inf),
                 std::invalid_argument);
    EXPECT_THROW(golden_section_minimize([](double x) { return x * x; }, 0.0, 5.0, nan),
                 std::invalid_argument);
}


// ---------------------------------------------------------------------------
// optimize_isp
// ---------------------------------------------------------------------------
class OptimizeIspTest : public ::testing::Test {
protected:
    const OptimizationResult result = optimize_isp(M0, DV, P, ETA);
};

TEST_F(OptimizeIspTest, IspWithinBounds) {
    const auto r = optimize_isp(M0, DV, P, ETA, 500.0, 5000.0);
    EXPECT_GE(r.optimal_isp, 500.0);
    EXPECT_LE(r.optimal_isp, 5000.0);
}

TEST_F(OptimizeIspTest, PropellantMassPositive) {
    EXPECT_GT(result.propellant_mass, 0.0);
}

TEST_F(OptimizeIspTest, PropellantMassLessThanInitial) {
    EXPECT_LT(result.propellant_mass, M0);
}

TEST_F(OptimizeIspTest, BurnTimePositive) {
    EXPECT_GT(result.burn_time, 0.0);
}

TEST_F(OptimizeIspTest, ObjectiveValueFinite) {
    EXPECT_TRUE(std::isfinite(result.objective_value));
}

TEST_F(OptimizeIspTest, OptimalIspPositive) {
    EXPECT_GT(result.optimal_isp, 0.0);
}

TEST_F(OptimizeIspTest, HigherLambdaLowersIsp) {
    // Heavier burn-time penalty pushes optimizer toward lower Isp
    const auto r_low  = optimize_isp(M0, DV, P, ETA, 500.0, 5000.0, 1e-6);
    const auto r_high = optimize_isp(M0, DV, P, ETA, 500.0, 5000.0, 1e-1);
    EXPECT_LE(r_high.optimal_isp, r_low.optimal_isp);
}

TEST_F(OptimizeIspTest, PropMassConsistencyWithTsiolkovsky) {
    // Verify propellant_mass in result matches Tsiolkovsky formula
    const double ve = result.optimal_isp * constants::G0;
    const double mp_check = M0 * (1.0 - std::exp(-DV / ve));
    EXPECT_NEAR(result.propellant_mass, mp_check, TOL);
}

TEST_F(OptimizeIspTest, BurnTimeConsistencyWithFormula) {
    const double ve      = result.optimal_isp * constants::G0;
    const double tb_check = M0 * DV * ve / (2.0 * ETA * P);
    EXPECT_NEAR(result.burn_time, tb_check, TOL);
}


// ---------------------------------------------------------------------------
// optimize_isp - invalid inputs
// ---------------------------------------------------------------------------
TEST(OptimizeIspValidation, NegativeMassThrows) {
    EXPECT_THROW(static_cast<void>(optimize_isp(-100.0, DV, P, ETA)), std::invalid_argument);
}

TEST(OptimizeIspValidation, NegativeDvThrows) {
    EXPECT_THROW(static_cast<void>(optimize_isp(M0, -1.0, P, ETA)), std::invalid_argument);
}

TEST(OptimizeIspValidation, ZeroDvThrows) {
    // delta_v == 0 would make the burn-time normalisation 0/0 (NaN).
    EXPECT_THROW(static_cast<void>(optimize_isp(M0, 0.0, P, ETA)), std::invalid_argument);
}

TEST(OptimizeIspValidation, NegativePowerThrows) {
    EXPECT_THROW(static_cast<void>(optimize_isp(M0, DV, -100.0, ETA)), std::invalid_argument);
}

TEST(OptimizeIspValidation, EfficiencyZeroThrows) {
    EXPECT_THROW(static_cast<void>(optimize_isp(M0, DV, P, 0.0)), std::invalid_argument);
}

TEST(OptimizeIspValidation, EfficiencyOneThrows) {
    EXPECT_THROW(static_cast<void>(optimize_isp(M0, DV, P, 1.0)), std::invalid_argument);
}

TEST(OptimizeIspValidation, BadIspBoundsThrows) {
    EXPECT_THROW(static_cast<void>(optimize_isp(M0, DV, P, ETA, 3000.0, 1000.0)),
                 std::invalid_argument);
}

TEST(OptimizeIspValidation, EqualIspBoundsThrows) {
    EXPECT_THROW(static_cast<void>(optimize_isp(M0, DV, P, ETA, 1500.0, 1500.0)),
                 std::invalid_argument);
}

TEST(OptimizeIspValidation, NegativeLambdaThrows) {
    EXPECT_THROW(static_cast<void>(optimize_isp(M0, DV, P, ETA, 800.0, 4000.0, -0.1)),
                 std::invalid_argument);
}

TEST(OptimizeIspValidation, ZeroLambdaAllowed) {
    const auto r = optimize_isp(M0, DV, P, ETA, 800.0, 4000.0, 0.0);
    EXPECT_NEAR(r.optimal_isp, 4000.0, 1e-4);
}

TEST(OptimizeIspValidation, NonFiniteThrows) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    EXPECT_THROW(static_cast<void>(optimize_isp(nan, DV, P, ETA)), std::invalid_argument);
    EXPECT_THROW(static_cast<void>(optimize_isp(M0, inf, P, ETA)), std::invalid_argument);
    EXPECT_THROW(static_cast<void>(optimize_isp(M0, DV, P, ETA, 800.0, 4000.0, nan)),
                 std::invalid_argument);
}


// ---------------------------------------------------------------------------
// min_propellant_transfer
// ---------------------------------------------------------------------------
TEST(MinPropellantTransferTest, EndToEndLeoGeo) {
    const auto r = min_propellant_transfer(400.0, 35786.0, M0, P, ETA);
    EXPECT_GT(r.propellant_mass, 0.0);
    EXPECT_GT(r.optimal_isp, 0.0);
}

TEST(MinPropellantTransferTest, ConsistencyWithOptimizeIsp) {
    const auto transfer = compute_hohmann(400.0, 35786.0);
    const double dv_ms  = transfer.total_dv() * 1000.0;
    const auto r_direct = optimize_isp(M0, dv_ms, P, ETA);
    const auto r_wrap   = min_propellant_transfer(400.0, 35786.0, M0, P, ETA);
    EXPECT_NEAR(r_direct.optimal_isp,     r_wrap.optimal_isp,     TOL);
    EXPECT_NEAR(r_direct.propellant_mass, r_wrap.propellant_mass, TOL);
}

TEST(MinPropellantTransferTest, LambdaWeightForwarding) {
    const auto transfer = compute_hohmann(400.0, 35786.0);
    const double dv_ms  = transfer.total_dv() * 1000.0;
    const auto r_direct = optimize_isp(M0, dv_ms, P, ETA, 500.0, 5000.0, 0.1);
    const auto r_wrap   = min_propellant_transfer(
        400.0, 35786.0, M0, P, ETA, 500.0, 5000.0, 0.1);
    EXPECT_NEAR(r_direct.optimal_isp,     r_wrap.optimal_isp,     TOL);
    EXPECT_NEAR(r_direct.propellant_mass, r_wrap.propellant_mass, TOL);
}


// ---------------------------------------------------------------------------
// Golden reference values (shared verbatim with Python and MATLAB suites)
// ---------------------------------------------------------------------------
TEST(OptimizationGolden, InteriorOptimumValues) {
    const auto r = optimize_isp(M0, DV, P, ETA, 500.0, 5000.0, 0.1);
    EXPECT_NEAR(r.optimal_isp,     1187.9261754441786, 1e-4);
    EXPECT_NEAR(r.propellant_mass, 281.5845876650631,  1e-4);
    EXPECT_NEAR(r.objective_value, 0.5191698227538988, 1e-6);
}
