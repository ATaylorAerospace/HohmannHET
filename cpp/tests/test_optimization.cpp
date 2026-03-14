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
    EXPECT_THROW(optimize_isp(-100.0, DV, P, ETA), std::invalid_argument);
}

TEST(OptimizeIspValidation, NegativeDvThrows) {
    EXPECT_THROW(optimize_isp(M0, -1.0, P, ETA), std::invalid_argument);
}

TEST(OptimizeIspValidation, NegativePowerThrows) {
    EXPECT_THROW(optimize_isp(M0, DV, -100.0, ETA), std::invalid_argument);
}

TEST(OptimizeIspValidation, EfficiencyZeroThrows) {
    EXPECT_THROW(optimize_isp(M0, DV, P, 0.0), std::invalid_argument);
}

TEST(OptimizeIspValidation, EfficiencyOneThrows) {
    EXPECT_THROW(optimize_isp(M0, DV, P, 1.0), std::invalid_argument);
}

TEST(OptimizeIspValidation, BadIspBoundsThrows) {
    EXPECT_THROW(optimize_isp(M0, DV, P, ETA, 3000.0, 1000.0),
                 std::invalid_argument);
}

TEST(OptimizeIspValidation, EqualIspBoundsThrows) {
    EXPECT_THROW(optimize_isp(M0, DV, P, ETA, 1500.0, 1500.0),
                 std::invalid_argument);
}
