#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <chrono>
#include <iomanip>

// Optimized constants with compile-time evaluation
namespace orbital_constants {
    constexpr double MU = 398600.4418;
    constexpr double R_EARTH = 6378.137;
    constexpr double INV_3600 = 1.0 / 3600.0;
    constexpr double PI = 3.14159265358979323846;
    constexpr double SQRT_MU = 631.348; // Pre-computed sqrt(MU)
}

class HohmannTransfer {
private:
    const double r1, r2, a_transfer, v1;
    const double delta_v_departure, delta_v_arrival, transfer_time;

    // Optimized calculation with fewer sqrt operations
    static constexpr std::tuple<double, double, double> calculateValues(
        double r1, double r2, double a_transfer, double v1) noexcept {
        
        using namespace orbital_constants;
        
        // Pre-compute common terms
        const double sqrt_mu_r1 = SQRT_MU / std::sqrt(r1);
        const double sqrt_mu_r2 = SQRT_MU / std::sqrt(r2);
        const double sqrt_mu_a = SQRT_MU / std::sqrt(a_transfer);
        
        // Optimized velocity calculations
        const double v_transfer_1 = sqrt_mu_r1 * std::sqrt(2.0 - r1 / a_transfer);
        const double v_transfer_2 = sqrt_mu_r2 * std::sqrt(2.0 - r2 / a_transfer);
        const double v2 = sqrt_mu_r2;
        
        // Optimized transfer time calculation
        const double a_cubed_half = a_transfer * std::sqrt(a_transfer);
        const double transfer_time = PI * a_cubed_half / SQRT_MU;
        
        return {v_transfer_1 - v1, v2 - v_transfer_2, transfer_time};
    }

    // Validate inputs at compile time when possible
    static constexpr void validateAltitudes(double initial_alt, double final_alt) {
        if (initial_alt <= 0.0 || final_alt <= 0.0) {
            throw std::invalid_argument("Altitudes must be positive");
        }
    }

public:
    constexpr HohmannTransfer(double initial_altitude, double final_altitude) 
        : r1(initial_altitude + orbital_constants::R_EARTH)
        , r2(final_altitude + orbital_constants::R_EARTH)
        , a_transfer((r1 + r2) * 0.5)
        , v1(orbital_constants::SQRT_MU / std::sqrt(r1))
        , delta_v_departure([&]() {
            validateAltitudes(initial_altitude, final_altitude);
            return std::get<0>(calculateValues(r1, r2, a_transfer, v1));
        }())
        , delta_v_arrival(std::get<1>(calculateValues(r1, r2, a_transfer, v1)))
        , transfer_time(std::get<2>(calculateValues(r1, r2, a_transfer, v1)))
    {}

    // Inline getters with noexcept and nodiscard
    [[nodiscard]] constexpr double getTransferTime() const noexcept { return transfer_time; }
    [[nodiscard]] constexpr double getDeltaVDeparture() const noexcept { return delta_v_departure; }
    [[nodiscard]] constexpr double getDeltaVArrival() const noexcept { return delta_v_arrival; }
    [[nodiscard]] constexpr double getTransferSemiMajorAxis() const noexcept { return a_transfer; }
    [[nodiscard]] constexpr double getTotalDeltaV() const noexcept { 
        return delta_v_departure + delta_v_arrival; 
    }

    // Optimized print function with single stream operation
    void printTransferDetails() const {
        using namespace orbital_constants;
        
        const double total_delta_v = getTotalDeltaV();
        const double transfer_hours = transfer_time * INV_3600;
        
        // Single formatted output for better performance
        std::cout << std::fixed << std::setprecision(2)
                  << "\nHohmann Transfer Orbit Details:\n"
                     "Initial Orbit Altitude: " << (r1 - R_EARTH) << " km\n"
                     "Final Orbit Altitude: " << (r2 - R_EARTH) << " km\n"
                     "Transfer Orbit Semi-Major Axis: " << a_transfer << " km\n"
                     "Delta-V (Departure Burn): " << delta_v_departure << " km/s\n"
                     "Delta-V (Arrival Burn): " << delta_v_arrival << " km/s\n"
                     "Total Delta-V: " << total_delta_v << " km/s\n"
                     "Transfer Time: " << transfer_hours << " hours\n";
    }
};

// Optimized test fixture with static constants
class HohmannTransferTest : public ::testing::Test {
protected:
    static constexpr double TOLERANCE = 1e-6;
    static constexpr double LEO_ALTITUDE = 400.0;
    static constexpr double GEO_ALTITUDE = 35786.0;

    // Inline helper for better performance
    static void ExpectNear(double actual, double expected, double tolerance = TOLERANCE) {
        EXPECT_NEAR(actual, expected, tolerance);
    }
};

// Optimized tests with reduced object creation
TEST_F(HohmannTransferTest, TestConstants) {
    using namespace orbital_constants;
    EXPECT_DOUBLE_EQ(MU, 398600.4418);
    EXPECT_DOUBLE_EQ(R_EARTH, 6378.137);
    EXPECT_DOUBLE_EQ(INV_3600, 1.0 / 3600.0);
}

TEST_F(HohmannTransferTest, ConstructorValidInputs) {
    EXPECT_NO_THROW(HohmannTransfer(LEO_ALTITUDE, GEO_ALTITUDE));
    EXPECT_NO_THROW(HohmannTransfer(100.0, 200.0));
    EXPECT_NO_THROW(HohmannTransfer(0.1, 0.2));
}

TEST_F(HohmannTransferTest, ConstructorInvalidInputs) {
    EXPECT_THROW(HohmannTransfer(-100.0, GEO_ALTITUDE), std::invalid_argument);
    EXPECT_THROW(HohmannTransfer(LEO_ALTITUDE, -100.0), std::invalid_argument);
    EXPECT_THROW(HohmannTransfer(0.0, GEO_ALTITUDE), std::invalid_argument);
    EXPECT_THROW(HohmannTransfer(LEO_ALTITUDE, 0.0), std::invalid_argument);
    EXPECT_THROW(HohmannTransfer(-100.0, -200.0), std::invalid_argument);
}

TEST_F(HohmannTransferTest, GetterMethodsPositiveValues) {
    const HohmannTransfer transfer(LEO_ALTITUDE, GEO_ALTITUDE);

    EXPECT_GT(transfer.getDeltaVDeparture(), 0.0);
    EXPECT_GT(transfer.getDeltaVArrival(), 0.0);
    EXPECT_GT(transfer.getTransferSemiMajorAxis(), orbital_constants::R_EARTH);
    EXPECT_GT(transfer.getTotalDeltaV(), 0.0);
    EXPECT_GT(transfer.getTransferTime(), 0.0);
}

TEST_F(HohmannTransferTest, TotalDeltaVCalculation) {
    const HohmannTransfer transfer(LEO_ALTITUDE, GEO_ALTITUDE);
    
    ExpectNear(transfer.getTotalDeltaV(), 
               transfer.getDeltaVDeparture() + transfer.getDeltaVArrival());
}

TEST_F(HohmannTransferTest, LEOToGEOKnownValues) {
    const HohmannTransfer transfer(LEO_ALTITUDE, GEO_ALTITUDE);

    ExpectNear(transfer.getDeltaVDeparture(), 2.44, 0.1);
    ExpectNear(transfer.getDeltaVArrival(), 1.47, 0.1);
    ExpectNear(transfer.getTotalDeltaV(), 3.91, 0.2);
    ExpectNear(transfer.getTransferTime() * orbital_constants::INV_3600, 5.25, 0.5);
}

TEST_F(HohmannTransferTest, SameAltitudeTransfer) {
    const HohmannTransfer transfer(LEO_ALTITUDE, LEO_ALTITUDE);

    ExpectNear(transfer.getDeltaVDeparture(), 0.0, TOLERANCE);
    ExpectNear(transfer.getDeltaVArrival(), 0.0, TOLERANCE);
    ExpectNear(transfer.getTotalDeltaV(), 0.0, TOLERANCE);
    EXPECT_GT(transfer.getTransferTime(), 0.0);
}

TEST_F(HohmannTransferTest, ReverseTransfer) {
    const HohmannTransfer transfer(GEO_ALTITUDE, LEO_ALTITUDE);

    EXPECT_GT(transfer.getDeltaVDeparture(), 0.0);
    EXPECT_GT(transfer.getDeltaVArrival(), 0.0);
    EXPECT_GT(transfer.getTotalDeltaV(), 0.0);
    EXPECT_GT(transfer.getTransferTime(), 0.0);
}

TEST_F(HohmannTransferTest, SemiMajorAxisCalculation) {
    constexpr double initial_alt = 200.0;
    constexpr double final_alt = 800.0;
    const HohmannTransfer transfer(initial_alt, final_alt);

    constexpr double r1 = initial_alt + orbital_constants::R_EARTH;
    constexpr double r2 = final_alt + orbital_constants::R_EARTH;
    constexpr double expected_a = (r1 + r2) * 0.5;

    ExpectNear(transfer.getTransferSemiMajorAxis(), expected_a, TOLERANCE);
}

TEST_F(HohmannTransferTest, PhysicalConstraints) {
    const HohmannTransfer transfer(LEO_ALTITUDE, GEO_ALTITUDE);

    constexpr double r1 = LEO_ALTITUDE + orbital_constants::R_EARTH;
    constexpr double r2 = GEO_ALTITUDE + orbital_constants::R_EARTH;
    
    EXPECT_GT(transfer.getTransferSemiMajorAxis(), r1);
    EXPECT_LT(transfer.getTransferSemiMajorAxis(), r2);
    ExpectNear(transfer.getTransferSemiMajorAxis(), (r1 + r2) * 0.5, TOLERANCE);
}

TEST_F(HohmannTransferTest, EdgeCases) {
    const HohmannTransfer small_transfer(400.0, 401.0);
    EXPECT_GT(small_transfer.getTotalDeltaV(), 0.0);
    EXPECT_LT(small_transfer.getTotalDeltaV(), 0.1);

    EXPECT_NO_THROW(HohmannTransfer(400.0, 100000.0));
}

TEST_F(HohmannTransferTest, MethodConsistency) {
    const HohmannTransfer transfer(LEO_ALTITUDE, GEO_ALTITUDE);

    const double deltaV1 = transfer.getTotalDeltaV();
    const double deltaV2 = transfer.getTotalDeltaV();
    EXPECT_DOUBLE_EQ(deltaV1, deltaV2);

    const double time1 = transfer.getTransferTime();
    const double time2 = transfer.getTransferTime();
    EXPECT_DOUBLE_EQ(time1, time2);
}

TEST_F(HohmannTransferTest, PrintFunctionOutput) {
    const HohmannTransfer transfer(LEO_ALTITUDE, GEO_ALTITUDE);

    std::ostringstream captured_output;
    const std::streambuf* orig = std::cout.rdbuf();
    std::cout.rdbuf(captured_output.rdbuf());

    transfer.printTransferDetails();

    std::cout.rdbuf(orig);
    const std::string output = captured_output.str();

    EXPECT_THAT(output, ::testing::HasSubstr("Hohmann Transfer Orbit Details"));
    EXPECT_THAT(output, ::testing::HasSubstr("Initial Orbit Altitude"));
    EXPECT_THAT(output, ::testing::HasSubstr("Final Orbit Altitude"));
    EXPECT_THAT(output, ::testing::HasSubstr("Delta-V"));
    EXPECT_THAT(output, ::testing::HasSubstr("Transfer Time"));
}

// Optimized performance test with better measurement
TEST_F(HohmannTransferTest, PerformanceTest) {
    constexpr int iterations = 10000; // Increased for better measurement
    
    const auto start = std::chrono::high_resolution_clock::now();
    
    double sum = 0.0; // Prevent optimization
    for (int i = 0; i < iterations; ++i) {
        const HohmannTransfer transfer(400.0 + i * 0.001, 35786.0);
        sum += transfer.getTotalDeltaV();
    }
    
    const auto end = std::chrono::high_resolution_clock::now();
    const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    // Prevent optimization of sum
    volatile double result = sum;
    (void)result;
    
    // Should be much faster now
    EXPECT_LT(duration.count(), 50000); // 50ms for 10k iterations
}

// Parameterized test remains the same but with const
class HohmannTransferParameterizedTest : public ::testing::TestWithParam<std::tuple<double, double, double, double>> {};

TEST_P(HohmannTransferParameterizedTest, VariousTransferScenarios) {
    const auto [initial_alt, final_alt, expected_min_deltaV, expected_max_deltaV] = GetParam();
    
    const HohmannTransfer transfer(initial_alt, final_alt);
    const double total_deltaV = transfer.getTotalDeltaV();
    
    EXPECT_GE(total_deltaV, expected_min_deltaV);
    EXPECT_LE(total_deltaV, expected_max_deltaV);
    EXPECT_GT(transfer.getTransferTime(), 0.0);
}

INSTANTIATE_TEST_SUITE_P(
    TransferScenarios,
    HohmannTransferParameterizedTest,
    ::testing::Values(
        std::make_tuple(200.0, 400.0, 0.0, 1.0),
        std::make_tuple(400.0, 35786.0, 3.0, 5.0),
        std::make_tuple(35786.0, 400.0, 3.0, 5.0),
        std::make_tuple(1000.0, 2000.0, 0.0, 2.0)
    )
);

TEST_F(HohmannTransferTest, MathematicalRelationships) {
    const HohmannTransfer transfer1(400.0, 800.0);
    const HohmannTransfer transfer2(800.0, 1600.0);
    
    EXPECT_GT(transfer2.getTotalDeltaV(), 0.0);
    EXPECT_GT(transfer1.getTotalDeltaV(), 0.0);
}

TEST_F(HohmannTransferTest, TransferSymmetry) {
    const HohmannTransfer forward(400.0, 35786.0);
    const HohmannTransfer reverse(35786.0, 400.0);
    
    ExpectNear(forward.getTotalDeltaV(), reverse.getTotalDeltaV(), 0.1);
    ExpectNear(forward.getTransferTime(), reverse.getTransferTime(), TOLERANCE);
}

TEST_F(HohmannTransferTest, BoundaryConditions) {
    EXPECT_NO_THROW(HohmannTransfer(0.001, 0.002));
    
    const HohmannTransfer large_diff(100.0, 100000.0);
    EXPECT_GT(large_diff.getTotalDeltaV(), 10.0);
}

TEST_F(HohmannTransferTest, NumericalStability) {
    const HohmannTransfer close_transfer(1000.0, 1000.001);
    
    EXPECT_GT(close_transfer.getDeltaVDeparture(), 0.0);
    EXPECT_GT(close_transfer.getDeltaVArrival(), 0.0);
    EXPECT_LT(close_transfer.getTotalDeltaV(), 0.001);
    
    EXPECT_TRUE(std::isfinite(close_transfer.getTotalDeltaV()));
    EXPECT_TRUE(std::isfinite(close_transfer.getTransferTime()));
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
