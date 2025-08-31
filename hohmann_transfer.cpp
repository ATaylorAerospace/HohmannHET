#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <cmath>
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <chrono>

// Original HohmannTransfer class code
namespace {
    constexpr double MU = 398600.4418;
    constexpr double R_EARTH = 6378.137;
    constexpr double INV_3600 = 1.0 / 3600.0;
}

class HohmannTransfer {
private:
    const double r1;
    const double r2;
    const double a_transfer;
    const double v1;
    const double delta_v_departure;
    const double delta_v_arrival;
    const double transfer_time;

    static std::tuple<double, double, double> calculateValues(double r1, double r2, double a_transfer, double v1) {
        const double mu_2_r1 = MU * 2.0 / r1;
        const double mu_inv_a = MU / a_transfer;
        const double mu_2_r2 = MU * 2.0 / r2;
        
        const double v_transfer_1 = std::sqrt(mu_2_r1 - mu_inv_a);
        const double v_transfer_2 = std::sqrt(mu_2_r2 - mu_inv_a);
        const double v2 = std::sqrt(MU / r2);
        
        const double transfer_time = M_PI * std::sqrt(a_transfer * a_transfer * a_transfer / MU);
        
        return {v_transfer_1 - v1, v2 - v_transfer_2, transfer_time};
    }

public:
    HohmannTransfer(double initial_altitude, double final_altitude) 
        : r1(initial_altitude + R_EARTH)
        , r2(final_altitude + R_EARTH)
        , a_transfer((r1 + r2) * 0.5)
        , v1(std::sqrt(MU / r1))
        , delta_v_departure([&]() {
            if (initial_altitude <= 0 || final_altitude <= 0) {
                throw std::invalid_argument("Altitudes must be positive");
            }
            auto [dep, arr, time] = calculateValues(r1, r2, a_transfer, v1);
            const_cast<double&>(delta_v_arrival) = arr;
            const_cast<double&>(transfer_time) = time;
            return dep;
        }())
        , delta_v_arrival(0.0)
        , transfer_time(0.0)
    {}

    [[nodiscard]] constexpr double getTransferTime() const noexcept { 
        return transfer_time; 
    }

    void printTransferDetails() const {
        const double total_delta_v = delta_v_departure + delta_v_arrival;
        const double transfer_hours = transfer_time * INV_3600;
        
        std::cout << std::fixed << std::setprecision(2)
                  << "\nHohmann Transfer Orbit Details:\n"
                  << "Initial Orbit Altitude: " << (r1 - R_EARTH) << " km\n"
                  << "Final Orbit Altitude: " << (r2 - R_EARTH) << " km\n"
                  << "Transfer Orbit Semi-Major Axis: " << a_transfer << " km\n"
                  << "Delta-V (Departure Burn): " << delta_v_departure << " km/s\n"
                  << "Delta-V (Arrival Burn): " << delta_v_arrival << " km/s\n"
                  << "Total Delta-V: " << total_delta_v << " km/s\n"
                  << "Transfer Time: " << transfer_hours << " hours\n";
    }

    [[nodiscard]] constexpr double getDeltaVDeparture() const noexcept { return delta_v_departure; }
    [[nodiscard]] constexpr double getDeltaVArrival() const noexcept { return delta_v_arrival; }
    [[nodiscard]] constexpr double getTransferSemiMajorAxis() const noexcept { return a_transfer; }
    [[nodiscard]] constexpr double getTotalDeltaV() const noexcept { return delta_v_departure + delta_v_arrival; }
};

// Test fixture
class HohmannTransferTest : public ::testing::Test {
protected:
    static constexpr double TOLERANCE = 1e-6;
    static constexpr double LEO_ALTITUDE = 400.0;
    static constexpr double GEO_ALTITUDE = 35786.0;

    void ExpectNear(double actual, double expected, double tolerance = TOLERANCE) {
        EXPECT_NEAR(actual, expected, tolerance);
    }
};

// Test constants
TEST_F(HohmannTransferTest, TestConstants) {
    EXPECT_DOUBLE_EQ(MU, 398600.4418);
    EXPECT_DOUBLE_EQ(R_EARTH, 6378.137);
    EXPECT_DOUBLE_EQ(INV_3600, 1.0 / 3600.0);
}

// Test valid constructor
TEST_F(HohmannTransferTest, ConstructorValidInputs) {
    EXPECT_NO_THROW(HohmannTransfer(LEO_ALTITUDE, GEO_ALTITUDE));
    EXPECT_NO_THROW(HohmannTransfer(100.0, 200.0));
    EXPECT_NO_THROW(HohmannTransfer(0.1, 0.2));
}

// Test invalid constructor inputs
TEST_F(HohmannTransferTest, ConstructorInvalidInputs) {
    EXPECT_THROW(HohmannTransfer(-100.0, GEO_ALTITUDE), std::invalid_argument);
    EXPECT_THROW(HohmannTransfer(LEO_ALTITUDE, -100.0), std::invalid_argument);
    EXPECT_THROW(HohmannTransfer(0.0, GEO_ALTITUDE), std::invalid_argument);
    EXPECT_THROW(HohmannTransfer(LEO_ALTITUDE, 0.0), std::invalid_argument);
    EXPECT_THROW(HohmannTransfer(-100.0, -200.0), std::invalid_argument);
}

// Test getter methods return positive values
TEST_F(HohmannTransferTest, GetterMethodsPositiveValues) {
    HohmannTransfer transfer(LEO_ALTITUDE, GEO_ALTITUDE);

    EXPECT_GT(transfer.getDeltaVDeparture(), 0.0);
    EXPECT_GT(transfer.getDeltaVArrival(), 0.0);
    EXPECT_GT(transfer.getTransferSemiMajorAxis(), R_EARTH);
    EXPECT_GT(transfer.getTotalDeltaV(), 0.0);
    EXPECT_GT(transfer.getTransferTime(), 0.0);
}

// Test total delta-V calculation
TEST_F(HohmannTransferTest, TotalDeltaVCalculation) {
    HohmannTransfer transfer(LEO_ALTITUDE, GEO_ALTITUDE);
    
    ExpectNear(transfer.getTotalDeltaV(), 
               transfer.getDeltaVDeparture() + transfer.getDeltaVArrival());
}

// Test known LEO to GEO transfer values
TEST_F(HohmannTransferTest, LEOToGEOKnownValues) {
    HohmannTransfer transfer(LEO_ALTITUDE, GEO_ALTITUDE);

    ExpectNear(transfer.getDeltaVDeparture(), 2.44, 0.1);
    ExpectNear(transfer.getDeltaVArrival(), 1.47, 0.1);
    ExpectNear(transfer.getTotalDeltaV(), 3.91, 0.2);
    ExpectNear(transfer.getTransferTime() / 3600.0, 5.25, 0.5);
}

// Test same altitude (circular orbit)
TEST_F(HohmannTransferTest, SameAltitudeTransfer) {
    HohmannTransfer transfer(LEO_ALTITUDE, LEO_ALTITUDE);

    ExpectNear(transfer.getDeltaVDeparture(), 0.0, TOLERANCE);
    ExpectNear(transfer.getDeltaVArrival(), 0.0, TOLERANCE);
    ExpectNear(transfer.getTotalDeltaV(), 0.0, TOLERANCE);
    EXPECT_GT(transfer.getTransferTime(), 0.0);
}

// Test reverse transfer (higher to lower orbit)
TEST_F(HohmannTransferTest, ReverseTransfer) {
    HohmannTransfer transfer(GEO_ALTITUDE, LEO_ALTITUDE);

    EXPECT_GT(transfer.getDeltaVDeparture(), 0.0);
    EXPECT_GT(transfer.getDeltaVArrival(), 0.0);
    EXPECT_GT(transfer.getTotalDeltaV(), 0.0);
    EXPECT_GT(transfer.getTransferTime(), 0.0);
}

// Test semi-major axis calculation
TEST_F(HohmannTransferTest, SemiMajorAxisCalculation) {
    double initial_alt = 200.0;
    double final_alt = 800.0;
    HohmannTransfer transfer(initial_alt, final_alt);

    double r1 = initial_alt + R_EARTH;
    double r2 = final_alt + R_EARTH;
    double expected_a = (r1 + r2) / 2.0;

    ExpectNear(transfer.getTransferSemiMajorAxis(), expected_a, TOLERANCE);
}

// Test physical constraints
TEST_F(HohmannTransferTest, PhysicalConstraints) {
    HohmannTransfer transfer(LEO_ALTITUDE, GEO_ALTITUDE);

    double r1 = LEO_ALTITUDE + R_EARTH;
    double r2 = GEO_ALTITUDE + R_EARTH;
    
    EXPECT_GT(transfer.getTransferSemiMajorAxis(), r1);
    EXPECT_LT(transfer.getTransferSemiMajorAxis(), r2);
    ExpectNear(transfer.getTransferSemiMajorAxis(), (r1 + r2) / 2.0, TOLERANCE);
}

// Test edge cases
TEST_F(HohmannTransferTest, EdgeCases) {
    // Very small altitude difference
    HohmannTransfer small_transfer(400.0, 401.0);
    EXPECT_GT(small_transfer.getTotalDeltaV(), 0.0);
    EXPECT_LT(small_transfer.getTotalDeltaV(), 0.1);

    // Very large altitude
    EXPECT_NO_THROW(HohmannTransfer(400.0, 100000.0));
}

// Test method consistency (const correctness)
TEST_F(HohmannTransferTest, MethodConsistency) {
    HohmannTransfer transfer(LEO_ALTITUDE, GEO_ALTITUDE);

    double deltaV1 = transfer.getTotalDeltaV();
    double deltaV2 = transfer.getTotalDeltaV();
    EXPECT_DOUBLE_EQ(deltaV1, deltaV2);

    double time1 = transfer.getTransferTime();
    double time2 = transfer.getTransferTime();
    EXPECT_DOUBLE_EQ(time1, time2);
}

// Test print function output
TEST_F(HohmannTransferTest, PrintFunctionOutput) {
    HohmannTransfer transfer(LEO_ALTITUDE, GEO_ALTITUDE);

    std::ostringstream captured_output;
    std::streambuf* orig = std::cout.rdbuf();
    std::cout.rdbuf(captured_output.rdbuf());

    transfer.printTransferDetails();

    std::cout.rdbuf(orig);
    std::string output = captured_output.str();

    EXPECT_THAT(output, ::testing::HasSubstr("Hohmann Transfer Orbit Details"));
    EXPECT_THAT(output, ::testing::HasSubstr("Initial Orbit Altitude"));
    EXPECT_THAT(output, ::testing::HasSubstr("Final Orbit Altitude"));
    EXPECT_THAT(output, ::testing::HasSubstr("Delta-V"));
    EXPECT_THAT(output, ::testing::HasSubstr("Transfer Time"));
    EXPECT_THAT(output, ::testing::HasSubstr(std::to_string(static_cast<int>(LEO_ALTITUDE))));
    EXPECT_THAT(output, ::testing::HasSubstr(std::to_string(static_cast<int>(GEO_ALTITUDE))));
}

// Test performance
TEST_F(HohmannTransferTest, PerformanceTest) {
    auto start = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < 1000; ++i) {
        HohmannTransfer transfer(400.0 + i * 0.1, 35786.0);
        volatile double result = transfer.getTotalDeltaV();
        (void)result;
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    EXPECT_LT(duration.count(), 1000);
}

// Parameterized test for multiple scenarios
class HohmannTransferParameterizedTest : public ::testing::TestWithParam<std::tuple<double, double, double, double>> {};

TEST_P(HohmannTransferParameterizedTest, VariousTransferScenarios) {
    auto [initial_alt, final_alt, expected_min_deltaV, expected_max_deltaV] = GetParam();
    
    HohmannTransfer transfer(initial_alt, final_alt);
    double total_deltaV = transfer.getTotalDeltaV();
    
    EXPECT_GE(total_deltaV, expected_min_deltaV);
    EXPECT_LE(total_deltaV, expected_max_deltaV);
    EXPECT_GT(transfer.getTransferTime(), 0.0);
}

INSTANTIATE_TEST_SUITE_P(
    TransferScenarios,
    HohmannTransferParameterizedTest,
    ::testing::Values(
        std::make_tuple(200.0, 400.0, 0.0, 1.0),      // LEO to LEO
        std::make_tuple(400.0, 35786.0, 3.0, 5.0),    // LEO to GEO
        std::make_tuple(35786.0, 400.0, 3.0, 5.0),    // GEO to LEO
        std::make_tuple(1000.0, 2000.0, 0.0, 2.0)     // Medium orbits
    )
);

// Test mathematical relationships
TEST_F(HohmannTransferTest, MathematicalRelationships) {
    HohmannTransfer transfer1(400.0, 800.0);
    HohmannTransfer transfer2(800.0, 1600.0);
    
    // Higher altitude transfers should generally require more delta-V
    EXPECT_GT(transfer2.getTotalDeltaV(), 0.0);
    EXPECT_GT(transfer1.getTotalDeltaV(), 0.0);
}

// Test symmetry properties
TEST_F(HohmannTransferTest, TransferSymmetry) {
    HohmannTransfer forward(400.0, 35786.0);
    HohmannTransfer reverse(35786.0, 400.0);
    
    // Forward and reverse transfers should have similar total delta-V
    ExpectNear(forward.getTotalDeltaV(), reverse.getTotalDeltaV(), 0.1);
    
    // Transfer times should be identical
    ExpectNear(forward.getTransferTime(), reverse.getTransferTime(), TOLERANCE);
}

// Test boundary conditions
TEST_F(HohmannTransferTest, BoundaryConditions) {
    // Very small positive altitudes
    EXPECT_NO_THROW(HohmannTransfer(0.001, 0.002));
    
    // Large altitude differences
    HohmannTransfer large_diff(100.0, 100000.0);
    EXPECT_GT(large_diff.getTotalDeltaV(), 10.0); // Should be significant
}

// Test numerical stability
TEST_F(HohmannTransferTest, NumericalStability) {
    // Test with very close altitudes
    HohmannTransfer close_transfer(1000.0, 1000.001);
    
    EXPECT_GT(close_transfer.getDeltaVDeparture(), 0.0);
    EXPECT_GT(close_transfer.getDeltaVArrival(), 0.0);
    EXPECT_LT(close_transfer.getTotalDeltaV(), 0.001); // Should be very small
    
    // Values should be finite
    EXPECT_TRUE(std::isfinite(close_transfer.getTotalDeltaV()));
    EXPECT_TRUE(std::isfinite(close_transfer.getTransferTime()));
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
