#include <gtest/gtest.h>
#include <cmath>
#include <stdexcept>

// Test fixture for Hohmann Transfer tests
class HohmannTransferTest : public ::testing::Test {
protected:
    // Known constants for verification
    const double MU = 398600.4418;
    const double R_EARTH = 6378.137;
    
    // Test tolerance for floating-point comparisons
    const double EPSILON = 1e-6;
    
    // Common test scenarios
    const double LEO_ALTITUDE = 400.0;      // Low Earth Orbit
    const double GEO_ALTITUDE = 35786.0;    // Geostationary Orbit
    const double MEO_ALTITUDE = 20200.0;    // Medium Earth Orbit (GPS)
    
    HohmannTransfer* leo_to_geo;
    HohmannTransfer* leo_to_meo;
    
    void SetUp() override {
        leo_to_geo = new HohmannTransfer(LEO_ALTITUDE, GEO_ALTITUDE);
        leo_to_meo = new HohmannTransfer(LEO_ALTITUDE, MEO_ALTITUDE);
    }
    
    void TearDown() override {
        delete leo_to_geo;
        delete leo_to_meo;
    }
    
    // Helper function to calculate expected circular orbit velocity
    double calculateCircularVelocity(double radius) {
        return std::sqrt(MU / radius);
    }
    
    // Helper function to calculate expected transfer orbit velocity at a point
    double calculateTransferVelocity(double r1, double r2, double r_current) {
        double a_transfer = (r1 + r2) / 2.0;
        return std::sqrt(MU * (2.0/r_current - 1.0/a_transfer));
    }
};

// Test constructor validation
TEST_F(HohmannTransferTest, ConstructorValidation) {
    EXPECT_THROW(HohmannTransfer(-100, 1000), std::invalid_argument);
    EXPECT_THROW(HohmannTransfer(0, 1000), std::invalid_argument);
    EXPECT_THROW(HohmannTransfer(1000, -100), std::invalid_argument);
    EXPECT_THROW(HohmannTransfer(1000, 0), std::invalid_argument);
}

// Test LEO to GEO transfer calculations
TEST_F(HohmannTransferTest, LEOtoGEODeltaV) {
    double r1 = LEO_ALTITUDE + R_EARTH;
    double r2 = GEO_ALTITUDE + R_EARTH;
    
    // Calculate expected values
    double v1 = calculateCircularVelocity(r1);
    double v2 = calculateCircularVelocity(r2);
    double v_transfer_1 = calculateTransferVelocity(r1, r2, r1);
    double v_transfer_2 = calculateTransferVelocity(r1, r2, r2);
    
    double expected_delta_v1 = v_transfer_1 - v1;
    double expected_delta_v2 = v2 - v_transfer_2;
    
    EXPECT_NEAR(leo_to_geo->getDeltaVDeparture(), expected_delta_v1, EPSILON);
    EXPECT_NEAR(leo_to_geo->getDeltaVArrival(), expected_delta_v2, EPSILON);
}

// Test semi-major axis calculation
TEST_F(HohmannTransferTest, TransferOrbitParameters) {
    double expected_sma_leo_geo = (LEO_ALTITUDE + R_EARTH + GEO_ALTITUDE + R_EARTH) / 2.0;
    EXPECT_NEAR(leo_to_geo->getTransferSemiMajorAxis(), expected_sma_leo_geo, EPSILON);
}

// Test transfer time calculation
TEST_F(HohmannTransferTest, TransferTime) {
    double r1 = LEO_ALTITUDE + R_EARTH;
    double r2 = GEO_ALTITUDE + R_EARTH;
    double a_transfer = (r1 + r2) / 2.0;
    
    double expected_time = M_PI * std::sqrt(std::pow(a_transfer, 3) / MU);
    double actual_time = leo_to_geo->calculateTransferTime();
    
    EXPECT_NEAR(actual_time, expected_time, EPSILON);
}

// Test realistic mission scenarios
TEST_F(HohmannTransferTest, RealisticMissionScenarios) {
    // Test LEO to MEO (typical GPS satellite transfer)
    EXPECT_GT(leo_to_meo->getDeltaVDeparture(), 0);
    EXPECT_GT(leo_to_meo->getDeltaVArrival(), 0);
    EXPECT_GT(leo_to_meo->calculateTransferTime(), 0);
    
    // Verify total delta-V is reasonable (should be less than direct transfer)
    double total_dv_leo_meo = leo_to_meo->getDeltaVDeparture() + 
                             leo_to_meo->getDeltaVArrival();
    double direct_v = std::sqrt(MU / (MEO_ALTITUDE + R_EARTH)) - 
                     std::sqrt(MU / (LEO_ALTITUDE + R_EARTH));
    EXPECT_LT(total_dv_leo_meo, std::abs(direct_v));
}

// Test edge cases
TEST_F(HohmannTransferTest, EdgeCases) {
    // Very close orbits
    HohmannTransfer small_change(400, 401);
    EXPECT_NEAR(small_change.getDeltaVDeparture(), 0.0, 0.1);
    EXPECT_NEAR(small_change.getDeltaVArrival(), 0.0, 0.1);
    
    // Large orbit change (but still realistic)
    HohmannTransfer large_change(400, 100000);
    EXPECT_GT(large_change.getDeltaVDeparture(), 0);
    EXPECT_GT(large_change.getDeltaVArrival(), 0);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
