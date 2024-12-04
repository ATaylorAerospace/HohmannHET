#include <iostream>
#include <cassert>
#include <cmath>

// Include your HohmannTransfer class code here or ensure it's accessible to this test.

void testHohmannTransfer() {
    // Test case: Transfer from a 500 km altitude orbit to a 2000 km altitude orbit
    double initial_altitude = 500.0;   // km
    double final_altitude = 2000.0;    // km

    HohmannTransfer transfer(initial_altitude, final_altitude);

    // Expected results (pre-calculated)
    const double expected_a_transfer = (6378.137 + 500.0 + 6378.137 + 2000.0) / 2.0; // km
    const double expected_v1 = std::sqrt(398600.4418 / (6378.137 + 500.0));          // km/s
    const double expected_v_transfer_1 = std::sqrt(398600.4418 * (2.0 / (6378.137 + 500.0) - 1.0 / expected_a_transfer));
    const double expected_delta_v_departure = expected_v_transfer_1 - expected_v1;   // km/s
    const double expected_v2 = std::sqrt(398600.4418 / (6378.137 + 2000.0));         // km/s
    const double expected_v_transfer_2 = std::sqrt(398600.4418 * (2.0 / (6378.137 + 2000.0) - 1.0 / expected_a_transfer));
    const double expected_delta_v_arrival = expected_v2 - expected_v_transfer_2;     // km/s
    const double expected_total_delta_v = expected_delta_v_departure + expected_delta_v_arrival; // km/s
    const double expected_transfer_time = M_PI * std::sqrt(std::pow(expected_a_transfer, 3) / 398600.4418); // seconds

    // Tolerance for floating-point comparisons
    const double TOLERANCE = 1e-3;

    // Assertions
    assert(std::abs(transfer.getTransferSemiMajorAxis() - expected_a_transfer) < TOLERANCE);
    assert(std::abs(transfer.getDeltaVDeparture() - expected_delta_v_departure) < TOLERANCE);
    assert(std::abs(transfer.getDeltaVArrival() - expected_delta_v_arrival) < TOLERANCE);
    assert(std::abs((transfer.getDeltaVDeparture() + transfer.getDeltaVArrival()) - expected_total_delta_v) < TOLERANCE);
    assert(std::abs(transfer.calculateTransferTime() - expected_transfer_time) < TOLERANCE);

    std::cout << "All HohmannTransfer tests passed successfully.\n";
}

int main() {
    try {
        testHohmannTransfer();
    } catch (const std::exception& e) {
        std::cerr << "A test failed: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
