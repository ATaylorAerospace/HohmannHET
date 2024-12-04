import unittest
import numpy as np

# Calculation Constants (should match those in your HohmannTransfer class)
MU = 398600.4418  # Earth's gravitational parameter (km^3/s^2)
R_EARTH = 6378.137  # Earth's equatorial radius (km)

# Include your HohmannTransfer class code here or ensure it's accessible to this test.
# For example, if it's in a file named hohmann_transfer.py, you can import it:
# from hohmann_transfer import HohmannTransfer

class TestHohmannTransfer(unittest.TestCase):
    def test_hohmann_transfer(self):
        # Test case: Transfer from a 500 km altitude orbit to a 2000 km altitude orbit
        initial_altitude = 500.0   # km
        final_altitude = 2000.0    # km

        # Create an instance of HohmannTransfer
        transfer = HohmannTransfer(initial_altitude, final_altitude)

        # Expected results (pre-calculated)
        r1 = R_EARTH + initial_altitude
        r2 = R_EARTH + final_altitude
        a_transfer_expected = (r1 + r2) / 2  # km

        v1_expected = np.sqrt(MU / r1)
        v_transfer_1_expected = np.sqrt(MU * (2 / r1 - 1 / a_transfer_expected))
        delta_v_departure_expected = v_transfer_1_expected - v1_expected

        v2_expected = np.sqrt(MU / r2)
        v_transfer_2_expected = np.sqrt(MU * (2 / r2 - 1 / a_transfer_expected))
        delta_v_arrival_expected = v2_expected - v_transfer_2_expected

        total_delta_v_expected = delta_v_departure_expected + delta_v_arrival_expected

        transfer_time_expected = np.pi * np.sqrt(a_transfer_expected ** 3 / MU)  # seconds

        # Tolerance for floating-point comparisons
        tolerance = 1e-6

        # Assertions to check if computed values match expected results
        self.assertAlmostEqual(transfer.a_transfer, a_transfer_expected, delta=tolerance)
        self.assertAlmostEqual(transfer.delta_v_departure, delta_v_departure_expected, delta=tolerance)
        self.assertAlmostEqual(transfer.delta_v_arrival, delta_v_arrival_expected, delta=tolerance)
        self.assertAlmostEqual(
            transfer.delta_v_departure + transfer.delta_v_arrival,
            total_delta_v_expected,
            delta=tolerance
        )
        self.assertAlmostEqual(transfer.calculate_transfer_time(), transfer_time_expected, delta=tolerance)

if __name__ == '__main__':
    unittest.main()
