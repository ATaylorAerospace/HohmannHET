import unittest
import numpy as np
import math

# Import the class to be tested
from hohmann_transfer import HohmannTransfer, R_EARTH, MU

class TestHohmannTransfer(unittest.TestCase):
    def setUp(self):
        """
        Set up test cases with some predefined orbital transfers
        """
        # Low Earth Orbit to Geostationary Orbit
        self.leo_to_geo = HohmannTransfer(initial_altitude=500, final_altitude=35786)
        
        # Medium Earth Orbit to High Earth Orbit
        self.meo_to_heo = HohmannTransfer(initial_altitude=5000, final_altitude=20000)
    
    def test_initialization(self):
        """
        Test initialization of HohmannTransfer class
        """
        # Test valid inputs
        transfer = HohmannTransfer(initial_altitude=500, final_altitude=1000)
        self.assertIsNotNone(transfer)
        
        # Test invalid inputs (negative altitudes)
        with self.assertRaises(ValueError):
            HohmannTransfer(initial_altitude=-100, final_altitude=1000)
        
        with self.assertRaises(ValueError):
            HohmannTransfer(initial_altitude=500, final_altitude=-1000)
    
    def test_orbital_radii(self):
        """
        Test calculation of orbital radii
        """
        # Verify orbital radii calculations
        self.assertAlmostEqual(self.leo_to_geo.r1, R_EARTH + 500, places=3)
        self.assertAlmostEqual(self.leo_to_geo.r2, R_EARTH + 35786, places=3)
    
    def test_transfer_orbit_semi_major_axis(self):
        """
        Test transfer orbit semi-major axis calculation
        """
        # Verify semi-major axis calculation
        expected_semi_major_axis = (self.leo_to_geo.r1 + self.leo_to_geo.r2) / 2
        self.assertAlmostEqual(
            self.leo_to_geo.a_transfer, 
            expected_semi_major_axis, 
            places=3
        )
    
    def test_delta_v_calculations(self):
        """
        Test delta-V calculations for departure and arrival burns
        """
        # Validate delta-V calculations are not zero
        self.assertNotAlmostEqual(self.leo_to_geo.delta_v_departure, 0)
        self.assertNotAlmostEqual(self.leo_to_geo.delta_v_arrival, 0)
        
        # Verify total delta-V is positive
        total_delta_v = (
            self.leo_to_geo.delta_v_departure + 
            self.leo_to_geo.delta_v_arrival
        )
        self.assertGreater(total_delta_v, 0)
    
    def test_transfer_time(self):
        """
        Test Hohmann transfer time calculation
        """
        # Transfer time should be reasonable
        transfer_time = self.leo_to_geo.calculate_transfer_time()
        self.assertGreater(transfer_time, 0)
        
        # Validate transfer time calculation using orbital mechanics formula
        expected_time = math.pi * math.sqrt(
            self.leo_to_geo.a_transfer**3 / MU
        )
        self.assertAlmostEqual(transfer_time, expected_time, places=3)
    
    def test_multiple_transfer_scenarios(self):
        """
        Test multiple transfer scenarios
        """
        # Test transfers with different altitude ranges
        scenarios = [
            (500, 1000),   # Low altitude transfer
            (5000, 20000), # Medium to high altitude
            (1000, 35786)  # LEO to Geostationary
        ]
        
        for initial_alt, final_alt in scenarios:
            with self.subTest(initial_alt=initial_alt, final_alt=final_alt):
                transfer = HohmannTransfer(initial_alt, final_alt)
                
                # Basic validation for each scenario
                self.assertGreater(transfer.delta_v_departure, 0)
                self.assertGreater(transfer.delta_v_arrival, 0)
                self.assertGreater(transfer.calculate_transfer_time(), 0)

def main():
    unittest.main()

if __name__ == '__main__':
    main()
