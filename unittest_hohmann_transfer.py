import unittest
import numpy as np
import math
from scipy.constants import pi, g

# Import the functions from the original script
from orbital_drag_script import (
    atmospheric_density,
    drag_force,
    orbital_lifetime,
    calculate_orbital_parameters
)

class TestOrbitalDragCalculations(unittest.TestCase):
    def setUp(self):
        # Constants for testing
        self.req = 6378.137  # Earth's equatorial radius (km)
        self.mu = 398600.4418  # Earth's gravitational parameter (km^3/s^2)

    # ... [previous test methods remain the same] ...

    def test_input_validation(self):
        """Comprehensive test for input validation across all functions"""
        # Atmospheric Density Function
        with self.assertRaises(ValueError, msg="Negative altitude should raise ValueError"):
            atmospheric_density(-10)
        with self.assertRaises(ValueError, msg="Non-numeric altitude should raise ValueError"):
            atmospheric_density('invalid')

        # Drag Force Function
        with self.assertRaises(ValueError, msg="Negative altitude should raise ValueError"):
            drag_force(-300, 7.7, 10)
        with self.assertRaises(ValueError, msg="Negative velocity should raise ValueError"):
            drag_force(300, -7.7, 10)
        with self.assertRaises(ValueError, msg="Negative area should raise ValueError"):
            drag_force(300, 7.7, -10)
        with self.assertRaises(ValueError, msg="Negative drag coefficient should raise ValueError"):
            drag_force(300, 7.7, 10, -1.0)

        # Orbital Lifetime Function
        with self.assertRaises(ValueError, msg="Negative altitude should raise ValueError"):
            orbital_lifetime(-300, 7.7, 1000, 10)
        with self.assertRaises(ValueError, msg="Negative velocity should raise ValueError"):
            orbital_lifetime(300, -7.7, 1000, 10)
        with self.assertRaises(ValueError, msg="Negative mass should raise ValueError"):
            orbital_lifetime(300, 7.7, -1000, 10)
        with self.assertRaises(ValueError, msg="Negative area should raise ValueError"):
            orbital_lifetime(300, 7.7, 1000, -10)

        # Calculate Orbital Parameters Function
        with self.assertRaises(ValueError, msg="Negative altitude should raise ValueError"):
            calculate_orbital_parameters(-300)

    def test_type_checking(self):
        """Test that functions handle incorrect input types"""
        # List of test cases with invalid types
        invalid_inputs = [
            {'func': atmospheric_density, 'args': ['invalid']},
            {'func': drag_force, 'args': [-300, 'invalid', 10]},
            {'func': orbital_lifetime, 'args': [300, 7.7, 'invalid', 10]},
            {'func': calculate_orbital_parameters, 'args': ['invalid']}
        ]

        # Test each function with invalid input types
        for test_case in invalid_inputs:
            with self.assertRaises((ValueError, TypeError), 
                                   msg=f"Invalid input type should raise error for {test_case['func'].__name__}"):
                test_case['func'](*test_case['args'])

def main():
    unittest.main()

if __name__ == '__main__':
    main()
