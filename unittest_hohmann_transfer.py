import unittest
import numpy as np
from hohmann_transfer import hohmfunc

class TestHohmannTransfer(unittest.TestCase):
    def test_hohmfunc(self):
        # Example test values
        dinc1, v1, hn1, hn2, hn3, dinc = 10, 7.8, 1.5, 0.9, 1.2, 20
        # Expected result calculated independently
        expected_result = 11.892394
        result = hohmfunc(dinc1, v1, hn1, hn2, hn3, dinc)
        # Assert that the actual result is close to the expected result within an acceptable margin
        self.assertAlmostEqual(result, expected_result, places=5)

    def test_invalid_input(self):
        # Test with invalid input values
        with self.assertRaises(ValueError):
            hohmfunc(0, 7.8, 1.5, 0.9, 1.2, 20)  # dinc1 cannot be zero
        with self.assertRaises(ValueError):
            hohmfunc(10, 0, 1.5, 0.9, 1.2, 20)  # v1 cannot be zero
        with self.assertRaises(ValueError):
            hohmfunc(10, 7.8, 0, 0.9, 1.2, 20)  # hn1 cannot be zero
        with self.assertRaises(ValueError):
            hohmfunc(10, 7.8, 1.5, 0, 1.2, 20)  # hn2 cannot be zero
        with self.assertRaises(ValueError):
            hohmfunc(10, 7.8, 1.5, 0.9, 0, 20)  # hn3 cannot be zero
        with self.assertRaises(ValueError):
            hohmfunc(10, 7.8, 1.5, 0.9, 1.2, 0)  # dinc cannot be zero

if __name__ == '__main__':
    unittest.main()
