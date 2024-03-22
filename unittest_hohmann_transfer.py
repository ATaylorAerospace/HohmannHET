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

if __name__ == '__main__':
    unittest.main()
