
import pytest
from numpy.testing import assert_approx_equal

# See for C to assrt_approx_equal:
# https://stackoverflow.com/questions/16839658/printf-width-specifier-to-maintain-precision-of-floating-point-value

import pyq3c

# ra1, dec1, ra2, dec2, distance - all degrees
expected_results = [
	(0.0000, 0.0000, 30.0000, 0.0000, 30.0),
	(15.0000, 88.0000, 350.0000, 87.5000, 1.089251492615371130768),
	(15.0000, 88.0000, 15.0000, -88.0000, 176.0),
	(45.4500, 48.8000, 127.8900, -27.7000, 105.844400629850568407164)
]

@pytest.mark.parametrize("ra1, ra2, dec1, dec2, distance", expected_results)
def test_distance(ra1, ra2, dec1, dec2, distance):
	'''
	Test Q3C dist, the angular distance between two points on a sphere.
	'''
	assert_approx_equal(distance, pyq3c.distance(ra1, ra2, dec1, dec2))
