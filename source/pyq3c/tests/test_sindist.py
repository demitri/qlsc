
import pytest
from numpy.testing import assert_approx_equal

# See for C to assrt_approx_equal:
# https://stackoverflow.com/questions/16839658/printf-width-specifier-to-maintain-precision-of-floating-point-value

import pyq3c

# ra1, dec1, ra2, dec2, distance - all degrees
expected_results = [
	(0.0000, 0.0000, 30.0000, 0.0000, 0.0669873),
	(15.0000, 88.0000, 350.0000, 87.5000, 0.000090352045911869147),
	(15.0000, 88.0000, 15.0000, -88.0000, 0.998782025129912209849),
	(45.4500, 48.8000, 127.8900, -27.7000, 0.636512911854655527577)
]

@pytest.mark.parametrize("ra1, ra2, dec1, dec2, distance", expected_results)
def test_sindist(ra1, ra2, dec1, dec2, distance):
	'''
	Test Q3C sindist, the sine of the angular distance between two points on a sphere.
	'''
	assert_approx_equal(distance, pyq3c.sindist(ra1, ra2, dec1, dec2))
