
import math

import pytest
import numpy as np
from numpy.testing import assert_approx_equal

# See for C to assert_approx_equal:
# https://stackoverflow.com/questions/16839658/printf-width-specifier-to-maintain-precision-of-floating-point-value

import qlsc

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
	Test QLSC sindist against recorded values.
	'''
	assert_approx_equal(distance, qlsc.q3c.sindist(ra1, ra2, dec1, dec2))

# What sindist actually computes. The name says "sine of the distance" and the
# docstrings used to as well, but the value is the HAVERSINE, sin^2(d/2). The
# difference is not academic: sin(d) turns over at 90 degrees, so anything that read
# the old description and used the value as a distance ordering would have been wrong
# for every separation past a quadrant.

# separation in degrees
separations = [0., 1e-6, 0.5, 1., 30., 60., 89.999, 90., 90.001, 120., 179., 180.]

@pytest.mark.parametrize("separation", separations)
def test_sindist_is_the_haversine_not_the_sine(separation):
	'''
	Test that sindist(d) is sin^2(d/2), not sin(d).
	'''
	# a separation of exactly `separation` degrees, due north along a meridian
	value = qlsc.sindist(10., 0., 10., separation)

	assert value == pytest.approx(math.sin(math.radians(separation) / 2) ** 2, abs=1e-15)
	if separation not in (0., 180.):
		assert value != pytest.approx(math.sin(math.radians(separation)), abs=1e-9)

def test_sindist_increases_all_the_way_to_180_degrees():
	'''
	Test the property that makes the haversine usable where sin(d) is not: it is
	monotonic over the whole range, where sin(d) peaks at 90 degrees and falls back to
	zero at 180.
	'''
	values = [qlsc.sindist(10., 0., 10., separation) for separation in range(0, 181)]
	assert all(earlier < later for earlier, later in zip(values, values[1:]))
	assert values[0] == 0.
	assert values[-1] == pytest.approx(1.)
	assert values[90] == pytest.approx(0.5)

def test_sindist_round_trips_to_distance():
	'''
	Test the conversion the docstring gives, against the independent `distance`
	function: d = 2*asin(sqrt(sindist)).
	'''
	rng = np.random.default_rng(4)
	for _ in range(200):
		ra1, ra2 = rng.uniform(0, 360, 2)
		dec1, dec2 = np.degrees(np.arcsin(rng.uniform(-1, 1, 2)))
		recovered = math.degrees(2 * math.asin(math.sqrt(qlsc.sindist(ra1, dec1, ra2, dec2))))
		assert recovered == pytest.approx(qlsc.distance(ra1, dec1, ra2, dec2), abs=1e-9)

def test_the_documentation_no_longer_claims_a_sine():
	'''
	Test that neither the Python wrapper nor the C function still describes the value
	as the sine of the distance. This is a docstring, but it is the docstring the API
	documentation is generated from, and it sent a reader to the wrong formula.
	'''
	for documentation in (qlsc.sindist.__doc__, qlsc.q3c.sindist.__doc__):
		assert documentation is not None
		# note the leading "the ": "haversine of the angular distance" is fine and is
		# what both now say, so a bare "sine of..." substring would match it
		assert "the sine of the angular distance" not in documentation
		assert "haversine" in documentation.lower()
