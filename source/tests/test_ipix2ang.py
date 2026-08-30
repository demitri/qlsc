
import pytest
from numpy.testing import assert_approx_equal

from qlsc import QLSC

# As of qlsc v1.1 (following the Q3C PostgreSQL plugin v2.0.1+),
# ipix2ang returns the center of the pixel, not the lower left corner.

# depth, ipix, ra, dec
# The center values were computed independently from the pixel corner (ipix2xy)
# plus half a bin length, projected via xy2ang. Away from ra=0 and the poles
# they agree with the corner to well past seven significant figures.
expected_results_30 = [
	(30, 294541929365396512, 12.3, 45.6),
	(30, 6438254810785524898, 311.76, -84.33),
	(30, 2017612633061982208, 5.336085289072462e-08, 5.336085289072462e-08),
	(30, 864691128455135232, 135.0, 89.99999992453637) # at the pole pixel, the center's ra differs greatly from the corner's
]

# depth, ipix, ra, dec - pixel centers (same values as test_ipixcenter2ang.py)
expected_results_2 = [
	(2, 0, 315.0, 43.31385665828306),
	(2, 3, 315.0, 70.52877936550931),
	(2, 5, 45.0, 43.31385665828306),
	(2, 9, 225.0, 70.52877936550931),
	(2, 15, 135.0, 43.31385665828306),
	(2, 17, 345.96375653207355, -36.039893430303856),
	(2, 22, 14.036243467926479, -13.63302222536641),
	(2, 25, 345.96375653207355, 13.63302222536641),
	(2, 30, 14.036243467926479, 36.039893430303856),
	(2, 82, 251.56505117707798, -51.67118189854412),
	(2, 83, 225.0, -70.52877936550931),
	(2, 87, 108.43494882292202, -51.67118189854412),
	(2, 90, 315.0, -43.31385665828306),
	(2, 95, 45.0, -43.31385665828306)
]

expected_results = expected_results_30 + expected_results_2

@pytest.mark.parametrize("depth, ipix, ra, dec", expected_results)
def test_ipix2ang(depth, ipix, ra, dec):
	'''
	Test QLSC ipix2ang.
	'''
	q = QLSC(depth=depth) # use default nside value
	ra_out, dec_out = q.ipix2ang(ipix)

	assert_approx_equal(ra, ra_out)
	assert_approx_equal(dec, dec_out)

@pytest.mark.parametrize("depth", [2, 26, 30])
def test_ipix2ang_matches_ipix2ang_center(depth):
	'''
	Test that ipix2ang and ipix2ang_center are synonyms.
	'''
	q = QLSC(depth=depth)
	for ipix in [0, 1, q.nbins//2, q.nbins-1]:
		assert q.ipix2ang(ipix) == q.ipix2ang_center(ipix)

@pytest.mark.parametrize("depth", [26, 30])
def test_ipix2ang_face_boundaries(depth):
	'''
	Test the first and last ipix of every face at high depths.

	Regression test: the face number was previously computed via floating point,
	which rounded ipix values near face boundaries onto the wrong face at depths >= 26
	(the maximum valid ipix raised a ValueError).
	'''
	q = QLSC(depth=depth)
	bins_per_face = q.nside * q.nside
	for face in range(6):
		for ipix in (face * bins_per_face, (face + 1) * bins_per_face - 1):
			facenum, x, y = q.ipix2xy(ipix)
			assert facenum == face, f"ipix2xy({ipix}) at depth {depth} returned face {facenum}; expected {face}"
			assert -1 <= x <= 1 and -1 <= y <= 1
			ra, dec = q.ipix2ang(ipix) # must not raise
			assert 0 <= ra <= 360 and -90 <= dec <= 90
