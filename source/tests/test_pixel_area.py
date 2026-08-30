
import math

import pytest
import numpy as np
from numpy.testing import assert_approx_equal

from qlsc import QLSC

# Tests for QLSC.ipix_area, which used to ignore both of its arguments and return the
# mean pixel area, 4*pi/nbins, with a NOTE in the docstring saying so. QLSC pixels are
# not equal area, so the mean is wrong for every pixel by up to a factor of 1.91 high
# or 2.72 low.
#
# These expectations are computed from Q3C's own q3c_pixarea, but the properties
# asserted - the sum over the sphere, the shape of the variation, the ancestor
# relationship - are geometry, not a re-statement of the implementation.

# The largest pixel of a face (at its centre) over the smallest (at a cube vertex),
# in the limit of fine subdivision. Measured at depth 30 as 5.196151, and approached
# from below as the depth increases: 2.47 at depth 2, 5.03 at depth 6, 5.18 at depth 9.
EXTREME_AREA_RATIO = 3 * math.sqrt(3)

def test_areas_sum_to_the_whole_sphere():
	'''
	Test that the pixels of a scheme tile the sphere: their areas sum to 4*pi.

	This is the property the old implementation got right (it returned 4*pi/nbins) and
	the one a per-pixel calculation could most easily get wrong.
	'''
	for depth in (0, 1, 2, 3, 4, 5, 6):
		qlsc = QLSC(depth=depth)
		total = sum(qlsc.ipix_area(ipix) for ipix in range(qlsc.nbins))
		assert total == pytest.approx(4 * math.pi, rel=1e-12), f"depth={depth} areas sum to {total}, not 4*pi"

def test_pixels_are_not_equal_area():
	'''
	Test the fact the old docstring denied ("they are very close"): the areas vary by a
	factor of five, and the mean is a poor answer for any particular pixel.
	'''
	qlsc = QLSC(depth=6)
	areas = np.array([qlsc.ipix_area(ipix) for ipix in range(qlsc.nbins)])
	mean = 4 * math.pi / qlsc.nbins

	assert areas.max() / areas.min() > 5
	assert areas.max() / mean > 1.9
	assert mean / areas.min() > 2.6

def test_the_area_spread_approaches_three_root_three():
	'''
	Test that the largest pixel of a face is 3*sqrt(3) times the smallest as the
	subdivision gets fine, and that the extremes sit where the geometry says they do:
	the largest at the centre of a face, the smallest at a cube vertex.
	'''
	# the ratio approaches the limit from below as the pixels shrink, so the tolerance
	# tightens with depth: 5.19361 at depth 12, 5.196143 at 20, 5.196151 at 30
	for depth, tolerance in ((12, 3e-3), (20, 1e-5), (30, 1e-5)):
		qlsc = QLSC(depth=depth)
		centre = qlsc.ipix_area(qlsc.ang2ipix(0., 0.))                 # centre of face 1
		vertex = qlsc.ipix_area(qlsc.ang2ipix(44.9999999, 35.2643896)) # a cube vertex
		assert centre / vertex == pytest.approx(EXTREME_AREA_RATIO, abs=tolerance), \
			f"depth={depth}: face centre / cube vertex area ratio is {centre/vertex}"

	# and the spread is already large at a coarse depth, where "within a few percent"
	# might otherwise seem plausible
	coarse = QLSC(depth=2)
	areas = [coarse.ipix_area(ipix) for ipix in range(coarse.nbins)]
	assert max(areas) / min(areas) == pytest.approx(2.472, abs=1e-3)

def test_depth_zero_and_one_are_uniform():
	'''
	Test that the variation starts at depth 2: a whole face, and a face quartered, are
	all the same size as each other.
	'''
	for depth in (0, 1):
		qlsc = QLSC(depth=depth)
		areas = [qlsc.ipix_area(ipix) for ipix in range(qlsc.nbins)]
		assert max(areas) == pytest.approx(min(areas), rel=1e-12)
		assert max(areas) == pytest.approx(4 * math.pi / qlsc.nbins, rel=1e-12)

def test_the_default_depth_is_the_schemes_own():
	'''
	Test that omitting depth measures the pixel itself.
	'''
	qlsc = QLSC(depth=8)
	ipix = int(qlsc.ang2ipix(12.3, 45.6))
	assert qlsc.ipix_area(ipix) == qlsc.ipix_area(ipix, qlsc.depth)

def test_a_coarser_depth_measures_the_ancestor_pixel():
	'''
	Test the ``depth`` parameter: it names an absolute depth, and measures the ancestor
	pixel containing this one at that depth.

	Two things are checked. Each step up is close to four times the area, since a pixel
	is four of the next depth down - only close, because the four are not equal to each
	other. And the ancestor's area is exactly what the coarser scheme reports for the
	pixel that contains the same coordinate, which is what "ancestor" has to mean.
	'''
	qlsc = QLSC(depth=8)
	ra, dec = 12.3, 45.6
	ipix = int(qlsc.ang2ipix(ra, dec))

	for depth in range(0, qlsc.depth + 1):
		coarse = QLSC(depth=depth)
		assert qlsc.ipix_area(ipix, depth) == pytest.approx(coarse.ipix_area(coarse.ang2ipix(ra, dec)), rel=1e-12), \
			f"the depth={depth} ancestor's area does not match what a depth={depth} scheme reports"

	for depth in range(1, qlsc.depth + 1):
		ratio = qlsc.ipix_area(ipix, depth - 1) / qlsc.ipix_area(ipix, depth)
		assert 3.0 < ratio < 5.0, f"going up one depth changed the area by {ratio}, not roughly 4"

def test_depth_zero_is_a_whole_cube_face():
	'''
	Test the coarsest case: at depth 0 a pixel is a sixth of the sphere.
	'''
	for scheme_depth in (0, 3, 9, 30):
		qlsc = QLSC(depth=scheme_depth)
		ipix = int(qlsc.ang2ipix(200., -30.))
		assert qlsc.ipix_area(ipix, 0) == pytest.approx(4 * math.pi / 6, rel=1e-12)

# The expected values below were carried over from the disabled xtest_pixelarea.py,
# whose third argument was the RAW Q3C offset - the number of levels up from the
# scheme's resolution - rather than a depth. On a depth-30 scheme an offset of n is an
# absolute depth of 30-n, which is the translation ipix_area now performs.
#
# ipix, depth, area (steradians)
q3c_reference_values = [
	(2882303761517117440,  0, 2.0943951024e+00), # ra=45, dec=0, whole face
	(2017612633061982208,  0, 2.0943951024e+00), # ra=0,  dec=0, whole face
	(2882303761517117440, 29, 4.9065384358e-18), # ra=45, dec=0
	(2882303761517117440, 20, 1.2862215819e-12), # ra=45, dec=0
	(480243908601203291,  29, 1.3871447804e-17), # ra=10, dec=89
	(6227361032241724404, 29, 1.0580610731e-17), # ra=134, dec=-66
	(6227361032241724404, 26, 6.7715908182e-16), # ra=134, dec=-66
	(6227361032241724404, 19, 1.1094554873e-11)  # ra=134, dec=-66
]

@pytest.mark.parametrize("ipix, depth, area", q3c_reference_values)
def test_against_q3c_reference_values(ipix, depth, area):
	'''
	Test against values recorded from Q3C itself, at the default depth of 30.
	'''
	assert_approx_equal(area, QLSC(depth=30).ipix_area(ipix, depth), significant=9)

def test_a_finer_depth_than_the_scheme_is_refused():
	'''
	Test that a depth this scheme does not resolve raises.

	Without the check the underlying q3c_pixarea shifts the pixel's face coordinates
	right out of the face and returns an area LARGER than a whole cube face, silently:
	on a depth-6 scheme the raw offsets 7, 10 and 30 give 3.11, 3.62 and 3.67 steradians
	against a face's 2.09.
	'''
	qlsc = QLSC(depth=6)
	ipix = int(qlsc.ang2ipix(12., 34.))
	for depth in (7, 8, 30, 31, -1):
		with pytest.raises(ValueError):
			qlsc.ipix_area(ipix, depth)

def test_no_area_exceeds_a_cube_face():
	'''
	Test the sanity bound the missing range check used to break: no pixel of any scheme,
	measured at any legal depth, is bigger than a sixth of the sphere.
	'''
	face = 4 * math.pi / 6
	for scheme_depth in (0, 2, 5, 11, 30):
		qlsc = QLSC(depth=scheme_depth)
		ipix = int(qlsc.ang2ipix(77., -12.))
		for depth in range(0, qlsc.depth + 1):
			assert 0 < qlsc.ipix_area(ipix, depth) <= face * (1 + 1e-12)

# depth
invalid_depth_values = [1.5, "3", None if False else "x", True]

@pytest.mark.parametrize("depth", invalid_depth_values)
def test_a_non_integer_depth_is_refused(depth):
	'''
	Test that a depth that is not an integer raises rather than being truncated.
	'''
	with pytest.raises(ValueError):
		QLSC(depth=6).ipix_area(0, depth)

def test_an_invalid_ipix_is_refused():
	'''
	Test that the pixel number is validated, which it previously was not - both
	arguments were ignored entirely.
	'''
	qlsc = QLSC(depth=6)
	for ipix in (-1, qlsc.nbins, 1.5, "4"):
		with pytest.raises(ValueError):
			qlsc.ipix_area(ipix)
