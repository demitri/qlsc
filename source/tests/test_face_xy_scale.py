
import pytest
import numpy as np

from qlsc import QLSC

# The face (x,y) coordinate systems in this package are NOT all the same scale.
#
#   ipix2xy, xy2ang, xy2ipix   -> [-1,1] edge to edge
#   ang2ipix_xy                -> [-0.5,0.5] edge to edge
#
# exactly a factor of two apart. Q3C's internal convention is the half-scale one and it
# surfaces through ang2ipix_xy; everything else in the package doubles it. Nothing
# raises when the two are mixed, because a value on [-0.5,0.5] is a perfectly legal
# value on [-1,1] - it just names a different point, at half the tangent of the angle.
#
# These tests pin the discrepancy rather than hide it. It is documented in
# ang2ipix_xy's docstring with the factor of two spelled out, and it is queued as a
# breaking change for a 2.0 release (AI/FEATURE_IDEAS.md). When that happens these
# tests must be updated deliberately, which is the point of having them.

depths = [1, 5, 12, 20, 30]

# ra, dec
positions = [
	(10., 20.),
	(0., 0.),
	(44., 30.),
	(200., -70.),
	(0., 89.),
	(359.5, -0.25),
	(135.5, 44.9)
]

@pytest.mark.parametrize("depth", depths)
@pytest.mark.parametrize("ra, dec", positions)
def test_ang2ipix_xy_is_half_the_scale_of_ipix2xy(depth, ra, dec):
	'''
	Test that the two are related by exactly a factor of two.

	ipix2xy returns the lower left corner of the pixel rather than the point, so at a
	coarse depth the two differ by up to a pixel as well; the comparison is made at
	depth 30, where a pixel is ~1e-9 across, and the tolerance is one pixel wide.
	'''
	qlsc = QLSC(depth=depth)

	position = qlsc.ang2ipix_xy(ra=ra, dec=dec)
	face, x, y = qlsc.ipix2xy(position["ipix"])

	assert face == position["facenum"]

	pixel = 2. / qlsc.nside # the width of a pixel in [-1,1] face coordinates
	assert 2 * position["x"] == pytest.approx(x, abs=pixel)
	assert 2 * position["y"] == pytest.approx(y, abs=pixel)

@pytest.mark.parametrize("ra, dec", positions)
def test_doubling_is_what_makes_the_round_trip_work(ra, dec):
	'''
	Test the practical consequence: ang2ipix_xy's values must be doubled before
	xy2ang will return the coordinate they came from, and passing them undoubled
	returns some other point without any complaint.
	'''
	qlsc = QLSC(depth=30)
	position = qlsc.ang2ipix_xy(ra=ra, dec=dec)

	doubled = qlsc.xy2ang(facenum=position["facenum"], x=2 * position["x"], y=2 * position["y"])
	assert doubled[0] == pytest.approx(ra % 360., abs=1e-9)
	assert doubled[1] == pytest.approx(dec, abs=1e-9)

	# the undoubled values are accepted and answer for a different point on the sky
	undoubled = qlsc.xy2ang(facenum=position["facenum"], x=position["x"], y=position["y"])
	if not (position["x"] == 0. and position["y"] == 0.):
		separation = max(abs(undoubled[0] - doubled[0]), abs(undoubled[1] - doubled[1]))
		assert separation > 1e-6, "the undoubled values should name a different point, silently"

@pytest.mark.parametrize("ra, dec", positions)
def test_the_documented_recipe_round_trips_to_the_same_pixel(ra, dec):
	'''
	Test the snippet the docstring gives, end to end: doubling the values and feeding
	them through xy2ipix returns the pixel ang2ipix_xy reported.
	'''
	qlsc = QLSC(depth=20)
	position = qlsc.ang2ipix_xy(ra=ra, dec=dec)

	assert qlsc.xy2ipix(facenum=position["facenum"],
						x=2 * position["x"], y=2 * position["y"]) == position["ipix"]

def test_ang2ipix_xy_stays_within_the_half_scale_face():
	'''
	Test that ang2ipix_xy's values really do span [-0.5,0.5] and not more, over a
	spread of positions covering all six faces.
	'''
	qlsc = QLSC(depth=10)
	rng = np.random.default_rng(1975)

	extreme = 0.
	for _ in range(2000):
		ra = rng.uniform(0., 360.)
		dec = np.degrees(np.arcsin(rng.uniform(-1., 1.)))
		position = qlsc.ang2ipix_xy(ra=ra, dec=dec)
		extreme = max(extreme, abs(position["x"]), abs(position["y"]))

	assert extreme <= 0.5 + 1e-12
	assert extreme > 0.45, "the sample should reach the edges of a face"

def test_the_discrepancy_is_documented():
	'''
	Test that the warning is still in the docstring the API documentation is built
	from. This is the only thing standing between a user and a silently wrong
	coordinate, so it is not left to chance.
	'''
	documentation = QLSC.ang2ipix_xy.__doc__
	assert documentation is not None
	assert "[-0.5,0.5]" in documentation
	assert "[-1,1]" in documentation
	assert "factor of two" in documentation

	for method in (QLSC.xy2ang, QLSC.xy2ipix, QLSC.ipix2xy):
		assert "ang2ipix_xy" in method.__doc__, \
			f"{method.__name__} should point at the one method that uses the other scale"
