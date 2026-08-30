
import pytest
import numpy as np

from qlsc import QLSC, QLSCIndex, distance, sindist
from qlsc import q3c
from qlsc.generate import sunflower_points_on_sphere

# Right ascension is periodic in 360 degrees, so every entry point that takes one
# must answer identically for ra, ra+360k and ra-360k. Several vendored Q3C
# functions carry the documented precondition "strictly 0<=ra<360" and do nothing
# to enforce it - q3c_get_facenum computes the cube face as
# (char)fmod((ra+45)/90, 4) + 1, which for a negative ra truncates toward zero and
# returned -2 for ra=-315 - so the wrapper reduces ra before handing it over.
#
# This module enumerates the class rather than the one function that was found:
# every public API taking an ra is listed and pinned here, and so is every C
# function that reaches a face-selection routine. An entry point added later
# belongs in these lists.

q30 = QLSC(depth=30)

SHIFTS = [-1080, -720, -360, 360, 720, 1080]

# ra values chosen to land on every face and on the face boundaries at 45+90k
RA_VALUES = [0., 10., 44.9, 45., 45.1, 90., 134.9, 135., 180., 224.9, 225., 270.,
			 314.9, 315., 350., 359.9]
DEC_VALUES = [-89.9, -46., -44., 0., 20., 44., 46., 89.9]

# ra=0 is excluded from the ipix tests only: it sits exactly on a *pixel* boundary,
# where the two representations land on either side of it. That is pinned on its own
# in test_ang2ipix_at_the_ra_zero_pixel_boundary below, with the reason. It is not
# excluded from the face tests, because the face boundaries are at 45+90k.
RA_VALUES_OFF_PIXEL_BOUNDARY = [value for value in RA_VALUES if value != 0.]


def _polygon(centre_ra, half_width=3.):
	'''
	A small square polygon centered on the given right ascension.
	'''
	return np.array([[centre_ra - half_width, -half_width],
					 [centre_ra + half_width, -half_width],
					 [centre_ra + half_width, half_width],
					 [centre_ra - half_width, half_width]])


def _catalog_index():
	index = QLSCIndex(qlsc=q30)
	points = sunflower_points_on_sphere(n=4000)
	index.add_points(ra=points[:,0], dec=points[:,1])
	return index


# ------------------------------------------------------- the reported failure

@pytest.mark.parametrize("dec", DEC_VALUES)
@pytest.mark.parametrize("ra", RA_VALUES)
def test_face_number_is_periodic(ra, dec):
	'''
	Test that face_number gives the same answer for ra, ra+360k and ra-360k.

	Regression test: q3c_get_facenum requires 0<=ra<360 and the wrapper did not
	reduce it, so face_number(-315, 20) returned -2 - not even a face number - and
	face_number(-350, 20) returned 0 where ang2ipix_xy said 1.
	'''
	expected = q30.face_number(ra, dec)
	assert 0 <= expected <= 5

	for shift in SHIFTS:
		assert q30.face_number(ra + shift, dec) == expected


@pytest.mark.parametrize("dec", DEC_VALUES)
@pytest.mark.parametrize("ra", RA_VALUES)
def test_face_number_agrees_with_ang2ipix_xy(ra, dec):
	'''
	Test that face_number and ang2ipix_xy name the same face for the same point,
	at every representation of the right ascension. These two disagreed for any
	negative ra: ang2ipix_xy reaches q3c_ang2ipix_xy, which unwraps internally.
	'''
	expected = q30.ang2ipix_xy(ra, dec)["facenum"]

	for shift in [0] + SHIFTS:
		assert q30.face_number(ra + shift, dec) == expected


def test_face_number_no_longer_rejects_out_of_range_ra():
	'''
	Test that a wrapped right ascension is accepted rather than rejected.

	face_number used to raise for ra outside [-360,360], which hid the fact that
	everything in [-360,0) was quietly answered wrongly.
	'''
	assert q30.face_number(-350., 20.) == q30.face_number(10., 20.)
	assert q30.face_number(1000., 20.) == q30.face_number(1000. % 360., 20.)

	# dec is not periodic and is still rejected
	with pytest.raises(ValueError):
		q30.face_number(10., 91.)


# ------------------------------------------------------- the rest of the class

@pytest.mark.parametrize("ra", RA_VALUES_OFF_PIXEL_BOUNDARY)
def test_ang2ipix_is_periodic(ra):
	'''
	Test ang2ipix, scalar and array form (it was already correct: the vendored
	q3c_ang2ipix_xy unwraps on entry).
	'''
	for dec in DEC_VALUES:
		expected = q30.ang2ipix(ra, dec)
		for shift in SHIFTS:
			assert q30.ang2ipix(ra + shift, dec) == expected

	shifted = np.array([ra + shift for shift in [0] + SHIFTS])
	dec_array = np.full(shifted.shape, 20.)
	assert len(set(q30.ang2ipix(ra=shifted, dec=dec_array).tolist())) == 1


@pytest.mark.parametrize("ra", RA_VALUES)
def test_ang2ipix_xy_is_periodic(ra):
	'''
	Test ang2ipix_xy, which returns the face, the ipix and the position on the face.

	The face and the ipix must be identical; the position is compared to a
	tolerance, because reducing ra by whole turns is not exact in floating point
	((44.9 - 360) + 360 is 44.900000000000034) and the projection carries those
	last few bits through. The difference is ~1e-16 of a coordinate that runs
	-1..1, a ten-millionth of a pixel at depth 30.
	'''
	expected = q30.ang2ipix_xy(ra, 20.)

	for shift in SHIFTS:
		result = q30.ang2ipix_xy(ra + shift, 20.)
		assert result["facenum"] == expected["facenum"]
		if ra != 0.: # see test_ang2ipix_at_the_ra_zero_pixel_boundary
			assert result["ipix"] == expected["ipix"]
		assert result["x"] == pytest.approx(expected["x"], abs=1e-12)
		assert result["y"] == pytest.approx(expected["y"], abs=1e-12)


@pytest.mark.parametrize("ra", RA_VALUES_OFF_PIXEL_BOUNDARY)
def test_ang2ipix_at_depth_is_periodic(ra):
	'''
	Test ang2ipix_at_depth, which goes through the same conversion.
	'''
	expected = q30.ang2ipix_at_depth(ra, 20., 10)
	for shift in SHIFTS:
		assert q30.ang2ipix_at_depth(ra + shift, 20., 10) == expected


@pytest.mark.parametrize("depth", [1, 2, 10, 30])
def test_ang2ipix_at_the_ra_zero_pixel_boundary(depth):
	'''
	Pin the one place a wrapped right ascension does not give the identical ipix.

	ra=0 lies exactly on a pixel boundary (x = tan(0) = 0 is a bin edge whenever
	nside is even). Q3C's UNWRAP_RA macro maps every negative whole turn to exactly
	360.0 rather than 0.0 — fmod(-360, 360) is -0.0, and -0.0 + 360 is 360.0 — and
	tan(2*pi) evaluates to -2.4e-16 rather than 0, so ra=-360 lands a hair to the
	*left* of the boundary and ra=0 lands on it.

	The two pixels are adjacent: the point is on their shared edge and either answer
	describes it. This is upstream behavior, unchanged here — the macro lives in the
	vendored q3c_ang2ipix_xy, and changing it would make qlsc disagree with the Q3C
	PostgreSQL plugin about the ipix of the same coordinate. It is pinned rather
	than papered over so that a future change to it is noticed.
	'''
	q = QLSC(depth=depth)

	at_zero = q.ang2ipix(0., 20.)
	at_360 = q.ang2ipix(360., 20.)

	assert at_zero != at_360
	assert at_360 in q.ipix_neighbors(at_zero), "the two must at least be adjacent pixels"

	# every negative whole turn agrees with 360, not with 0
	for shift in [-360., -720., -1080.]:
		assert q.ang2ipix(shift, 20.) == at_360

	# and away from the boundary, whole turns are exact
	assert q.ang2ipix(-350., 20.) == q.ang2ipix(10., 20.)


def test_ang2ipix_at_ra_zero_is_stable_at_depth_zero():
	'''
	Test the counterpart: at depth 0 a face is a single pixel, so ra=0 is not on a
	pixel boundary and every representation gives the same answer.
	'''
	q = QLSC(depth=0)
	expected = q.ang2ipix(0., 20.)
	for ra in (360., -360., -1080., 720.):
		assert q.ang2ipix(ra, 20.) == expected


@pytest.mark.parametrize("ra", RA_VALUES)
def test_distance_and_sindist_are_periodic(ra):
	'''
	Test the two separation functions. ra enters them only as sin((ra1-ra2)/2)**2,
	which is periodic in 360 for each argument separately.
	'''
	expected_distance = distance(ra, 20., 100., -30.)
	expected_sindist = sindist(ra, 20., 100., -30.)

	for shift in SHIFTS:
		assert distance(ra + shift, 20., 100., -30.) == pytest.approx(expected_distance, abs=1e-9)
		assert sindist(ra + shift, 20., 100., -30.) == pytest.approx(expected_sindist, abs=1e-12)
		assert distance(ra, 20., 100. + shift, -30.) == pytest.approx(expected_distance, abs=1e-9)


@pytest.mark.parametrize("centre", [10., 45., 100., 180., 265., 310., 350.])
def test_point_in_polygon_is_periodic_in_the_vertices(centre):
	'''
	Test that a polygon written with wrapped vertex right ascensions - the natural
	way to describe one straddling ra=0 - behaves like the plain one.

	Regression test: both polygon entry points reach q3c_get_facenum through
	q3c_get_facenum_poly, which reads the FIRST vertex, so a polygon at ra=-140
	picked the wrong cube face and was then rejected as "too large".
	'''
	plain = _polygon(centre)

	for shift in SHIFTS:
		wrapped = _polygon(centre + shift)
		for offset in (0., 2., -2., 5.):
			ra = (centre + offset) % 360.
			assert q30.point_in_polygon(ra=ra, dec=0., polygon=wrapped) == \
				   q30.point_in_polygon(ra=ra, dec=0., polygon=plain)


@pytest.mark.parametrize("centre", [10., 45., 100., 180., 265., 310., 350.])
def test_point_in_polygon_is_periodic_in_the_test_point(centre):
	'''
	Test the same for the point being tested, rather than the polygon.
	'''
	polygon = _polygon(centre)

	for shift in [0] + SHIFTS:
		assert q30.point_in_polygon(ra=centre + shift, dec=0., polygon=polygon) is True
		assert q30.point_in_polygon(ra=centre + 10. + shift, dec=0., polygon=polygon) is False


@pytest.mark.parametrize("centre", [10., 100., 180., 265., 350.])
def test_point_in_ellipse_is_periodic(centre):
	'''
	Test point_in_ellipse in both its right ascensions. q3c_in_ellipse uses them
	only through cos and sin of (alpha1 - alpha), so it is periodic already; this
	pins it.
	'''
	def call(ra, ra_cen):
		return q30.point_in_ellipse(ra=ra, dec=0., ra_cen=ra_cen, dec_cen=0.,
									semi_major_axis=5., axis_ratio=0.5, position_angle=30.)

	for shift in SHIFTS:
		assert call(centre + 2. + shift, centre) == call(centre + 2., centre)
		assert call(centre + 2., centre + shift) == call(centre + 2., centre)
		assert call(centre + 20. + shift, centre) == call(centre + 20., centre)


@pytest.mark.parametrize("centre", [10., 45., 100., 180., 265., 310., 350.])
def test_index_queries_are_periodic(centre):
	'''
	Test the three region queries on a populated index: the same region described
	with a wrapped right ascension must return the same points.
	'''
	with _catalog_index() as index:
		radial = index.radial_query(centre, 15., 8.)
		polygon = index.polygon_query(polygon=_polygon(centre, half_width=6.))
		ellipse = index.ellipse_query(ra=centre, dec=15., semi_major_axis=8.,
									  axis_ratio=0.6, position_angle=20.)

		assert radial.size > 0 and ellipse.size > 0, "the comparison would be vacuous"

		for shift in SHIFTS:
			assert np.array_equal(index.radial_query(centre + shift, 15., 8.), radial)
			assert np.array_equal(index.polygon_query(polygon=_polygon(centre + shift, half_width=6.)), polygon)
			assert np.array_equal(index.ellipse_query(ra=centre + shift, dec=15., semi_major_axis=8.,
													  axis_ratio=0.6, position_angle=20.), ellipse)


def test_points_added_with_a_wrapped_ra_are_found():
	'''
	Test that a point stored with a wrapped right ascension is found by a query
	written the ordinary way, and the other way round.
	'''
	with QLSCIndex(qlsc=q30) as index:
		index.add_point(-10., 30., key="negative")
		index.add_points(ra=np.array([730.]), dec=np.array([-40.]), keys=["large"])

		assert len(index.radial_query(350., 30., 0.5)) > 0
		assert len(index.radial_query(10., -40., 0.5)) > 0
		assert len(index.radial_query(-10., 30., 0.5)) > 0


# ------------------------------------------------------- the C boundary itself

@pytest.mark.parametrize("ra", RA_VALUES)
def test_c_level_facenum_is_periodic(ra):
	'''
	Test the exported C function directly: qlsc.q3c is importable on its own, so
	the boundary has to hold there and not only through QLSC.
	'''
	expected = q3c.facenum(q30._hprm, ra, 20.)
	assert 0 <= expected <= 5

	for shift in SHIFTS:
		assert q3c.facenum(q30._hprm, ra + shift, 20.) == expected


def test_c_level_facenum_rejects_non_finite():
	'''
	Test that non-finite values are refused rather than reaching the undefined
	float-to-int conversion in q3c_get_facenum's (char) cast (upstream issue #52).
	'''
	for bad in (float("nan"), float("inf"), -float("inf")):
		with pytest.raises(ValueError):
			q3c.facenum(q30._hprm, bad, 20.)
		with pytest.raises(ValueError):
			q3c.facenum(q30._hprm, 20., bad)


@pytest.mark.parametrize("centre", [10., 100., 265., 350.])
def test_c_level_polygon_functions_are_periodic(centre):
	'''
	Test the two exported polygon functions directly.
	'''
	plain = _polygon(centre)

	for shift in SHIFTS:
		wrapped = _polygon(centre + shift)

		assert q3c.check_sphere_point_in_poly(q30._hprm, centre, 0.,
											  np.ascontiguousarray(wrapped[:,0]),
											  np.ascontiguousarray(wrapped[:,1])) == \
			   q3c.check_sphere_point_in_poly(q30._hprm, centre, 0.,
											  np.ascontiguousarray(plain[:,0]),
											  np.ascontiguousarray(plain[:,1]))

		wrapped_fulls, wrapped_partials = q3c.poly_query(q30._hprm,
														 np.ascontiguousarray(wrapped[:,0]),
														 np.ascontiguousarray(wrapped[:,1]))
		plain_fulls, plain_partials = q3c.poly_query(q30._hprm,
													 np.ascontiguousarray(plain[:,0]),
													 np.ascontiguousarray(plain[:,1]))
		assert np.array_equal(wrapped_fulls, plain_fulls)
		assert np.array_equal(wrapped_partials, plain_partials)


def test_c_level_polygon_functions_do_not_modify_the_caller_array():
	'''
	Test that unwrapping the vertex right ascensions leaves the caller's array
	alone - the values are copied, not reduced in place.
	'''
	polygon = _polygon(-140.)
	poly_ra = np.ascontiguousarray(polygon[:,0])
	poly_dec = np.ascontiguousarray(polygon[:,1])
	original = poly_ra.copy()

	q3c.check_sphere_point_in_poly(q30._hprm, 220., 0., poly_ra, poly_dec)
	q3c.poly_query(q30._hprm, poly_ra, poly_dec)

	assert np.array_equal(poly_ra, original)


def test_c_level_polygon_functions_reject_non_finite_vertices():
	'''
	Test that a non-finite vertex is refused at the C boundary as well as in the
	Python layer.
	'''
	poly_ra = np.array([np.nan, 5., -5.])
	poly_dec = np.array([-5., 5., 5.])

	with pytest.raises(ValueError):
		q3c.check_sphere_point_in_poly(q30._hprm, 0., 0., poly_ra, poly_dec)
	with pytest.raises(ValueError):
		q3c.poly_query(q30._hprm, poly_ra, poly_dec)
