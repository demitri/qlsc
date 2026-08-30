
import pytest
import numpy as np

from qlsc import QLSC, QLSCIndex

# Tests for QLSC.point_in_ellipse, the exact containment test that filters
# QLSCIndex.ellipse_query results. The underlying q3c_in_ellipse takes the
# eccentricity; the wrapper takes the axis ratio, as q3c's own SQL interface does.

q30 = QLSC(depth=30)

def test_center_is_inside():
	'''
	Test that the center of an ellipse is inside it.
	'''
	assert q30.point_in_ellipse(ra=10., dec=20., ra_cen=10., dec_cen=20.,
								semi_major_axis=1., axis_ratio=0.5, position_angle=30.) is True

def test_axis_ratio_one_is_a_circle():
	'''
	Test that an axis ratio of 1 gives a circle: containment depends only on the
	angular separation, whatever the position angle.
	'''
	for position_angle in (0., 45., 90., 217.):
		# 2 degrees due north of the center
		assert q30.point_in_ellipse(ra=10., dec=22., ra_cen=10., dec_cen=20.,
									semi_major_axis=3., axis_ratio=1., position_angle=position_angle) is True
		assert q30.point_in_ellipse(ra=10., dec=24.5, ra_cen=10., dec_cen=20.,
									semi_major_axis=3., axis_ratio=1., position_angle=position_angle) is False

def test_position_angle_orients_the_major_axis():
	'''
	Test that the position angle rotates the ellipse: a point that is inside a
	north-south elongated ellipse is outside the same ellipse turned east-west.
	'''
	north = dict(ra=0., dec=3., ra_cen=0., dec_cen=0., semi_major_axis=4., axis_ratio=0.2)
	assert q30.point_in_ellipse(position_angle=0., **north) is True
	assert q30.point_in_ellipse(position_angle=90., **north) is False

	east = dict(ra=3., dec=0., ra_cen=0., dec_cen=0., semi_major_axis=4., axis_ratio=0.2)
	assert q30.point_in_ellipse(position_angle=90., **east) is True
	assert q30.point_in_ellipse(position_angle=0., **east) is False

def test_semi_minor_axis_bounds():
	'''
	Test the extent of the ellipse across its minor axis: with a semi-major axis of
	4 degrees and an axis ratio of 0.25 the semi-minor axis is 1 degree.
	'''
	inside = dict(ra=0.9, dec=0., ra_cen=0., dec_cen=0., semi_major_axis=4., axis_ratio=0.25, position_angle=0.)
	outside = dict(inside, ra=1.1)
	assert q30.point_in_ellipse(**inside) is True
	assert q30.point_in_ellipse(**outside) is False

def test_antipodal_point_is_outside():
	'''
	Test that the point opposite the center is not inside; q3c_in_ellipse rejects
	the far hemisphere before evaluating the conic.
	'''
	assert q30.point_in_ellipse(ra=190., dec=-20., ra_cen=10., dec_cen=20.,
								semi_major_axis=30., axis_ratio=0.9) is False

def test_equivalent_coordinate_representations():
	'''
	Test that ra=360 and ra=0, and any ra at a pole, describe the same point.
	'''
	assert q30.point_in_ellipse(ra=360., dec=20., ra_cen=0., dec_cen=20.,
								semi_major_axis=0.001, axis_ratio=1.) is True
	assert q30.point_in_ellipse(ra=123., dec=90., ra_cen=0., dec_cen=90.,
								semi_major_axis=0.001, axis_ratio=1.) is True

def test_missing_parameters_raise():
	'''
	Test that omitting any required parameter raises rather than being silently defaulted.
	'''
	full = dict(ra=1., dec=2., ra_cen=1., dec_cen=2., semi_major_axis=1., axis_ratio=0.5)
	for name in full:
		partial = dict(full)
		partial[name] = None
		with pytest.raises(ValueError):
			q30.point_in_ellipse(**partial)

# parameter, value
invalid_parameter_cases = [
	("axis_ratio", 0.),
	("axis_ratio", -0.5),
	("axis_ratio", 1.5),
	("axis_ratio", np.nan),
	("semi_major_axis", -1.),
	("semi_major_axis", np.inf),
	("position_angle", np.nan),
	("ra", np.nan),
	("dec", np.inf),
	("dec_cen", np.nan)
]

@pytest.mark.parametrize("parameter, value", invalid_parameter_cases)
def test_invalid_parameters_raise(parameter, value):
	'''
	Test that out of range and non-finite parameters raise ValueError rather than
	reaching the C code, where an axis ratio above 1 would become a NaN eccentricity
	and a non-finite coordinate an undefined float to int conversion.
	'''
	arguments = dict(ra=1., dec=2., ra_cen=1., dec_cen=2., semi_major_axis=1., axis_ratio=0.5, position_angle=0.)
	arguments[parameter] = value
	with pytest.raises(ValueError):
		q30.point_in_ellipse(**arguments)

def test_out_of_range_dec_is_normalized():
	'''
	Test that a dec beyond the pole is normalized rather than rejected, as
	point_in_polygon does: dec=91 is dec=89 on the opposite side of the sphere.
	'''
	assert q30.point_in_ellipse(ra=0., dec=91., ra_cen=180., dec_cen=89.,
								semi_major_axis=1., axis_ratio=1.) is True

def test_result_does_not_depend_on_depth():
	'''
	Test that the containment test is a property of the sphere, not of the scheme.
	'''
	arguments = dict(ra=10.5, dec=20.3, ra_cen=10., dec_cen=20., semi_major_axis=1., axis_ratio=0.4, position_angle=17.)
	assert QLSC(depth=0).point_in_ellipse(**arguments) == QLSC(depth=30).point_in_ellipse(**arguments)

def test_astropy_quantity_parameters():
	'''
	Test that the semi-major axis and position angle accept astropy Quantities.
	'''
	astropy_units = pytest.importorskip("astropy.units")
	degrees = q30.point_in_ellipse(ra=0., dec=1., ra_cen=0., dec_cen=0.,
								   semi_major_axis=2., axis_ratio=0.3, position_angle=0.)
	quantities = q30.point_in_ellipse(ra=0., dec=1., ra_cen=0., dec_cen=0.,
									  semi_major_axis=7200 * astropy_units.arcsec,
									  axis_ratio=0.3,
									  position_angle=0 * astropy_units.rad)
	assert degrees is True
	assert quantities == degrees

# ---------------------------------------------------------------------------
# equivalent representations of the center
#
# ra=0 and ra=360 name the same meridian, and at a pole the right ascension is
# arbitrary, so a point given either way IS the center. The angular separation
# between two such representations is a floating point residual rather than
# exactly zero, so a raw comparison against the radius excluded the center from
# its own circle whenever the radius was smaller than that residual - which the
# degenerate circle of radius 0 always is.
# ---------------------------------------------------------------------------

EQUIVALENT_CENTERS = [
	# (ra, dec, ra_cen, dec_cen, what it is)
	(360., 20., 0., 20., "ra=360 against ra=0"),
	(0., 20., 360., 20., "ra=0 against ra=360"),
	(-350., 20., 10., 20., "a negative ra against its positive form"),
	(137., 90., 0., 90., "two right ascensions at the north pole"),
	(0., -90., 212., -90., "two right ascensions at the south pole"),
	(10., 20., 10., 20., "the identical coordinate"),
]

@pytest.mark.parametrize("ra,dec,ra_cen,dec_cen,description", EQUIVALENT_CENTERS)
@pytest.mark.parametrize("semi_major_axis", [0., 1e-12, 1e-6, 1.])
def test_a_circle_contains_every_representation_of_its_center(ra, dec, ra_cen, dec_cen,
															  description, semi_major_axis):
	'''
	Test that a point which is another representation of the center is inside a circle
	of any size, the degenerate circle of radius 0 included.
	'''
	assert q30.point_in_ellipse(ra=ra, dec=dec, ra_cen=ra_cen, dec_cen=dec_cen,
								semi_major_axis=semi_major_axis, axis_ratio=1.) is True, description

@pytest.mark.parametrize("ra,dec,ra_cen,dec_cen,description", EQUIVALENT_CENTERS)
@pytest.mark.parametrize("semi_major_axis", [0., 1e-12, 1.])
def test_the_circle_agrees_with_the_cone_search(ra, dec, ra_cen, dec_cen,
												description, semi_major_axis):
	'''
	Test the docstring's claim that a circular point_in_ellipse and radial_query apply
	the identical test — on the coordinates where they used to disagree.

	This is the promise that matters: ellipse_query filters its candidates with the
	containment test, so a point the index reports and the predicate denies would be the
	two halves of one API contradicting each other.
	'''
	index = QLSCIndex(qlsc=q30)
	index.add_point(ra=ra, dec=dec)

	found = len(np.atleast_2d(index.radial_query(ra=ra_cen, dec=dec_cen,
												 radius=semi_major_axis))) > 0
	predicate = q30.point_in_ellipse(ra=ra, dec=dec, ra_cen=ra_cen, dec_cen=dec_cen,
									 semi_major_axis=semi_major_axis, axis_ratio=1.)

	assert predicate is found, description

	# ellipse_query hands the circular case to radial_query, so it must agree too
	from_ellipse = len(index.ellipse_query(ra=ra_cen, dec=dec_cen,
										   semi_major_axis=semi_major_axis, axis_ratio=1.)) > 0
	assert from_ellipse is found, description

@pytest.mark.parametrize("ra,dec,ra_cen,dec_cen,description", EQUIVALENT_CENTERS)
@pytest.mark.parametrize("semi_major_axis", [0., 1e-12, 1.])
def test_a_genuine_ellipse_answers_the_same_for_every_representation(ra, dec, ra_cen, dec_cen,
																	 description, semi_major_axis):
	'''
	Test the rule on the Q3C conic path: the answer is a property of the point, not of how
	its coordinates were written.

	The claim here is deliberately weaker than "the center is inside". Q3C's convention is
	that the edge of a genuinely elliptical region is outside it, and the center of a
	zero-size ellipse IS its edge — at a pole, where the tangent plane projection is
	degenerate, Q3C answers False for the center given identically. That convention is
	kept; what is fixed is that writing the same point a different way used to change the
	answer.
	'''
	expected = q30.point_in_ellipse(ra=ra_cen, dec=dec_cen, ra_cen=ra_cen, dec_cen=dec_cen,
									semi_major_axis=semi_major_axis, axis_ratio=0.5)

	assert q30.point_in_ellipse(ra=ra, dec=dec, ra_cen=ra_cen, dec_cen=dec_cen,
								semi_major_axis=semi_major_axis, axis_ratio=0.5) is expected, description

@pytest.mark.parametrize("ra,dec,ra_cen,dec_cen,description", EQUIVALENT_CENTERS)
@pytest.mark.parametrize("semi_major_axis", [0., 1e-12, 1.])
def test_the_ellipse_query_agrees_with_the_predicate_it_documents(ra, dec, ra_cen, dec_cen,
																  description, semi_major_axis):
	'''
	Test that ellipse_query returns exactly the points point_in_ellipse calls inside, on
	the coordinates where the two used to contradict each other.

	ellipse_query documents that it filters its candidates with the containment test, so a
	point the index reports and the predicate denies is the two halves of one API
	disagreeing. It did: the query matched any representation of its center outright while
	the predicate ran the projection.
	'''
	index = QLSCIndex(qlsc=q30)
	index.add_point(ra=ra, dec=dec)

	predicate = q30.point_in_ellipse(ra=ra, dec=dec, ra_cen=ra_cen, dec_cen=dec_cen,
									 semi_major_axis=semi_major_axis, axis_ratio=0.5)
	found = len(index.ellipse_query(ra=ra_cen, dec=dec_cen,
									semi_major_axis=semi_major_axis, axis_ratio=0.5)) > 0

	assert found is predicate, description

def test_a_point_that_is_not_the_center_is_still_outside_a_zero_radius_circle():
	'''
	Test that the canonicalization did not turn the degenerate circle into one that
	contains everything: a genuinely different point is still outside.
	'''
	assert q30.point_in_ellipse(ra=10.0001, dec=20., ra_cen=10., dec_cen=20.,
								semi_major_axis=0., axis_ratio=1.) is False
	# a pole is one point, but the other pole is not
	assert q30.point_in_ellipse(ra=0., dec=-90., ra_cen=0., dec_cen=90.,
								semi_major_axis=0., axis_ratio=1.) is False
	# and the same meridian at a different declination is a different place
	assert q30.point_in_ellipse(ra=360., dec=20.5, ra_cen=0., dec_cen=20.,
								semi_major_axis=0., axis_ratio=1.) is False
