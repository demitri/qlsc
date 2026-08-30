
import numpy as np
import pytest

import qlsc.qlsc_base
from qlsc import QLSC, QLSCIndex

# Tests for first-class astropy support: a SkyCoord (scalar or array, in any frame)
# accepted wherever an ra,dec pair goes, and an opt-in SkyCoord return from the methods
# that hand coordinates back.
#
# astropy is a soft dependency, so every test that needs it skips loudly via
# importorskip, the pattern the rest of the suite uses. The one test that exercises the
# absence of astropy fakes it by flipping the module's availability flag, which is the
# only way to reach that branch on a machine that has astropy installed.

astropy_coordinates = pytest.importorskip("astropy.coordinates")
astropy_units = pytest.importorskip("astropy.units")
astropy_table = pytest.importorskip("astropy.table")

SkyCoord = astropy_coordinates.SkyCoord
u = astropy_units
QTable = astropy_table.QTable

q30 = QLSC(depth=30)

# a position with no special symmetry, in ICRS degrees (roughly M31)
M31_RA = 10.6847083
M31_DEC = 41.26875

def index_with_points():
	'''
	Build a small in-memory index with keyed points, one of them at M31.
	'''
	index = QLSCIndex(qlsc=QLSC(depth=30))
	index.add_point(ra=M31_RA, dec=M31_DEC, key="m31")
	index.add_point(ra=M31_RA + 0.01, dec=M31_DEC + 0.01, key="near")
	index.add_point(ra=200., dec=-30., key="far")
	return index

# ---------------------------------------------------------------- input: QLSC methods

def test_ang2ipix_scalar_skycoord_matches_degrees():
	'''
	Test that a scalar SkyCoord given as the single positional argument, or as 'ra',
	gives the same ipix as the equivalent degrees.
	'''
	expected = q30.ang2ipix(M31_RA, M31_DEC)
	coord = SkyCoord(ra=M31_RA * u.deg, dec=M31_DEC * u.deg, frame="icrs")

	assert q30.ang2ipix(coord) == expected
	assert q30.ang2ipix(ra=coord) == expected

def test_ang2ipix_array_skycoord_matches_degrees():
	'''
	Test that an array SkyCoord gives the same ipix array as the equivalent degrees,
	whether passed as 'ra' or as 'points'.
	'''
	ra = np.array([M31_RA, 200., 359.9, 0.0])
	dec = np.array([M31_DEC, -30., 0.0, 89.9])
	expected = q30.ang2ipix(ra=ra, dec=dec)

	coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")

	assert np.array_equal(q30.ang2ipix(coord), expected)
	assert np.array_equal(q30.ang2ipix(points=coord), expected)

def test_ang2ipix_length_one_skycoord_returns_an_array():
	'''
	Test that a length-1 array SkyCoord returns an array of one ipix, not a bare
	integer - the same scalar/array convention the degrees interface follows.
	'''
	coord = SkyCoord(ra=[M31_RA] * u.deg, dec=[M31_DEC] * u.deg, frame="icrs")
	result = q30.ang2ipix(coord)

	assert isinstance(result, np.ndarray)
	assert result.shape == (1,)
	assert result[0] == q30.ang2ipix(M31_RA, M31_DEC)

def test_galactic_skycoord_is_transformed_not_reinterpreted():
	'''
	Test that a SkyCoord in another frame is converted to ICRS: the galactic
	representation of a known ICRS position yields the ipix of that position, and
	NOT the ipix of its galactic angles read as ra,dec.
	'''
	expected = q30.ang2ipix(M31_RA, M31_DEC)

	galactic = SkyCoord(ra=M31_RA * u.deg, dec=M31_DEC * u.deg, frame="icrs").galactic

	assert q30.ang2ipix(galactic) == expected

	# the check that the transformation actually happened: taking the galactic angles
	# at face value gives a different pixel entirely
	naive = q30.ang2ipix(galactic.l.deg, galactic.b.deg)
	assert naive != expected

def test_fk5_and_frame_objects_are_accepted():
	'''
	Test that an FK5 SkyCoord and a bare coordinate frame object both work and agree
	with the ICRS position they represent.
	'''
	expected = q30.ang2ipix(M31_RA, M31_DEC)
	icrs = SkyCoord(ra=M31_RA * u.deg, dec=M31_DEC * u.deg, frame="icrs")

	assert q30.ang2ipix(icrs.fk5) == expected
	# a BaseCoordinateFrame instance rather than a SkyCoord
	assert q30.ang2ipix(icrs.galactic.frame) == expected

def test_quantity_angles_are_converted():
	'''
	Test that individual angles may be given as Quantities in any angular unit.
	'''
	expected = q30.ang2ipix(45., 60.)

	assert q30.ang2ipix(45. * u.deg, 60. * u.deg) == expected
	assert q30.ang2ipix(np.deg2rad(45.) * u.rad, 3600. * 60. * u.arcsec) == expected

	ra = np.array([45., 200.])
	dec = np.array([60., -30.])
	assert np.array_equal(q30.ang2ipix(ra * u.deg, dec * u.deg), q30.ang2ipix(ra, dec))

def test_skycoord_in_the_dec_slot_raises():
	'''
	Test that a SkyCoord given as the second of the two angles is refused, rather than
	being coerced into a meaningless number.
	'''
	coord = SkyCoord(ra=M31_RA * u.deg, dec=M31_DEC * u.deg, frame="icrs")

	with pytest.raises(ValueError, match="passed as 'ra'"):
		q30.ang2ipix(10., coord)

def test_skycoord_with_a_second_angle_raises():
	'''
	Test that giving both a SkyCoord and a declination is an error: the SkyCoord already
	carries both angles, so one of the two values would otherwise be silently discarded.
	'''
	coord = SkyCoord(ra=M31_RA * u.deg, dec=M31_DEC * u.deg, frame="icrs")

	with pytest.raises(ValueError, match="already carries both angles"):
		q30.ang2ipix(coord, 20.)

def test_skycoord_and_points_together_raise():
	'''
	Test that a SkyCoord and a 'points' array cannot both be given.
	'''
	coord = SkyCoord(ra=[M31_RA] * u.deg, dec=[M31_DEC] * u.deg, frame="icrs")

	with pytest.raises(ValueError, match="should not also be specified"):
		q30.ang2ipix(ra=coord, points=np.array([[1., 2.]]))

def test_ang2ipix_xy_accepts_skycoord():
	'''
	Test that ang2ipix_xy takes a scalar SkyCoord.
	'''
	coord = SkyCoord(ra=M31_RA * u.deg, dec=M31_DEC * u.deg, frame="icrs")

	assert q30.ang2ipix_xy(coord) == q30.ang2ipix_xy(ra=M31_RA, dec=M31_DEC)

def test_face_number_accepts_skycoord():
	'''
	Test that face_number takes a scalar SkyCoord, in any frame.
	'''
	coord = SkyCoord(ra=M31_RA * u.deg, dec=M31_DEC * u.deg, frame="icrs")

	assert q30.face_number(coord) == q30.face_number(M31_RA, M31_DEC)
	assert q30.face_number(coord.galactic) == q30.face_number(M31_RA, M31_DEC)

def test_ang2ipix_at_depth_accepts_skycoord():
	'''
	Test that ang2ipix_at_depth takes a scalar SkyCoord.
	'''
	coord = SkyCoord(ra=M31_RA * u.deg, dec=M31_DEC * u.deg, frame="icrs")

	assert q30.ang2ipix_at_depth(coord, depth=10) == q30.ang2ipix_at_depth(M31_RA, M31_DEC, depth=10)

def test_point_in_polygon_accepts_skycoords():
	'''
	Test that both the point and the polygon may be given as SkyCoords.
	'''
	ipix = q30.ang2ipix(M31_RA, M31_DEC)
	polygon = q30.ipix2polygon(ipix)
	polygon_coord = q30.ipix2polygon(ipix, return_skycoord=True)
	point = SkyCoord(ra=M31_RA * u.deg, dec=M31_DEC * u.deg, frame="icrs")

	assert q30.point_in_polygon(ra=M31_RA, dec=M31_DEC, polygon=polygon) is True
	assert q30.point_in_polygon(point, polygon=polygon) is True
	assert q30.point_in_polygon(point, polygon=polygon_coord) is True
	assert q30.point_in_polygon(point.galactic, polygon=polygon_coord) is True

def test_point_in_ellipse_accepts_skycoords():
	'''
	Test that the point and the ellipse center may each be given as a SkyCoord.
	'''
	point = SkyCoord(ra=10.5 * u.deg, dec=20. * u.deg, frame="icrs")
	center = SkyCoord(ra=10. * u.deg, dec=20. * u.deg, frame="icrs")

	expected = q30.point_in_ellipse(ra=10.5, dec=20., ra_cen=10., dec_cen=20.,
									semi_major_axis=1., axis_ratio=0.5, position_angle=90.)

	assert q30.point_in_ellipse(point, ra_cen=center, semi_major_axis=1.,
								axis_ratio=0.5, position_angle=90.) is expected
	assert q30.point_in_ellipse(point.galactic, ra_cen=center.fk5, semi_major_axis=1.,
								axis_ratio=0.5, position_angle=90.) is expected

# --------------------------------------------------------------- output: return_skycoord

def test_ipix2ang_return_skycoord_scalar():
	'''
	Test that ipix2ang returns a scalar ICRS SkyCoord agreeing with the degrees it
	otherwise returns, and that the default return value is unchanged.
	'''
	ipix = q30.ang2ipix(M31_RA, M31_DEC)
	ra, dec = q30.ipix2ang(ipix)

	coord = q30.ipix2ang(ipix, return_skycoord=True)

	assert isinstance(coord, SkyCoord)
	assert coord.isscalar
	assert coord.frame.name == "icrs"
	assert coord.ra.deg == pytest.approx(ra)
	assert coord.dec.deg == pytest.approx(dec)

	# the default is untouched
	assert isinstance(q30.ipix2ang(ipix), tuple)

def test_ipix2ang_return_skycoord_array():
	'''
	Test that array input gives an array SkyCoord of the same length.
	'''
	ipix = q30.ang2ipix(np.array([M31_RA, 200., 0.]), np.array([M31_DEC, -30., 89.]))
	ra, dec = q30.ipix2ang(ipix)

	coord = q30.ipix2ang(ipix, return_skycoord=True)

	assert not coord.isscalar
	assert len(coord) == 3
	assert np.allclose(coord.ra.deg, ra)
	assert np.allclose(coord.dec.deg, dec)

def test_ipix2ang_center_return_skycoord():
	'''
	Test that the ipix2ang_center synonym passes the flag through.
	'''
	ipix = q30.ang2ipix(M31_RA, M31_DEC)
	coord = q30.ipix2ang_center(ipix, return_skycoord=True)

	assert isinstance(coord, SkyCoord)
	assert coord.ra.deg == pytest.approx(q30.ipix2ang(ipix)[0])

def test_xy2ang_return_skycoord():
	'''
	Test that xy2ang returns a SkyCoord for both a single (x,y) and an array of them.
	'''
	scalar = q30.xy2ang(facenum=0, x=0.1, y=0.2, return_skycoord=True)
	ra, dec = q30.xy2ang(facenum=0, x=0.1, y=0.2)

	assert isinstance(scalar, SkyCoord)
	assert scalar.isscalar
	assert scalar.ra.deg == pytest.approx(ra)
	assert scalar.dec.deg == pytest.approx(dec)

	points = np.array([[0.1, 0.2], [-0.5, 0.5]])
	array_coord = q30.xy2ang(facenum=0, points=points, return_skycoord=True)
	array_degrees = q30.xy2ang(facenum=0, points=points)

	assert len(array_coord) == 2
	assert np.allclose(array_coord.ra.deg, array_degrees[:, 0])
	assert np.allclose(array_coord.dec.deg, array_degrees[:, 1])

def test_ipix2polygon_return_skycoord():
	'''
	Test that ipix2polygon returns the corners as a SkyCoord, in the same order, with
	duplicate_endpoint honored.
	'''
	ipix = q30.ang2ipix(M31_RA, M31_DEC)
	polygon = q30.ipix2polygon(ipix)
	coord = q30.ipix2polygon(ipix, return_skycoord=True)

	assert len(coord) == 4
	assert np.allclose(coord.ra.deg, polygon[:, 0])
	assert np.allclose(coord.dec.deg, polygon[:, 1])

	closed = q30.ipix2polygon(ipix, duplicate_endpoint=True, return_skycoord=True)
	assert len(closed) == 5
	assert closed.ra.deg[0] == pytest.approx(closed.ra.deg[-1])

# ------------------------------------------------------------------------- QLSCIndex

def test_index_add_point_accepts_skycoord():
	'''
	Test that add_point stores the ICRS position of a SkyCoord given in any frame.
	'''
	index = QLSCIndex(qlsc=QLSC(depth=30))
	coord = SkyCoord(ra=M31_RA * u.deg, dec=M31_DEC * u.deg, frame="icrs")

	index.add_point(coord.galactic, key="m31")

	stored = index.find_by_key("m31")
	assert stored["ra"][0] == pytest.approx(M31_RA)
	assert stored["dec"][0] == pytest.approx(M31_DEC)

def test_index_add_points_accepts_skycoord():
	'''
	Test that add_points takes an array SkyCoord as 'ra' or as 'points', and a scalar
	SkyCoord as a single point.
	'''
	ra = np.array([M31_RA, 200., 5.])
	dec = np.array([M31_DEC, -30., 60.])
	coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")

	for keyword in ("ra", "points"):
		index = QLSCIndex(qlsc=QLSC(depth=30))
		index.add_points(**{keyword: coord}, keys=["a", "b", "c"])
		assert len(index) == 3
		assert index.find_by_key("a")["ra"][0] == pytest.approx(M31_RA)

	index = QLSCIndex(qlsc=QLSC(depth=30))
	index.add_points(coord[0], keys=["only"])
	assert len(index) == 1

def test_index_add_points_skycoord_and_points_together_raise():
	'''
	Test that a SkyCoord and a 'points' array cannot both be given to add_points.
	'''
	coord = SkyCoord(ra=[1.] * u.deg, dec=[2.] * u.deg, frame="icrs")
	index = QLSCIndex(qlsc=QLSC(depth=30))

	with pytest.raises(ValueError, match="Only 'ra','dec' OR 'points'"):
		index.add_points(ra=coord, points=np.array([[1., 2.]]))

def test_index_remove_point_accepts_skycoord():
	'''
	Test that remove_point takes a SkyCoord for the coordinate to remove.
	'''
	index = index_with_points()
	coord = SkyCoord(ra=M31_RA * u.deg, dec=M31_DEC * u.deg, frame="icrs")

	assert index.remove_point(coord) == 1
	assert len(index) == 2

def test_index_radial_query_accepts_skycoord_center():
	'''
	Test that radial_query takes a SkyCoord center, in any frame, with the radius by
	keyword, and finds the same points as the degrees form.
	'''
	index = index_with_points()
	coord = SkyCoord(ra=M31_RA * u.deg, dec=M31_DEC * u.deg, frame="icrs")

	expected = index.radial_query(ra=M31_RA, dec=M31_DEC, radius=1.)

	assert np.allclose(index.radial_query(coord, radius=1.), expected)
	assert np.allclose(index.radial_query(coord.galactic, radius=1.), expected)
	assert np.allclose(index.radial_query(coord, radius=60. * u.arcmin), expected)

def test_index_radial_query_missing_radius_raises():
	'''
	Test that omitting the radius is reported, not left to fail obscurely: 'radius' is
	only formally optional so a SkyCoord can be the single positional coordinate.
	'''
	index = index_with_points()

	with pytest.raises(ValueError, match="'radius'"):
		index.radial_query(ra=M31_RA, dec=M31_DEC)

def test_index_radial_query_missing_dec_raises():
	'''
	Test that omitting the declination is reported for the same reason.
	'''
	index = index_with_points()

	with pytest.raises(ValueError, match="Both 'ra' and 'dec'"):
		index.radial_query(ra=M31_RA, radius=1.)

def test_index_radial_query_return_skycoord():
	'''
	Test that return_skycoord gives an array SkyCoord of the matches - never squeezed,
	so a single match is a length-1 SkyCoord, unlike the default squeezed array.
	'''
	index = index_with_points()

	matches = index.radial_query(ra=M31_RA, dec=M31_DEC, radius=1., return_skycoord=True)
	assert isinstance(matches, SkyCoord)
	assert matches.frame.name == "icrs"
	assert len(matches) == 2

	one = index.radial_query(ra=200., dec=-30., radius=0.5, return_skycoord=True)
	assert len(one) == 1
	assert one.ra.deg[0] == pytest.approx(200.)

	# the default return value is unchanged
	assert index.radial_query(ra=200., dec=-30., radius=0.5).shape == (2,)

def test_index_radial_query_return_skycoord_no_matches():
	'''
	Test that a query with no matches returns an empty SkyCoord rather than failing.
	'''
	index = index_with_points()
	matches = index.radial_query(ra=300., dec=0., radius=0.1, return_skycoord=True)

	assert isinstance(matches, SkyCoord)
	assert len(matches) == 0

def test_index_radial_query_skycoord_with_key_and_distance_gives_a_qtable():
	'''
	Test that combining return_skycoord with the other flags gives a QTable carrying a
	"coord" SkyCoord column plus the requested fields, since a SkyCoord cannot be a
	field of a recarray.
	'''
	index = index_with_points()

	table = index.radial_query(ra=M31_RA, dec=M31_DEC, radius=1.,
							   return_key=True, return_distance=True, return_skycoord=True)

	assert isinstance(table, QTable)
	assert table.colnames == ["coord", "key", "distance"]
	assert isinstance(table["coord"], SkyCoord)
	assert len(table) == 2
	assert set(table["key"]) == {"m31", "near"}
	assert table["distance"].unit == u.deg

	recarray = index.radial_query(ra=M31_RA, dec=M31_DEC, radius=1.,
								  return_key=True, return_distance=True)
	assert np.allclose(np.sort(table["distance"].to_value(u.deg)), np.sort(recarray["distance"]))

def test_index_radial_query_skycoord_with_key_only():
	'''
	Test the QTable with just the key column, including a point stored without one.
	'''
	index = QLSCIndex(qlsc=QLSC(depth=30))
	index.add_point(ra=10., dec=20.)   # no key
	index.add_point(ra=10.001, dec=20., key="labelled")

	table = index.radial_query(ra=10., dec=20., radius=1., return_key=True, return_skycoord=True)

	assert table.colnames == ["coord", "key"]
	assert sorted(table["key"], key=lambda value: (value is not None, value)) == [None, "labelled"]

def test_index_polygon_query_skycoord_in_and_out():
	'''
	Test that polygon_query takes a SkyCoord polygon and returns a SkyCoord on request.
	'''
	index = index_with_points()
	scheme = QLSC(depth=6)
	ipix = scheme.ang2ipix(M31_RA, M31_DEC)
	polygon = scheme.ipix2polygon(ipix)
	polygon_coord = scheme.ipix2polygon(ipix, return_skycoord=True)

	expected = index.polygon_query(polygon=polygon)
	assert len(expected) > 0

	assert np.allclose(index.polygon_query(polygon=polygon_coord), expected)

	matches = index.polygon_query(polygon=polygon_coord, return_skycoord=True)
	assert isinstance(matches, SkyCoord)
	assert len(matches) == len(expected)

	table = index.polygon_query(polygon=polygon_coord, return_key=True, return_skycoord=True)
	assert isinstance(table, QTable)
	assert table.colnames == ["coord", "key"]

def test_index_ellipse_query_skycoord_in_and_out():
	'''
	Test the SkyCoord center and return for both ellipse_query paths: the genuinely
	elliptical one and the axis_ratio=1 circle, which delegates to radial_query.
	'''
	index = index_with_points()
	coord = SkyCoord(ra=M31_RA * u.deg, dec=M31_DEC * u.deg, frame="icrs")

	for axis_ratio in (0.4, 1.0):
		expected = index.ellipse_query(ra=M31_RA, dec=M31_DEC, semi_major_axis=1.,
									   axis_ratio=axis_ratio)
		by_coord = index.ellipse_query(coord.galactic, semi_major_axis=1., axis_ratio=axis_ratio)
		assert np.allclose(by_coord, expected)

		matches = index.ellipse_query(coord, semi_major_axis=1., axis_ratio=axis_ratio,
									  return_skycoord=True)
		assert isinstance(matches, SkyCoord)
		assert len(matches) == len(expected)

		table = index.ellipse_query(coord, semi_major_axis=1., axis_ratio=axis_ratio,
									return_distance=True, return_skycoord=True)
		assert isinstance(table, QTable)
		assert table.colnames == ["coord", "distance"]

def test_index_find_by_key_return_skycoord():
	'''
	Test that find_by_key hands back a SkyCoord on request.
	'''
	index = index_with_points()

	coord = index.find_by_key("m31", return_skycoord=True)
	assert isinstance(coord, SkyCoord)
	assert len(coord) == 1
	assert coord.ra.deg[0] == pytest.approx(M31_RA)

	# the default is unchanged
	assert isinstance(index.find_by_key("m31"), np.recarray)

# ------------------------------------------------------------------ astropy not installed

def test_return_skycoord_without_astropy_raises_importerror(monkeypatch):
	'''
	Test that asking for a SkyCoord back without astropy installed raises a clear
	ImportError naming the package, rather than a NameError from the guarded import.

	astropy is present here (this module skips otherwise), so its absence is simulated
	by flipping the flag the package's guarded import sets - the same flag every
	astropy path is gated on.
	'''
	monkeypatch.setattr(qlsc.qlsc_base, "astropy_available", False)

	index = index_with_points()
	ipix = q30.ang2ipix(M31_RA, M31_DEC)

	for call in (lambda: q30.ipix2ang(ipix, return_skycoord=True),
				 lambda: q30.ipix2ang_center(ipix, return_skycoord=True),
				 lambda: q30.xy2ang(facenum=0, x=0.1, y=0.2, return_skycoord=True),
				 lambda: q30.ipix2polygon(ipix, return_skycoord=True),
				 lambda: index.find_by_key("m31", return_skycoord=True),
				 lambda: index.radial_query(ra=M31_RA, dec=M31_DEC, radius=1., return_skycoord=True)):
		with pytest.raises(ImportError, match="astropy"):
			call()

def test_coordinate_object_detection_is_gated_on_availability(monkeypatch):
	'''
	Test that the SkyCoord input path is gated on the same availability flag, so that a
	build without astropy never reaches an astropy name.

	This state cannot arise in practice - a SkyCoord cannot exist without astropy - but
	it is what the gate promises, and asserting it keeps the guard from decaying into a
	branch nothing checks.
	'''
	coord = SkyCoord(ra=M31_RA * u.deg, dec=M31_DEC * u.deg, frame="icrs")
	assert qlsc.qlsc_base._is_coordinate_object(coord) is True

	monkeypatch.setattr(qlsc.qlsc_base, "astropy_available", False)
	assert qlsc.qlsc_base._is_coordinate_object(coord) is False

	# and with the gate closed the object is refused loudly, not read as a number
	with pytest.raises((TypeError, ValueError)):
		q30.ang2ipix(coord)
