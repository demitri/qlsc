
import numpy as np
import pytest

from qlsc import QLSC, QLSCIndex, GeoQLSC, GeoQLSCIndex, wrap_longitude
from qlsc.generate import sunflower_points_on_sphere

# Tests for the geographic (longitude/latitude) veneer: GeoQLSC, GeoQLSCIndex and the
# wrap_longitude helper.
#
# The mapping under test, as documented on GeoQLSC:
#   latitude  == declination, unchanged in both directions;
#   longitude == right ascension on the way in (any range is accepted, as it already was);
#   longitude == right ascension wrapped into [-180,180) on the way out.
# So every assertion below is either "the geographic call equals the astronomical one" or
# "the value that came back is in [-180,180)".

DEPTH = 20

geo = GeoQLSC(depth=DEPTH)
astro = QLSC(depth=DEPTH)

# a spread of positions that exercises the seam, the poles and the ordinary sky
POSITIONS = [(0., 0.), (-74.0060, 40.7128), (2.3522, 48.8566), (179.9, -12.5),
			 (-179.9, -12.5), (180., 0.), (-180., 0.), (45., 90.), (100., -90.),
			 (359.9, 5.), (286., 40.7128)]

# ------------------------------------------------------------------- wrap_longitude

def test_wrap_longitude_scalars():
	'''
	Test the longitude convention on the values that define it, including the half-open
	upper bound: +180 and -180 are the same meridian and both come back as -180.
	'''
	assert wrap_longitude(0.) == 0.
	assert wrap_longitude(74.) == 74.
	assert wrap_longitude(-74.) == -74.
	assert wrap_longitude(286.) == pytest.approx(-74.)
	assert wrap_longitude(180.) == -180.
	assert wrap_longitude(-180.) == -180.
	assert wrap_longitude(360.) == 0.
	assert wrap_longitude(-360.) == 0.
	assert wrap_longitude(540.) == -180.
	assert wrap_longitude(179.999) == pytest.approx(179.999)

def test_wrap_longitude_is_idempotent_and_in_range():
	'''
	Test that wrapping an already wrapped value changes nothing, and that the result is
	always in [-180,180).
	'''
	values = np.linspace(-1000., 1000., 501)
	wrapped = wrap_longitude(values)

	assert isinstance(wrapped, np.ndarray)
	assert ((-180. <= wrapped) & (wrapped < 180.)).all()
	assert np.allclose(wrap_longitude(wrapped), wrapped)

def test_wrap_longitude_preserves_the_position():
	'''
	Test that wrapping only changes the representation: the pixel a longitude falls in is
	the same before and after.
	'''
	for longitude, latitude in POSITIONS:
		assert astro.ang2ipix(wrap_longitude(longitude), latitude) == astro.ang2ipix(longitude, latitude)

# ------------------------------------------------------------------------- GeoQLSC

def test_lonlat2ipix_agrees_with_ang2ipix():
	'''
	Test that the geographic spelling gives the same pixel as the astronomical one, for
	longitudes on both sides of zero.
	'''
	for longitude, latitude in POSITIONS:
		expected = astro.ang2ipix(longitude, latitude)
		assert geo.lonlat2ipix(lon=longitude, lat=latitude) == expected
		assert geo.lonlat2ipix(longitude, latitude) == expected
		assert geo.ang2ipix(lon=longitude, lat=latitude) == expected
		assert geo.ang2ipix(longitude, latitude) == expected
		# a longitude 360 degrees away is the same place and takes the same path through
		# the C layer's unwrapping, so the two spellings agree there too (compared against
		# the astronomical answer for the same value, not against the wrapped one: at an
		# exact cube corner such as ra=0/360 on the equator, Q3C picks between the pixels
		# meeting there by the sign of a floating point zero, which the geographic veneer
		# neither causes nor cures)
		assert geo.lonlat2ipix(lon=longitude + 360., lat=latitude) == astro.ang2ipix(longitude + 360., latitude)

def test_lonlat2ipix_arrays_and_points():
	'''
	Test the array and 'points' forms of the geographic spelling.
	'''
	lon = np.array([point[0] for point in POSITIONS])
	lat = np.array([point[1] for point in POSITIONS])
	expected = astro.ang2ipix(ra=lon, dec=lat)

	assert np.array_equal(geo.lonlat2ipix(lon=lon, lat=lat), expected)
	assert np.array_equal(geo.lonlat2ipix(points=np.column_stack((lon, lat))), expected)

def test_ipix2lonlat_is_in_range_and_matches_ipix2ang():
	'''
	Test that every longitude reported by a GeoQLSC is in [-180,180) and is the wrapped
	form of what the astronomical class reports for the same pixel.
	'''
	ipix_values = [astro.ang2ipix(longitude, latitude) for longitude, latitude in POSITIONS]

	for ipix in ipix_values:
		ra, dec = astro.ipix2ang(ipix)
		lon, lat = geo.ipix2lonlat(ipix)

		assert -180. <= lon < 180.
		assert lon == pytest.approx(wrap_longitude(ra))
		assert lat == pytest.approx(dec)
		assert geo.ipix2lonlat_center(ipix) == (lon, lat)
		# the inherited astronomical name follows the same convention on this object
		assert geo.ipix2ang(ipix) == (lon, lat)

def test_ipix2lonlat_arrays():
	'''
	Test the array form of ipix2lonlat.
	'''
	ipix = np.array([astro.ang2ipix(longitude, latitude) for longitude, latitude in POSITIONS])
	lon, lat = geo.ipix2lonlat(ipix)
	ra, dec = astro.ipix2ang(ipix)

	assert ((-180. <= lon) & (lon < 180.)).all()
	assert np.allclose(lon, wrap_longitude(ra))
	assert np.allclose(lat, dec)

def test_lonlat_round_trip():
	'''
	Test that a pixel survives the round trip through geographic coordinates.
	'''
	points = sunflower_points_on_sphere(n=500)
	for ra, dec in points:
		ipix = geo.lonlat2ipix(lon=wrap_longitude(ra), lat=dec)
		lon, lat = geo.ipix2lonlat(ipix)
		assert geo.lonlat2ipix(lon=lon, lat=lat) == ipix

def test_xy2ang_and_ipix2polygon_follow_the_convention():
	'''
	Test that every method computing a longitude follows the object's convention, not
	just the ones with a geographic name: xy2ang and the polygon corners it builds.
	'''
	for facenum in range(6):
		for x, y in ((0., 0.), (-0.9, 0.9), (0.9, -0.9)):
			lon, lat = geo.xy2ang(facenum=facenum, x=x, y=y)
			ra, dec = astro.xy2ang(facenum=facenum, x=x, y=y)
			assert -180. <= lon < 180.
			assert lon == pytest.approx(wrap_longitude(ra))
			assert lat == pytest.approx(dec)

	# a pixel on the seam, where the two conventions actually differ
	ipix = astro.ang2ipix(359.9, 10.)
	polygon = geo.ipix2polygon(ipix)
	assert ((-180. <= polygon[:,0]) & (polygon[:,0] < 180.)).all()
	assert np.allclose(polygon[:,0], wrap_longitude(astro.ipix2polygon(ipix)[:,0]))

def test_polygon_from_geo_is_accepted_back():
	'''
	Test that a polygon reported in the geographic convention can be handed straight back
	to point_in_polygon, seam included - longitudes need no conversion on the way in.
	'''
	for longitude, latitude in ((10., 20.), (179.95, -5.), (-179.95, -5.)):
		ipix = geo.lonlat2ipix(lon=longitude, lat=latitude)
		polygon = geo.ipix2polygon(ipix)
		center_lon, center_lat = geo.ipix2lonlat(ipix)
		assert geo.point_in_polygon(lon=center_lon, lat=center_lat, polygon=polygon) is True

def test_other_geo_methods_agree_with_the_astronomical_spelling():
	'''
	Test face_number, ang2ipix_xy and point_in_ellipse under both spellings.
	'''
	for longitude, latitude in POSITIONS:
		assert geo.face_number(lon=longitude, lat=latitude) == astro.face_number(longitude, latitude)
		if abs(latitude) < 90:
			assert geo.ang2ipix_xy(lon=longitude, lat=latitude) == astro.ang2ipix_xy(ra=longitude, dec=latitude)

	assert geo.point_in_ellipse(lon=-74.5, lat=40.7, lon_cen=-74.0, lat_cen=40.7,
								semi_major_axis=1., axis_ratio=0.5, position_angle=90.) is \
		   astro.point_in_ellipse(ra=-74.5, dec=40.7, ra_cen=-74.0, dec_cen=40.7,
								  semi_major_axis=1., axis_ratio=0.5, position_angle=90.)

def test_both_spellings_at_once_raise():
	'''
	Test that giving a coordinate in both spellings is refused rather than one of them
	being silently dropped.
	'''
	with pytest.raises(ValueError, match="two spellings"):
		geo.ang2ipix(ra=10., dec=20., lon=30.)
	with pytest.raises(ValueError, match="two spellings"):
		geo.face_number(ra=10., dec=20., lat=20.)
	with pytest.raises(ValueError, match="two spellings"):
		geo.point_in_ellipse(lon=1., lat=2., ra_cen=3., dec_cen=4., lon_cen=5.,
							 semi_major_axis=1., axis_ratio=0.5)

def test_points_and_lonlat_together_raise():
	'''
	Test that 'points' and the geographic spelling cannot both be given.
	'''
	with pytest.raises(ValueError, match="'points'"):
		geo.ang2ipix(lon=10., lat=20., points=np.array([[1., 2.]]))

# -------------------------------------------------------------------- GeoQLSCIndex

def test_index_requires_a_geo_scheme():
	'''
	Test that a GeoQLSCIndex refuses a plain QLSC: the longitude convention is carried by
	the scheme object, so accepting one would produce an index that reports [0,360).
	'''
	with pytest.raises(ValueError, match="must be a GeoQLSC"):
		GeoQLSCIndex(qlsc=QLSC(depth=DEPTH))

	index = GeoQLSCIndex()
	assert isinstance(index.qlsc, GeoQLSC)
	assert index.qlsc.depth == 30

def test_index_stores_and_returns_wrapped_longitudes():
	'''
	Test that a point is stored with its longitude wrapped, whichever equivalent value it
	was given with, and that a query returns it in that form.
	'''
	index = GeoQLSCIndex(qlsc=GeoQLSC(depth=DEPTH))
	index.add_point(lon=286., lat=40.7128, key="new york")

	stored = index.find_by_key("new york")
	assert stored["ra"][0] == pytest.approx(-74.)

	matches = index.radial_query(lon=-74., lat=40.7128, radius=1.)
	assert matches.shape == (2,)
	assert matches[0] == pytest.approx(-74.)

def test_index_recognizes_equivalent_longitudes_as_one_point():
	'''
	Test that the same place given with either longitude is one point, not two: the
	wrapping happens before the index's uniqueness constraint sees the coordinate.
	'''
	index = GeoQLSCIndex(qlsc=GeoQLSC(depth=DEPTH))
	index.add_point(lon=-74., lat=40.)
	index.add_point(lon=286., lat=40.)

	assert len(index) == 1

def test_index_add_points_forms():
	'''
	Test the array and 'points' forms of add_points, with keys.
	'''
	lon = np.array([-74., 286., 2.3522])
	lat = np.array([40.7128, 40.7128, 48.8566])

	by_arrays = GeoQLSCIndex(qlsc=GeoQLSC(depth=DEPTH))
	by_arrays.add_points(lon=lon, lat=lat, keys=["a", "b", "c"])
	# the first two are the same place, so the index holds two points
	assert len(by_arrays) == 2

	by_points = GeoQLSCIndex(qlsc=GeoQLSC(depth=DEPTH))
	by_points.add_points(points=np.column_stack((lon, lat)), keys=["a", "b", "c"])
	assert len(by_points) == 2

	assert by_points.find_by_key("a")["ra"][0] == pytest.approx(-74.)

def test_index_add_points_rejects_both_forms():
	'''
	Test that 'points' and the coordinate arrays cannot both be given.
	'''
	index = GeoQLSCIndex(qlsc=GeoQLSC(depth=DEPTH))
	with pytest.raises(ValueError, match="OR 'points'"):
		index.add_points(lon=np.array([1.]), lat=np.array([2.]), points=np.array([[1., 2.]]))

def test_index_rejects_a_latitude_out_of_range():
	'''
	Test that a latitude outside [-90,90] is an error here, rather than being read as a
	point on the far side of a pole with its longitude moved by 180 degrees - which is what
	the astronomical classes do with an unnormalized declination.
	'''
	index = GeoQLSCIndex(qlsc=GeoQLSC(depth=DEPTH))

	with pytest.raises(ValueError, match=r"\[-90,90\]"):
		index.add_point(lon=10., lat=100.)
	with pytest.raises(ValueError, match=r"\[-90,90\]"):
		index.add_points(lon=np.array([10., 20.]), lat=np.array([10., 91.]))
	with pytest.raises(ValueError, match=r"\[-90,90\]"):
		index.remove_point(lon=10., lat=-91.)

	# the astronomical index does relocate it, which is the behaviour being diverged from
	astronomical = QLSCIndex(qlsc=QLSC(depth=DEPTH))
	astronomical.add_point(ra=10., dec=100.)
	assert len(astronomical) == 1

def test_index_remove_point_accepts_either_longitude():
	'''
	Test that a point can be removed with any equivalent longitude, since the coordinate
	given is put into the stored convention first.
	'''
	index = GeoQLSCIndex(qlsc=GeoQLSC(depth=DEPTH))
	index.add_point(lon=-74., lat=40.)
	assert index.remove_point(lon=286., lat=40.) == 1
	assert len(index) == 0

def test_index_query_across_the_antimeridian():
	'''
	Test a cone search centered on the antimeridian, where the two conventions differ and
	the longitudes of the matches have opposite signs.
	'''
	index = GeoQLSCIndex(qlsc=GeoQLSC(depth=DEPTH))
	index.add_point(lon=179.9, lat=-12.5, key="west")
	index.add_point(lon=-179.9, lat=-12.5, key="east")
	index.add_point(lon=0., lat=-12.5, key="far")

	matches = index.radial_query(lon=180., lat=-12.5, radius=0.5, return_key=True)

	assert sorted(matches["key"]) == ["east", "west"]
	assert sorted(np.round(matches["ra"], 4)) == [-179.9, 179.9]

def test_index_agrees_with_the_astronomical_index():
	'''
	Test that the two indices answer the same question identically under the documented
	mapping: the same points, the same cone, the same matches once the longitudes are put
	in the same convention.
	'''
	points = sunflower_points_on_sphere(n=2000)

	astronomical = QLSCIndex(qlsc=QLSC(depth=DEPTH))
	astronomical.add_points(points=points)

	geographic = GeoQLSCIndex(qlsc=GeoQLSC(depth=DEPTH))
	geographic.add_points(points=np.column_stack((wrap_longitude(points[:,0]), points[:,1])))

	assert len(geographic) == len(astronomical)

	for longitude, latitude, radius in ((-74., 40., 8.), (0., 0., 5.), (180., -30., 10.), (30., 89., 15.)):
		expected = astronomical.radial_query(ra=longitude, dec=latitude, radius=radius)
		found = geographic.radial_query(lon=longitude, lat=latitude, radius=radius)

		expected = np.asarray(expected, dtype=np.double).reshape(-1, 2)
		found = np.asarray(found, dtype=np.double).reshape(-1, 2)
		assert len(found) == len(expected)
		assert len(found) > 0

		expected_sorted = expected[np.lexsort((expected[:,1], wrap_longitude(expected[:,0])))]
		found_sorted = found[np.lexsort((found[:,1], found[:,0]))]

		assert np.allclose(wrap_longitude(expected_sorted[:,0]), found_sorted[:,0])
		assert np.allclose(expected_sorted[:,1], found_sorted[:,1])

def test_index_ellipse_and_polygon_queries():
	'''
	Test that the elliptical and polygonal searches work in the geographic spelling and
	report longitudes in the stored convention.
	'''
	index = GeoQLSCIndex(qlsc=GeoQLSC(depth=DEPTH))
	index.add_points(lon=np.array([-74.0, -74.2, 100.]), lat=np.array([40.0, 40.1, -20.]),
					 keys=["a", "b", "far"])

	ellipse = index.ellipse_query(lon=-74.1, lat=40.05, semi_major_axis=1., axis_ratio=0.5,
								  return_key=True)
	assert sorted(ellipse["key"]) == ["a", "b"]
	assert ((-180. <= ellipse["ra"]) & (ellipse["ra"] < 180.)).all()

	# polygon_query is inherited unchanged: its vertices are [lon,lat] pairs under the
	# same mapping, so a box written in negative longitudes needs no conversion
	polygon = np.array([[-74.4, 39.8], [-73.8, 39.8], [-73.8, 40.3], [-74.4, 40.3]])
	matches = index.polygon_query(polygon=polygon, return_key=True)
	assert sorted(matches["key"]) == ["a", "b"]
	assert ((-180. <= matches["ra"]) & (matches["ra"] < 180.)).all()

def test_index_file_keeps_the_geographic_scheme(tmp_path):
	'''
	Test that reopening an index file does not quietly replace the GeoQLSC with a plain
	QLSC when the file's depth differs from the scheme handed to the constructor - which
	would silently change the convention the index reports in.
	'''
	filepath = tmp_path / "geo.qlscidx"

	written = GeoQLSCIndex(qlsc=GeoQLSC(depth=10), filepath=filepath)
	written.add_point(lon=286., lat=40., key="here")
	del written

	reopened = GeoQLSCIndex(qlsc=GeoQLSC(depth=DEPTH), filepath=filepath)

	assert reopened.qlsc.depth == 10
	assert isinstance(reopened.qlsc, GeoQLSC)
	assert reopened.find_by_key("here")["ra"][0] == pytest.approx(-74.)
	assert reopened.qlsc.ipix2lonlat(reopened.qlsc.lonlat2ipix(lon=-74., lat=40.))[0] < 0.

def test_poles_and_seam_boundaries():
	'''
	Test the boundary cases the geographic convention is most likely to get wrong: a
	longitude of exactly +/-180, and the poles, where the longitude is arbitrary.
	'''
	index = GeoQLSCIndex(qlsc=GeoQLSC(depth=DEPTH))
	index.add_point(lon=180., lat=0., key="antimeridian")
	index.add_point(lon=0., lat=90., key="north pole")
	index.add_point(lon=0., lat=-90., key="south pole")

	# +180 is stored as -180, the same meridian
	assert index.find_by_key("antimeridian")["ra"][0] == pytest.approx(-180.)
	# and is found from either side
	assert len(index.radial_query(lon=-180., lat=0., radius=0.01)) > 0
	assert len(index.radial_query(lon=180., lat=0., radius=0.01)) > 0

	# A pole is the point where four pixels meet, and which of them Q3C answers with
	# depends on the sign of a floating point zero derived from the longitude - the same
	# in both spellings, since the geographic layer changes nothing on the way in.
	for longitude in (0., -179., 137., 359.):
		assert geo.lonlat2ipix(lon=longitude, lat=90.) == astro.ang2ipix(longitude, 90.)

	north = index.radial_query(lon=137., lat=90., radius=0.01, return_key=True)
	assert list(north["key"]) == ["north pole"]

	lon, lat = geo.ipix2lonlat(geo.lonlat2ipix(lon=0., lat=-90.))
	assert -180. <= lon < 180.
	assert lat < -89.
