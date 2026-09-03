
import inspect

import numpy as np
import pytest

from qlsc import QLSC, GeoQLSC, GeoQLSCIndex

# A gate over the geographic classes: every public entry point that takes a latitude must
# REFUSE one outside [-90,90], not normalize it.
#
# This exists because the invariant is mechanically decidable and was broken. The
# geographic classes were written to reject an impossible latitude — GeoQLSCIndex did,
# face_number and ang2ipix_xy did — but GeoQLSC.ang2ipix delegated to the astronomical
# normalization, so lonlat2ipix(lon=10, lat=100) quietly answered for the point at
# (190, 80): a real place, on the other side of the Earth, with no error. point_in_polygon
# and point_in_ellipse did the same. Fixing those three would have left the rule
# unenforced, so the rule is enforced here instead.
#
# The list of methods to check is not written out. It is DERIVED from the signatures, so a
# new geographic method that takes a latitude fails this gate until it is added to the
# calls below and passes it. That is the whole point: the next method to be written is the
# one this catches.
#
# Blind spot, stated rather than assumed covered: a parameter that carries a latitude
# under some other name (a polygon's vertices, say) is not found by the signature scan.
# Polygon vertices are checked by the shared _polygon_ra_dec, which rejects any dec
# outside [-90,90] on both the astronomical and geographic paths, and that is pinned by
# test_geo_polygon_vertices_are_rejected below.

# the parameter names that carry a latitude
LATITUDE_PARAMETERS = {"dec", "lat", "dec_cen", "lat_cen"}

OFF_THE_SPHERE = (100., -100., 90.001)

def latitude_taking_methods(cls) -> set:
	'''
	Return the names of every public method of a class whose signature takes a latitude.
	'''
	found = set()
	for name, member in inspect.getmembers(cls, predicate=inspect.isfunction):
		if name.startswith("_"):
			continue
		parameters = set(inspect.signature(member).parameters)
		if parameters & LATITUDE_PARAMETERS:
			found.add(name)
	return found

def geo_scheme_calls(latitude:float) -> dict:
	'''
	One call per latitude-taking GeoQLSC method, each given a latitude off the sphere.

	The good latitude used alongside it is 80 degrees, so that nothing else about the call
	is out of range and the only reason to raise is the latitude under test.
	'''
	scheme = GeoQLSC(depth=10)
	polygon = scheme.ipix2polygon(scheme.lonlat2ipix(lon=10., lat=80.))

	return {
		"ang2ipix"           : lambda: scheme.ang2ipix(lon=10., lat=latitude),
		"ang2ipix_xy"        : lambda: scheme.ang2ipix_xy(lon=10., lat=latitude),
		"ang2ipix_at_depth"  : lambda: scheme.ang2ipix_at_depth(10., latitude, depth=5),
		"lonlat2ipix"        : lambda: scheme.lonlat2ipix(lon=10., lat=latitude),
		"face_number"        : lambda: scheme.face_number(lon=10., lat=latitude),
		"nearby_ipix_ranges" : lambda: scheme.nearby_ipix_ranges(lon=10., lat=latitude, radius=1.),
		"point_in_polygon"   : lambda: scheme.point_in_polygon(lon=10., lat=latitude, polygon=polygon),
		"point_in_ellipse"   : lambda: scheme.point_in_ellipse(lon=10., lat=latitude, lon_cen=10., lat_cen=80.,
															   semi_major_axis=1., axis_ratio=0.5),
	}

def geo_index_calls(latitude:float) -> dict:
	'''
	One call per latitude-taking GeoQLSCIndex method, each given a latitude off the sphere.
	'''
	index = GeoQLSCIndex(qlsc=GeoQLSC(depth=10))
	index.add_point(lon=10., lat=80., key="somewhere")

	return {
		"add_point"    : lambda: index.add_point(lon=10., lat=latitude),
		"add_points"   : lambda: index.add_points(lon=np.array([10.]), lat=np.array([latitude])),
		"remove_point" : lambda: index.remove_point(lon=10., lat=latitude),
		"radial_query" : lambda: index.radial_query(lon=10., lat=latitude, radius=1.),
		"ellipse_query": lambda: index.ellipse_query(lon=10., lat=latitude, semi_major_axis=1., axis_ratio=0.5),
	}

@pytest.mark.parametrize("cls,calls", [(GeoQLSC, geo_scheme_calls), (GeoQLSCIndex, geo_index_calls)])
def test_every_latitude_taking_method_is_covered(cls, calls):
	'''
	Test that the calls below cover every method the signature scan finds.

	A geographic method added later that takes a latitude fails here until it is given a
	call, which is what keeps the rest of this file from going stale.
	'''
	expected = latitude_taking_methods(cls)
	covered = set(calls(100.))

	assert covered == expected, (
		f"the latitude gate does not match {cls.__name__}'s methods.\n"
		f"  not covered by this file: {sorted(expected - covered)}\n"
		f"  covered but not found:    {sorted(covered - expected)}")

@pytest.mark.parametrize("latitude", OFF_THE_SPHERE)
@pytest.mark.parametrize("calls", [geo_scheme_calls, geo_index_calls])
def test_a_latitude_off_the_sphere_is_refused(calls, latitude):
	'''
	Test that every geographic entry point refuses a latitude outside [-90,90].
	'''
	failures = []
	for name, call in calls(latitude).items():
		try:
			call()
		except ValueError:
			continue
		except Exception as e:
			failures.append(f"{name} raised {type(e).__name__} rather than ValueError: {e}")
		else:
			failures.append(f"{name} accepted lat={latitude} and answered for some other point")

	assert not failures, "latitudes off the sphere were not refused:\n  " + "\n  ".join(failures)

@pytest.mark.parametrize("calls", [geo_scheme_calls, geo_index_calls])
def test_the_gate_passes_for_a_latitude_on_the_sphere(calls):
	'''
	Test that the calls above succeed with a good latitude; a gate whose calls raise for
	their own reasons would pass the test above while checking nothing.
	'''
	for name, call in calls(80.).items():
		if name == "remove_point":
			continue # nothing is stored at (10,80) in the scheme fixture except by add_point
		call()

def test_both_spellings_are_refused():
	'''
	Test that the strictness follows the meaning, not the parameter name: on a geographic
	object 'dec' is only another name for 'lat', so it is refused too.
	'''
	scheme = GeoQLSC(depth=10)
	index = GeoQLSCIndex(qlsc=GeoQLSC(depth=10))

	for call in (lambda: scheme.ang2ipix(ra=10., dec=100.),
				 lambda: scheme.lonlat2ipix(10., 100.),
				 lambda: index.add_point(ra=10., dec=100.)):
		with pytest.raises(ValueError, match=r"\[-90,90\]"):
			call()

def test_points_arrays_are_refused():
	'''
	Test the 'points' form, which does not pass through the shared coordinate resolver and
	so needs its own check.
	'''
	scheme = GeoQLSC(depth=10)
	index = GeoQLSCIndex(qlsc=GeoQLSC(depth=10))
	points = np.array([[10., 80.], [20., 100.]])

	with pytest.raises(ValueError, match=r"\[-90,90\]"):
		scheme.ang2ipix(points=points)
	with pytest.raises(ValueError, match=r"\[-90,90\]"):
		scheme.lonlat2ipix(points=points)
	with pytest.raises(ValueError, match=r"\[-90,90\]"):
		index.add_points(points=points)

	# and the same array with every latitude on the sphere is accepted
	good = np.array([[10., 80.], [20., -80.]])
	assert len(scheme.ang2ipix(points=good)) == 2
	index.add_points(points=good)
	assert len(index) == 2

def test_geo_polygon_vertices_are_rejected():
	'''
	Test the blind spot named at the top of this file: a latitude carried as a polygon
	vertex rather than as a named parameter is refused as well.
	'''
	scheme = GeoQLSC(depth=10)
	index = GeoQLSCIndex(qlsc=GeoQLSC(depth=10))
	polygon = np.array([[9., 99.], [11., 99.], [11., 101.]])

	with pytest.raises(ValueError, match=r"\[-90,90\]"):
		scheme.point_in_polygon(lon=10., lat=80., polygon=polygon)
	with pytest.raises(ValueError, match=r"\[-90,90\]"):
		index.polygon_query(polygon=polygon)

def test_the_astronomical_classes_still_normalize():
	'''
	Test that this strictness is confined to the geographic classes: an unnormalized
	declination remains a legitimate celestial coordinate, and QLSC keeps reading it as
	the point past the pole.
	'''
	astronomical = QLSC(depth=10)

	assert astronomical.ang2ipix(10., 100.) == astronomical.ang2ipix(190., 80.)
	assert astronomical.ang2ipix_at_depth(10., 100., depth=5) == \
		   astronomical.ang2ipix_at_depth(190., 80., depth=5)

def test_ang2ipix_at_depth_agrees_with_ang2ipix():
	'''
	Test that ang2ipix_at_depth reads a coordinate the way ang2ipix does.

	It used to call the C function directly, which truncated a declination past the pole
	to the pole itself instead of normalizing it — a different pixel from the one
	ang2ipix names for the same coordinate — and let a non-finite value through to the
	undefined float-to-int conversion of upstream q3c issue #52.
	'''
	astronomical = QLSC(depth=30)

	for ra, dec in ((10., 20.), (350., -85.), (10., 100.), (-350., 20.)):
		expected = astronomical.ang2ipix(ra, dec)
		at_depth = astronomical.ang2ipix_at_depth(ra, dec, depth=1)
		# depth 1 keeps all but the lowest two bits, plus q3c's centering term
		assert at_depth >> 2 == expected >> 2

	for dec in (float("nan"), float("inf")):
		with pytest.raises(ValueError):
			astronomical.ang2ipix_at_depth(10., dec, depth=5)
