#!/usr/bin/env python

'''
Tests for the region and image queries qlsc answers through cornish/Starlink AST:
:py:func:`QLSCIndex.region_query` and :py:func:`QLSCIndex.image_query`.

The point of every agreement test here is the same one: the fast paths hand the whole
query to one of qlsc's own region queries, on the claim that qlsc and AST describe the
same shape, so the fast-path answer must equal the answer AST gives when it is the sole
arbiter of membership (``fast_path=False``, where AST's ``pointInRegion`` decides every
candidate). The catalogues are the deterministic sunflower point set, so a disagreement
is reproducible rather than a flake.

Regions crossing ra=0 and containing the poles get their own cases: that is where the
coordinate seams are, and where a coarse step that quietly loses points would show up.

The tests that need cornish skip with a named reason when it is not installed; the ones
that check the missing-package path run always.
'''

import numpy as np
import pytest

import qlsc
from qlsc import QLSC, QLSCIndex, q3c
from qlsc import cornish_interop
from qlsc.generate import sunflower_points_on_sphere

try:
	import cornish
	from cornish import (ASTBox, ASTCircle, ASTPolygon, ASTICRSFrame, ASTCompoundRegion,
						 ASTRegion, RegionOperation)
	cornish_available = True
except ImportError:
	cornish_available = False

CORNISH_SKIP_REASON = "the optional 'cornish' package is not installed"
requires_cornish = pytest.mark.skipif(not cornish_available, reason=CORNISH_SKIP_REASON)

try:
	from astropy.io import fits
	astropy_available = True
except ImportError:
	astropy_available = False

requires_astropy = pytest.mark.skipif(not astropy_available, reason="astropy is not installed")


@pytest.fixture(scope="module")
def index():
	'''An index of evenly spread points over the whole sphere, at a modest depth.'''
	index = QLSCIndex(qlsc=QLSC(depth=10))
	index.add_points(points=sunflower_points_on_sphere(n=20000))
	yield index
	index.close()

@pytest.fixture(scope="module")
def keyed_index():
	'''A smaller index whose points carry keys, for the return_key paths.'''
	index = QLSCIndex(qlsc=QLSC(depth=10))
	points = sunflower_points_on_sphere(n=2000)
	index.add_points(points=points, keys=[f"star-{i}" for i in range(len(points))])
	yield index
	index.close()

@pytest.fixture(scope="module")
def frame():
	return ASTICRSFrame() if cornish_available else None


def sorted_pairs(matches:np.ndarray) -> np.ndarray:
	'''Order a (n,2) match array so two answers can be compared regardless of row order.'''
	matches = np.asarray(matches, dtype=np.double).reshape(-1, 2)
	return matches[np.lexsort((matches[:,1], matches[:,0]))]

def assert_paths_agree(index, region, minimum_matches:int=1):
	'''
	Assert the fast path and the generic AST-filtered path return the same matches.

	:param index: the index to query
	:param region: the region to query it with
	:param minimum_matches: the number of matches below which the case is not testing anything
	'''
	fast = index.region_query(region)
	generic = index.region_query(region, fast_path=False)

	assert len(fast) >= minimum_matches, \
		f"the region matched {len(fast)} points, too few for this comparison to mean anything"
	np.testing.assert_array_equal(sorted_pairs(fast), sorted_pairs(generic))


# ---------------------------------------------------------------------------
# The missing-package path: these run whether or not cornish is installed.
# ---------------------------------------------------------------------------

def test_region_query_method_without_cornish_raises_import_error(monkeypatch):
	'''
	With cornish absent, the index method raises an ImportError naming the package and
	how to install it - never a silent fallback to some other answer.
	'''
	original = ImportError("No module named 'cornish'")
	monkeypatch.setattr(cornish_interop, "cornish_available", False)
	monkeypatch.setattr(cornish_interop, "cornish_import_error", original)

	index = QLSCIndex(qlsc=QLSC(depth=5))
	with pytest.raises(ImportError) as e:
		index.region_query(object())

	assert "cornish" in str(e.value)
	assert "pip install cornish" in str(e.value)
	assert e.value.__cause__ is original

def test_image_query_method_without_cornish_raises_import_error(monkeypatch):
	'''The same for image_query, which is where a FITS-header user meets the bridge first.'''
	original = ImportError("No module named 'cornish'")
	monkeypatch.setattr(cornish_interop, "cornish_available", False)
	monkeypatch.setattr(cornish_interop, "cornish_import_error", original)

	index = QLSCIndex(qlsc=QLSC(depth=5))
	with pytest.raises(ImportError) as e:
		index.image_query({"NAXIS": 2})

	assert "cornish" in str(e.value)
	assert e.value.__cause__ is original

def test_region_query_import_is_deferred():
	'''
	qlsc_base must not import cornish_interop (and so cornish) at module import time:
	cornish is optional, and cornish_interop imports qlsc_base in turn.
	'''
	import inspect
	source = inspect.getsource(qlsc.qlsc_base)
	# the delegation imports live inside the two methods, not at the top of the file
	assert "\nfrom . import cornish_interop" not in source
	assert source.count("\t\tfrom . import cornish_interop") == 2

def test_polygon_vertex_limit_is_read_from_the_extension():
	'''
	The vertex limit the fast path checks is q3c's own, published by the extension module
	rather than copied into Python where the two could drift apart.
	'''
	assert cornish_interop.MAXIMUM_POLYGON_VERTICES == q3c.Q3C_MAX_N_POLY_VERTEX
	assert isinstance(q3c.Q3C_MAX_N_POLY_VERTEX, int)


# ---------------------------------------------------------------------------
# The assumption the polygon fast path rests on.
# ---------------------------------------------------------------------------

def test_polygon_query_is_winding_agnostic():
	'''
	The polygon fast path hands AST's vertices to polygon_query without reversing them,
	which is only correct because Q3C's polygon test answers about the same side of the
	sphere whichever way round the vertices are given. AST does not: it fills the region
	to the left of the vertex path, which is why ipix2region reverses on the way out.
	'''
	index = QLSCIndex(qlsc=QLSC(depth=8))
	index.add_points(points=sunflower_points_on_sphere(n=20000))
	try:
		vertices = np.array([[44.,13.],[46.,13.],[46.,11.],[44.,11.]])
		forwards = index.polygon_query(polygon=vertices)
		backwards = index.polygon_query(polygon=vertices[::-1])
		assert len(forwards) > 0
		np.testing.assert_array_equal(sorted_pairs(forwards), sorted_pairs(backwards))
	finally:
		index.close()


# ---------------------------------------------------------------------------
# Exact (AST-filtered) versus fast-path agreement, region type by region type.
# ---------------------------------------------------------------------------

@requires_cornish
@pytest.mark.parametrize("centre,radius", [
	((45., 12.), 2.0),      # ordinary position, away from every seam
	((0., 0.), 3.0),        # centred on ra=0 on the equator
	((359.5, -0.2), 2.0),   # straddling ra=0
	((123., 89.5), 3.0),    # containing the north pole
	((260., -89.9), 5.0),   # containing the south pole
	((45., 12.), 1.5),      # a small cone, a few pixels across at this depth
	((10., 20.), 60.0),     # large: a sixth of the sky
])
def test_circle_fast_path_agrees_with_ast(index, frame, centre, radius):
	'''An AST circle is a geodesic disc, which is exactly what radial_query searches.'''
	assert_paths_agree(index, ASTCircle(frame=frame, center=list(centre), radius=radius))

@requires_cornish
@pytest.mark.parametrize("vertices", [
	[[44.,13.],[46.,13.],[46.,11.],[44.,11.]],          # ordinary square
	[[359.,1.],[1.,1.],[1.,-1.],[359.,-1.]],            # crossing ra=0
	[[270.,88.],[180.,88.],[90.,88.],[0.,88.]],         # a cap around the north pole
	[[90.,-88.],[180.,-88.],[270.,-88.],[0.,-88.]],     # a cap around the south pole
	[[40.,15.],[42.5,12.],[45.,15.],[45.,10.],[40.,10.]], # concave
])
def test_polygon_fast_path_agrees_with_ast(index, frame, vertices):
	'''
	An AST polygon joins its vertices with great circles, and so does Q3C's, so
	polygon_query answers about the same shape. The concave case matters because Q3C's
	coarse coverage step is generous for a concave polygon; the exact filter removes
	the extras, and the two paths must still land in the same place.
	'''
	assert_paths_agree(index, ASTPolygon(frame=frame, points=vertices))

@requires_cornish
@pytest.mark.parametrize("corners", [
	(([44.,11.],[46.,13.])),      # ordinary
	(([359.,-1.],[1.,1.])),       # crossing ra=0
	(([0.,86.],[360.,90.])),      # a band up to the north pole
	(([0.,-90.],[360.,-84.])),    # a band down to the south pole
])
def test_box_takes_the_generic_path_and_is_still_right(index, frame, corners):
	'''
	An AST box on a sky frame has edges parallel to the frame axes - parallels of
	declination, not great circles - so there is no polygon fast path for it and both
	calls run the generic path. The case is kept anyway: it pins that a box is answered
	correctly, and it is the regression test if anyone ever adds a box fast path by
	assuming its corners describe a spherical polygon.
	'''
	region = ASTBox.fromCorners(frame=frame, corners=corners)
	assert_paths_agree(index, region)

@requires_cornish
def test_ast_box_edges_are_parallels_of_declination(frame):
	'''
	The measurement behind the decision above: the midpoint of the box's lower edge is
	inside at exactly the box's declination and outside just below it, which is the
	parallel of declination and not the great circle joining the two corners (that arc
	would bow polewards, i.e. above the edge, at its midpoint).
	'''
	box = ASTBox.fromCorners(frame=frame, corners=([10.,20.],[20.,30.]))
	assert box.pointInRegion([15., 20.0])
	assert not box.pointInRegion([15., 19.98])
	assert box.pointInRegion([15., 20.02])

@requires_cornish
def test_compound_region_agrees(index, frame):
	'''A union of two circles has no fast path; both calls take the generic path.'''
	a = ASTCircle(frame=frame, center=[44., 12.], radius=2.0)
	b = ASTCircle(frame=frame, center=[47., 12.], radius=2.0)
	union = ASTCompoundRegion(regions=[a, b], operation=RegionOperation.OR)
	fast = index.region_query(union)
	generic = index.region_query(union, fast_path=False)
	np.testing.assert_array_equal(sorted_pairs(fast), sorted_pairs(generic))

@requires_cornish
def test_region_query_method_matches_the_module_function(index, frame):
	'''The method is the module function; a caller can use either.'''
	region = ASTCircle(frame=frame, center=[45., 12.], radius=2.0)
	np.testing.assert_array_equal(index.region_query(region),
								  cornish_interop.region_query(index, region))


# ---------------------------------------------------------------------------
# The fast paths' guards.
# ---------------------------------------------------------------------------

@requires_cornish
def test_a_polygon_with_too_many_vertices_takes_the_generic_path(index, frame):
	'''
	Q3C carries polygon vertices in fixed-size buffers; a polygon with more than it
	accepts is answered by the generic path rather than refused.
	'''
	n = q3c.Q3C_MAX_N_POLY_VERTEX + 20
	angle = np.linspace(0., 2*np.pi, n, endpoint=False)
	# a small circle of vertices, wound the way AST wants for a bounded interior
	vertices = np.column_stack((45. + 1.0*np.cos(angle), 12. + 1.0*np.sin(angle)))[::-1]
	region = ASTPolygon(frame=frame, points=vertices)

	assert len(region.points) > q3c.Q3C_MAX_N_POLY_VERTEX
	# polygon_query itself refuses it, so a fast path would have raised
	with pytest.raises(ValueError):
		index.polygon_query(polygon=region.points)

	assert_paths_agree(index, region)

@requires_cornish
def test_a_polygon_larger_than_the_fast_path_limit_still_answers(index, frame):
	'''
	A polygon whose bounding circle is beyond POLYGON_FAST_PATH_MAXIMUM_RADIUS takes the
	generic path; it is not refused, and the answer is the same one AST gives.
	'''
	region = ASTPolygon(frame=frame, points=[[20.,20.],[40.,20.],[40.,0.],[20.,0.]])
	assert region.boundingCircle().radius.to("deg").value > cornish_interop.POLYGON_FAST_PATH_MAXIMUM_RADIUS
	assert_paths_agree(index, region)

@requires_cornish
def test_the_polygon_fast_path_limit_is_inside_q3c_projection_limit(index, frame):
	'''
	The fast path's size limit is only useful if Q3C really can project everything under
	it, wherever it sits. The worst position measured for the projection is a face edge;
	a polygon of the limiting size there must go through polygon_query without raising.
	'''
	radius = cornish_interop.POLYGON_FAST_PATH_MAXIMUM_RADIUS
	for centre_ra, centre_dec in [(45., 0.), (0., 0.), (45., 44.), (0., 60.), (0., 75.)]:
		half = radius / np.sqrt(2.)
		scale = 1. / np.cos(np.radians(abs(centre_dec) + half))
		vertices = np.array([[centre_ra - half*scale, centre_dec + half],
							 [centre_ra + half*scale, centre_dec + half],
							 [centre_ra + half*scale, centre_dec - half],
							 [centre_ra - half*scale, centre_dec - half]])
		# raises ValueError if q3c cannot project it onto a single cube face
		index.polygon_query(polygon=vertices)

@requires_cornish
def test_a_range_overflow_falls_back_to_the_generic_path(index, frame, monkeypatch, caplog):
	'''
	polygon_query lets q3c's range-capacity overflow through - a polygon has no single
	size to grow and retry with. The region query does have another route to the same
	answer, so it takes it, says so in the log, and returns the matches rather than the
	error.
	'''
	region = ASTPolygon(frame=frame, points=[[44.,13.],[46.,13.],[46.,11.],[44.,11.]])
	expected = index.region_query(region, fast_path=False)

	def overflow(*args, **kwargs):
		raise q3c.RangeOverflowError("too many ipix ranges for the query")

	monkeypatch.setattr(QLSCIndex, "polygon_query", overflow)

	with caplog.at_level("WARNING", logger="qlsc"):
		matches = index.region_query(region)

	np.testing.assert_array_equal(sorted_pairs(matches), sorted_pairs(expected))
	assert any("overflow" in record.message for record in caplog.records)

@requires_cornish
def test_all_sky_uses_the_generic_path(index, frame):
	'''
	all_sky replaces the coarse step with a whole-index scan; the fast paths are about
	the coarse step, so they do not apply, and the answer is unchanged.
	'''
	region = ASTCircle(frame=frame, center=[45., 12.], radius=2.0)
	np.testing.assert_array_equal(sorted_pairs(index.region_query(region)),
								  sorted_pairs(index.region_query(region, all_sky=True)))

@requires_cornish
def test_a_polygon_wound_the_other_way_is_refused(index, frame):
	'''
	AST fills the region to the left of the vertex path, so reversing a polygon's
	vertices hands AST the rest of the sky - a region whose bounding circle bounds its
	complement. Q3C would happily answer about the small side, so this must be caught
	before any fast path runs, not after.

	Which of the two guards refuses it depends on the AST build: some report the
	reversed polygon as a plain polygon containing the rest of the sky (caught by the
	antipodal-point test), others as a negated region (starlink-pyast 4.2.0 does,
	caught by the isNegated test). Both refuse before the fast paths; either is right.
	'''
	region = ASTPolygon(frame=frame, points=[[44.,11.],[46.,11.],[46.,13.],[44.,13.]])
	assert region.pointInRegion([225., -12.]) and not region.pointInRegion([45., 12.]), \
		"the reversed polygon is expected to be the rest of the sky, not the small square"
	with pytest.raises(ValueError) as e:
		index.region_query(region)
	message = str(e.value)
	assert ("antipodal" in message) == (not region.isNegated)
	assert ("negated region cannot be searched" in message) == region.isNegated


# ---------------------------------------------------------------------------
# Result shape and flags, on the fast paths as well as the generic one.
# ---------------------------------------------------------------------------

@requires_cornish
@pytest.mark.parametrize("fast_path", [True, False])
def test_a_single_match_is_not_squeezed(keyed_index, frame, fast_path):
	'''
	region_query is MatchArrayShape.STABLE: one match is (1,2), not (2,), even when
	radial_query - which squeezes - supplied the answer.
	'''
	point = keyed_index.region_query(ASTCircle(frame=frame, center=[45., 12.], radius=30.0))[0]
	region = ASTCircle(frame=frame, center=list(point), radius=1e-4)
	matches = keyed_index.region_query(region, fast_path=fast_path)
	assert matches.shape == (1, 2)

@requires_cornish
@pytest.mark.parametrize("fast_path", [True, False])
def test_return_key_gives_the_same_records_either_way(keyed_index, frame, fast_path):
	'''The key field survives the fast path, as a recarray with the same three fields.'''
	region = ASTCircle(frame=frame, center=[45., 12.], radius=10.0)
	matches = keyed_index.region_query(region, return_key=True, fast_path=fast_path)
	assert matches.dtype.names == ("ra", "dec", "key")
	assert len(matches) > 0
	assert all(key.startswith("star-") for key in matches["key"])

@requires_cornish
def test_return_key_agrees_between_the_paths(keyed_index, frame):
	region = ASTPolygon(frame=frame, points=[[40.,20.],[50.,20.],[50.,10.],[40.,10.]])
	fast = keyed_index.region_query(region, return_key=True)
	generic = keyed_index.region_query(region, return_key=True, fast_path=False)
	assert len(fast) > 0
	assert sorted(fast["key"]) == sorted(generic["key"])

@requires_cornish
@requires_astropy
@pytest.mark.parametrize("fast_path", [True, False])
def test_return_skycoord_on_both_paths(keyed_index, frame, fast_path):
	'''
	The SkyCoord is built by the shared result helper, not inherited from whichever
	query supplied the matches, so it is an array SkyCoord on both paths.
	'''
	from astropy.coordinates import SkyCoord
	region = ASTCircle(frame=frame, center=[45., 12.], radius=10.0)
	matches = keyed_index.region_query(region, return_skycoord=True, fast_path=fast_path)
	assert isinstance(matches, SkyCoord)
	assert matches.isscalar is False
	assert len(matches) > 0

@requires_cornish
@pytest.mark.parametrize("fast_path", [True, False])
def test_no_matches_is_an_empty_array_of_pairs(index, frame, fast_path):
	'''An empty answer keeps the (0,2) shape rather than collapsing to something else.'''
	# a tiny circle in a gap between the seeded points
	point = index.region_query(ASTCircle(frame=frame, center=[45., 12.], radius=10.0))[0]
	region = ASTCircle(frame=frame, center=[point[0] + 0.2, point[1] + 0.2], radius=1e-6)
	matches = index.region_query(region, fast_path=fast_path)
	assert matches.shape == (0, 2)


# ---------------------------------------------------------------------------
# image_query
# ---------------------------------------------------------------------------

def tan_header(crval1:float, crval2:float, size:int=512, scale:float=0.001, radesys:str="ICRS"):
	'''
	Build a minimal gnomonic (TAN) FITS header for a square image.

	:param crval1: the right ascension of the reference pixel (degrees)
	:param crval2: the declination of the reference pixel (degrees)
	:param size: the image size in pixels, on both axes
	:param scale: the pixel scale in degrees
	:param radesys: the value of the RADESYS card
	'''
	header = fits.Header()
	header["NAXIS"] = 2
	header["NAXIS1"] = size
	header["NAXIS2"] = size
	header["CTYPE1"] = "RA---TAN"
	header["CTYPE2"] = "DEC--TAN"
	header["CRPIX1"] = size/2. + 0.5
	header["CRPIX2"] = size/2. + 0.5
	header["CRVAL1"] = crval1
	header["CRVAL2"] = crval2
	header["CDELT1"] = -scale
	header["CDELT2"] = scale
	header["RADESYS"] = radesys
	header["EQUINOX"] = 2000.0
	return header

@requires_cornish
@requires_astropy
@pytest.mark.parametrize("crval1,crval2", [
	(31., 44.),      # ordinary field
	(0., 0.),        # centred on ra=0
	(359.9, 0.1),    # straddling ra=0
	(120., 89.7),    # containing the north pole
	(300., -89.7),   # containing the south pole
])
def test_image_query_agrees_with_the_generic_path(frame, crval1, crval2):
	'''
	The sources on an image are found by building the field's region from the header and
	querying it, so the fast path and AST's own filter must agree here too. The field is
	made large enough (a degree across) that a 200k-point catalogue lands several sources
	on it.
	'''
	index = QLSCIndex(qlsc=QLSC(depth=14))
	try:
		index.add_points(points=sunflower_points_on_sphere(n=200000))
		header = tan_header(crval1, crval2, size=1024, scale=0.002)
		fast = index.image_query(header)
		generic = index.image_query(header, fast_path=False)
		assert len(fast) >= 1
		np.testing.assert_array_equal(sorted_pairs(fast), sorted_pairs(generic))
	finally:
		index.close()

@requires_cornish
@requires_astropy
def test_image_query_finds_the_sources_on_the_image_and_not_the_ones_beside_it():
	'''
	The question the method exists to answer, checked directly: points placed inside the
	field come back, points placed just outside it do not.
	'''
	header = tan_header(31., 44., size=512, scale=0.002)  # about a degree across
	index = QLSCIndex(qlsc=QLSC(depth=20))
	try:
		inside = [(31.0, 44.0), (31.2, 44.2), (30.8, 43.8)]
		outside = [(31.0, 46.0), (35.0, 44.0), (31.0, 42.0)]
		index.add_points(points=np.array(inside + outside),
						 keys=[f"in-{i}" for i in range(len(inside))] +
							  [f"out-{i}" for i in range(len(outside))])

		matches = index.image_query(header, return_key=True)
		assert sorted(matches["key"]) == sorted(f"in-{i}" for i in range(len(inside)))
	finally:
		index.close()

@requires_cornish
@requires_astropy
def test_image_query_refuses_a_frame_mismatch():
	'''
	A header declaring FK5 describes an FK5 field; qlsc stores bare numbers and cannot
	convert them, so the mismatch is refused rather than quietly mixing the two frames.
	'''
	index = QLSCIndex(qlsc=QLSC(depth=10))
	try:
		index.add_points(points=np.array([[31., 44.]]))
		with pytest.raises(ValueError) as e:
			index.image_query(tan_header(31., 44., radesys="FK5"))
		assert "FK5" in str(e.value)

		# ...and it is answerable by saying which frame the index is really in
		assert len(index.image_query(tan_header(31., 44., radesys="FK5"), frame="FK5")) == 1
	finally:
		index.close()

@requires_cornish
@requires_astropy
def test_image_query_passes_the_polygon_parameters_through():
	'''
	maxerr and maxvert reach cornish; a coarser maxerr gives a polygon with fewer
	vertices, which is the visible effect of the parameter arriving.
	'''
	header = tan_header(31., 44., size=2048, scale=0.001)
	fine = ASTRegion.fromFITSHeader(header, maxvert=200)
	coarse = ASTRegion.fromFITSHeader(header, maxvert=4)
	assert len(coarse.points) <= len(fine.points)

	index = QLSCIndex(qlsc=QLSC(depth=20))
	try:
		index.add_points(points=np.array([[31., 44.]]))
		assert len(index.image_query(header, maxvert=4)) == 1
	finally:
		index.close()

@requires_cornish
def test_image_query_refuses_a_missing_header():
	index = QLSCIndex(qlsc=QLSC(depth=5))
	try:
		with pytest.raises(ValueError) as e:
			index.image_query()
		assert "fits_header" in str(e.value)
	finally:
		index.close()

@requires_cornish
@requires_astropy
def test_image_query_method_matches_the_module_function():
	index = QLSCIndex(qlsc=QLSC(depth=14))
	try:
		index.add_points(points=sunflower_points_on_sphere(n=200000))
		header = tan_header(31., 44., size=1024, scale=0.002)
		np.testing.assert_array_equal(index.image_query(header),
									  cornish_interop.image_query(index, header))
	finally:
		index.close()

@requires_cornish
@requires_astropy
def test_image_query_forwards_all_sky():
	'''
	Both of region_query's refusals recommend all_sky=True, so image_query has to be
	able to take that advice. The field region is an ordinary bounded polygon, so the
	two answers are the same; what is checked is that the option reaches the query at
	all rather than being swallowed by the signature.
	'''
	index = QLSCIndex(qlsc=QLSC(depth=14))
	try:
		index.add_points(points=sunflower_points_on_sphere(n=200000))
		header = tan_header(31., 44., size=1024, scale=0.002)
		np.testing.assert_array_equal(sorted_pairs(index.image_query(header)),
									  sorted_pairs(index.image_query(header, all_sky=True)))
	finally:
		index.close()

@requires_cornish
def test_the_fast_path_returns_the_built_in_query_result_unchanged(index, frame):
	'''
	The fast path exists to avoid per-candidate Python work, so it returns what the
	built-in query built rather than taking it apart and reassembling it. Only
	radial_query's plain array needs anything done to it - its historical squeeze - and
	the shapes below are what says so.
	'''
	circle = ASTCircle(frame=frame, center=[45., 12.], radius=2.0)
	polygon = ASTPolygon(frame=frame, points=[[44.,13.],[46.,13.],[46.,11.],[44.,11.]])

	for region in (circle, polygon):
		plain = index.region_query(region)
		keyed = index.region_query(region, return_key=True)
		assert plain.ndim == 2 and plain.shape[1] == 2
		assert keyed.dtype.names == ("ra", "dec", "key")
		assert len(plain) == len(keyed)

@requires_cornish
def test_a_polygon_with_too_few_vertices_takes_the_generic_path(frame, monkeypatch):
	'''
	Q3C's polygon routines take at least three vertices, and polygon_query raises below
	that. AST will not build a polygon with fewer, so the guard in the dispatch cannot
	be reached through the public API - which is exactly why it is worth having, and why
	the branch is entered here directly with the vertex list replaced, rather than
	asserted about from the outside and left unexecuted.
	'''
	assert cornish_interop.MINIMUM_POLYGON_VERTICES == 3

	index = QLSCIndex(qlsc=QLSC(depth=8))
	try:
		# what polygon_query does with two vertices, which is what the guard prevents
		with pytest.raises(ValueError) as e:
			index.polygon_query(polygon=np.array([[10., 10.], [11., 10.]]))
		assert "3" in str(e.value)

		region = ASTPolygon(frame=frame, points=[[44.,13.],[46.,13.],[46.,11.],[44.,11.]])
		monkeypatch.setattr(type(region), "points",
							property(lambda self: np.array([[44., 13.], [46., 13.]])))
		assert cornish_interop._fast_path_matches(index, region, 1.0, False, False) \
			is cornish_interop._NO_FAST_PATH
	finally:
		index.close()
