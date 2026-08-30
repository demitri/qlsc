
'''
Tests for :py:mod:`qlsc.cornish_interop`, the bridge between qlsc's index and
cornish/Starlink AST regions.

Most of these tests need the optional ``cornish`` package. When it is absent they
are SKIPPED, not silently passed, and the skip reason names the missing package —
run pytest with ``-rs`` to see it. The tests that guard the missing-package
behaviour itself run unconditionally.
'''

import numpy as np
import pytest

import qlsc
from qlsc import QLSC, QLSCIndex
from qlsc.generate import sunflower_points_on_sphere
from qlsc import cornish_interop

try:
	import cornish
	from cornish import (ASTCircle, ASTBox, ASTPolygon, ASTMoc, ASTICRSFrame,
						 ASTSkyFrame, ASTCompoundRegion, RegionOperation)
	cornish_available = True
except ImportError:
	cornish_available = False

CORNISH_SKIP_REASON = ("the optional 'cornish' package is not installed; "
					   "qlsc.cornish_interop cannot be exercised without it "
					   "(pip install cornish)")

requires_cornish = pytest.mark.skipif(not cornish_available, reason=CORNISH_SKIP_REASON)

# ---------------------------------------------------------------------------
# behaviour when cornish is NOT installed (runs whether or not it is)
# ---------------------------------------------------------------------------

def test_missing_cornish_raises_import_error(monkeypatch):
	'''
	Every public function must fail with a descriptive ImportError, never quietly,
	when cornish is unavailable.
	'''
	original = ImportError("No module named 'cornish'")
	monkeypatch.setattr(cornish_interop, "cornish_available", False)
	monkeypatch.setattr(cornish_interop, "cornish_import_error", original)

	q = QLSC(depth=4)
	index = QLSCIndex(qlsc=q)

	for call in (lambda: cornish_interop.ipix2region(q, 0),
				 lambda: cornish_interop.ipix2moc(q, [0]),
				 lambda: cornish_interop.region_query(index, object())):
		with pytest.raises(ImportError) as excinfo:
			call()
		assert "cornish" in str(excinfo.value)
		assert "could not be imported" in str(excinfo.value)
		# the underlying reason must survive: "not installed" and "installed but broken"
		# are different problems and must not be reported as the same one
		assert excinfo.value.__cause__ is original

def test_module_imports_without_cornish():
	'''
	The module itself must import cleanly regardless; only calls fail.
	'''
	assert hasattr(cornish_interop, "cornish_available")
	assert isinstance(cornish_interop.cornish_available, bool)

# ---------------------------------------------------------------------------
# ipix -> AST polygon
# ---------------------------------------------------------------------------

@requires_cornish
def test_ipix2region_contains_pixel_center():
	'''
	The polygon for a pixel must contain that pixel's center — the check that the
	vertex winding was translated to AST's convention.
	'''
	q = QLSC(depth=5)
	for ra, dec in [(12.0, 5.0), (300.0, -80.0), (0.05, 0.0), (45.0, 89.9),
					(180.0, -45.0), (359.99, 0.5), (90.0, -89.99)]:
		ipix = int(q.ang2ipix(ra, dec))
		polygon = cornish_interop.ipix2region(q, ipix)
		center = q.ipix2ang_center(ipix)
		assert polygon.pointInRegion([center[0], center[1]]) is True
		assert polygon.points.shape == (4, 2)

@requires_cornish
def test_ipix2region_agrees_with_point_in_polygon():
	'''
	AST's containment test and Q3C's must agree on the same pixel: both define
	polygon edges as great circles, which is what makes the conversion exact.
	'''
	q = QLSC(depth=5)
	rng = np.random.default_rng(20260830)
	tested = 0

	for ipix in [int(q.ang2ipix(ra, dec)) for ra, dec in
				 [(12.0, 5.0), (128.0, 33.0), (222.0, -17.0), (300.0, -80.0), (45.0, 88.0)]]:
		vertices = q.ipix2polygon(ipix)
		polygon = cornish_interop.ipix2region(q, ipix)
		disc = polygon.boundingCircle()
		center = disc.centre
		radius = disc.radius.to("deg").value

		for _ in range(60):
			angle = rng.uniform(0, 2 * np.pi)
			offset = radius * np.sqrt(rng.uniform())
			dec = center[1] + offset * np.sin(angle)
			if abs(dec) > 89.5:
				continue # the local flat-sky offset below is not valid at the pole
			ra = (center[0] + offset * np.cos(angle) / np.cos(np.radians(center[1]))) % 360.0
			assert polygon.pointInRegion([ra, dec]) == q.point_in_polygon(ra, dec, vertices)
			tested += 1

	assert tested > 200 # the loop actually ran

@requires_cornish
def test_ipix2region_rejects_wrong_types():
	q = QLSC(depth=4)
	with pytest.raises(TypeError):
		cornish_interop.ipix2region("not a QLSC", 0)

# ---------------------------------------------------------------------------
# ipix set -> MOC
# ---------------------------------------------------------------------------

@requires_cornish
def test_ipix2moc_covers_the_pixels():
	'''
	A MOC built from a set of pixels contains their centers and has an area of
	the same order as the pixels' own.
	'''
	q = QLSC(depth=5)
	ipix = [int(q.ang2ipix(ra, 5.0)) for ra in (10.0, 12.0, 14.0, 16.0)]
	ipix = sorted(set(ipix))

	moc = cornish_interop.ipix2moc(q, ipix)
	assert isinstance(moc, ASTMoc)
	assert moc.maxOrder == q.depth + 5

	for value in ipix:
		center = q.ipix2ang_center(value)
		assert moc.pointInRegion([center[0], center[1]]) is True

	# the sphere is 41253 sq deg over 6*4**depth pixels; the MOC should be close to
	# that many pixels' worth of area (QLSC pixels are only approximately equal-area,
	# so this is an order-of-magnitude check, not an equality)
	nominal = 41253.0 / (6 * 4 ** q.depth) * len(ipix)
	assert 0.3 * nominal < moc.area.to("deg2").value < 3.0 * nominal

@requires_cornish
def test_ipix2moc_refuses_a_max_order_that_loses_coverage():
	'''
	AST returns an EMPTY MOC, without complaint, when max_order is too coarse for the
	region being added; ipix2moc must turn that into a loud failure. (A depth-5 pixel
	is silently dropped at order 7 on AST 9.3.)
	'''
	q = QLSC(depth=5)
	ipix = [int(q.ang2ipix(10.0, 5.0))]

	with pytest.raises(ValueError) as excinfo:
		cornish_interop.ipix2moc(q, ipix, max_order=7)
	assert "max_order" in str(excinfo.value)

@requires_cornish
def test_ipix2moc_refuses_order_zero():
	'''
	Adding a region to an order-0 MOC segfaults AST 9.3; the guard must run before AST does.
	'''
	q = QLSC(depth=5)
	with pytest.raises(ValueError) as excinfo:
		cornish_interop.ipix2moc(q, [int(q.ang2ipix(10.0, 5.0))], max_order=0)
	assert "SIGSEGV" in str(excinfo.value)

@requires_cornish
def test_ipix2moc_round_trips_through_the_ivoa_string():
	q = QLSC(depth=5)
	ipix = [int(q.ang2ipix(ra, 5.0)) for ra in (10.0, 12.0)]
	moc = cornish_interop.ipix2moc(q, ipix)

	restored = ASTMoc.fromString(moc.toString())
	assert restored.area.to("deg2").value == pytest.approx(moc.area.to("deg2").value, rel=1e-9)

@requires_cornish
def test_ipix2moc_rejects_an_empty_set():
	q = QLSC(depth=4)
	with pytest.raises(ValueError):
		cornish_interop.ipix2moc(q, [])

@requires_cornish
@pytest.mark.parametrize("bad", [1.9, 1.0, True, "3", np.float64(2.5)])
def test_ipix2moc_refuses_a_pixel_number_that_is_not_an_integer(bad):
	'''
	A value that is not an integral pixel number must raise, not be truncated to the
	pixel next door.

	int() used to do the truncating here, so ipix2moc(q, [1.9]) built a MOC covering
	pixel 1 and said nothing — a map of somewhere the caller did not ask about. Every
	other entry point in qlsc refuses these, and the bar is the same one ipix2ang sets.
	'''
	q = QLSC(depth=4)

	with pytest.raises(ValueError):
		cornish_interop.ipix2moc(q, [bad])

	# the same bar the core API applies
	with pytest.raises(ValueError):
		q.ipix2ang(bad)

@requires_cornish
@pytest.mark.parametrize("bad", [-1, 6144])
def test_ipix2moc_refuses_a_pixel_number_out_of_range(bad):
	'''
	A pixel number outside the scheme raises, and it raises before any region is built:
	a bad value at the end of a list must not be found only after the work is done.
	'''
	q = QLSC(depth=5) # 6144 pixels, so 6144 is one past the end

	with pytest.raises(ValueError):
		cornish_interop.ipix2moc(q, [bad])
	with pytest.raises(ValueError):
		cornish_interop.ipix2moc(q, [0, 1, bad])

@requires_cornish
def test_ipix2region_refuses_the_same_values():
	'''
	The single-pixel conversion refuses what the MOC conversion refuses; it reaches the
	same validator through ipix2polygon, and this pins that it keeps doing so.
	'''
	q = QLSC(depth=4)
	for bad in (1.9, True, -1, q.nbins):
		with pytest.raises(ValueError):
			cornish_interop.ipix2region(q, bad)

@requires_cornish
def test_ipix2moc_accepts_integer_types():
	'''
	Integers of any width are accepted, as they are everywhere else in qlsc; only the
	non-integral values are refused.
	'''
	q = QLSC(depth=4)
	for good in (3, np.int32(3), np.int64(3), np.uint16(3)):
		moc = cornish_interop.ipix2moc(q, [good])
		assert moc.pointInRegion(list(q.ipix2ang_center(3)))

	# and a NumPy array of pixel numbers, the natural output of a coverage query
	moc = cornish_interop.ipix2moc(q, np.array([3, 4], dtype=np.int64))
	assert moc.pointInRegion(list(q.ipix2ang_center(4)))

# ---------------------------------------------------------------------------
# region_query
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def index_with_points():
	'''
	An in-memory index of evenly distributed points. Module-scoped: building it is
	the slow part of these tests, and no test modifies it.
	'''
	q = QLSC(depth=6)
	index = QLSCIndex(qlsc=q)
	points = sunflower_points_on_sphere(50000)
	index.add_points(ra=points[:, 0], dec=points[:, 1])
	return index, points

@requires_cornish
def test_region_query_circle_equals_radial_query(index_with_points):
	'''
	An ASTCircle must return precisely what the equivalent radial_query returns.
	'''
	index, _ = index_with_points
	circle = ASTCircle(frame=ASTICRSFrame(), center=[12.0, 5.0], radius=2.0)

	from_region = cornish_interop.region_query(index, circle)
	from_cone = np.atleast_2d(index.radial_query(12.0, 5.0, 2.0))

	assert from_region.shape == from_cone.shape
	assert from_region.shape[0] > 10 # a meaningful number of matches
	assert np.allclose(np.sort(from_region, axis=0), np.sort(from_cone, axis=0))

@requires_cornish
def test_region_query_box_matches_brute_force(index_with_points):
	'''
	A box has no qlsc equivalent; check against AST's own answer over every point.
	'''
	index, points = index_with_points
	box = ASTBox.fromCentreAndCorner(frame=ASTICRSFrame(), centre=[12.0, 5.0], corner=[14.0, 7.0])

	matches = cornish_interop.region_query(index, box)
	expected = np.array([p for p in points if box.pointInRegion([p[0], p[1]])])

	assert matches.shape[0] == expected.shape[0]
	assert matches.shape[0] > 10
	assert np.allclose(np.sort(matches, axis=0), np.sort(expected, axis=0))

@requires_cornish
def test_region_query_compound_union(index_with_points):
	'''
	The union of two disjoint regions must return the union of their matches —
	the case that needs AST, since qlsc has no compound region.
	'''
	index, _ = index_with_points
	frame = ASTICRSFrame()
	circle = ASTCircle(frame=frame, center=[12.0, 5.0], radius=2.0)
	box = ASTBox.fromCentreAndCorner(frame=frame, centre=[40.0, 5.0], corner=[42.0, 7.0])
	union = ASTCompoundRegion(regions=[circle, box], operation=RegionOperation.OR)

	from_union = cornish_interop.region_query(index, union)
	from_parts = np.vstack([cornish_interop.region_query(index, circle),
							cornish_interop.region_query(index, box)])

	assert from_union.shape == from_parts.shape
	assert from_union.shape[0] > 20
	assert np.allclose(np.sort(from_union, axis=0), np.sort(from_parts, axis=0))

@requires_cornish
def test_region_query_polygon_equals_polygon_query(index_with_points):
	'''
	An ASTPolygon must return what QLSCIndex.polygon_query returns for the same
	vertices — the two libraries' polygons are the same shape.
	'''
	index, _ = index_with_points
	vertices = np.array([[10.0, 4.0], [14.0, 4.0], [14.0, 8.0], [10.0, 8.0]])
	# AST fills the region to the LEFT of the vertex path, hence the reversal
	polygon = ASTPolygon(frame=ASTICRSFrame(), points=vertices[::-1])

	from_region = cornish_interop.region_query(index, polygon)
	from_qlsc = np.atleast_2d(index.polygon_query(vertices))

	assert from_region.shape == from_qlsc.shape
	assert from_region.shape[0] > 10
	assert np.allclose(np.sort(from_region, axis=0), np.sort(from_qlsc, axis=0))

@requires_cornish
def test_region_query_refuses_a_negated_region(index_with_points):
	'''
	A negated region's bounding circle is the circle of the un-negated region, which
	would silently return nothing; the query must refuse it and say what to do instead.
	'''
	index, _ = index_with_points
	circle = ASTCircle(frame=ASTICRSFrame(), center=[12.0, 5.0], radius=2.0)
	circle.negate()

	with pytest.raises(ValueError) as excinfo:
		cornish_interop.region_query(index, circle)
	assert "all_sky=True" in str(excinfo.value)

@requires_cornish
def test_region_query_all_sky_handles_a_negated_region():
	'''
	With all_sky=True a negated region is answered correctly: everything outside the circle.
	'''
	q = QLSC(depth=5)
	index = QLSCIndex(qlsc=q)
	points = sunflower_points_on_sphere(2000)
	index.add_points(ra=points[:, 0], dec=points[:, 1])

	circle = ASTCircle(frame=ASTICRSFrame(), center=[12.0, 5.0], radius=20.0)
	inside = cornish_interop.region_query(index, circle)

	circle.negate()
	outside = cornish_interop.region_query(index, circle, all_sky=True)

	assert inside.shape[0] > 0
	assert inside.shape[0] + outside.shape[0] == points.shape[0]

@requires_cornish
def test_region_query_refuses_a_compound_region_its_disc_does_not_bound():
	'''
	A compound region with a negated component reports isNegated == False, so the
	antipodal falsification probe is what catches it.
	'''
	q = QLSC(depth=5)
	index = QLSCIndex(qlsc=q)
	points = sunflower_points_on_sphere(2000)
	index.add_points(ra=points[:, 0], dec=points[:, 1])

	frame = ASTICRSFrame()
	small = ASTCircle(frame=frame, center=[12.0, 5.0], radius=1.0)
	small.negate()
	other = ASTCircle(frame=frame, center=[200.0, -40.0], radius=1.0)
	# NOT(small) OR other covers nearly the whole sky
	union = ASTCompoundRegion(regions=[small, other], operation=RegionOperation.OR)
	assert union.isNegated is False

	with pytest.raises(ValueError) as excinfo:
		cornish_interop.region_query(index, union)
	assert "antipodal" in str(excinfo.value)

@requires_cornish
def test_region_query_refuses_a_frame_mismatch(index_with_points):
	'''
	A galactic region against an index declared to hold ICRS coordinates must fail loudly.
	'''
	index, _ = index_with_points
	galactic = ASTCircle(frame=ASTSkyFrame(system="Galactic"), center=[12.0, 5.0], radius=1.0)

	with pytest.raises(ValueError) as excinfo:
		cornish_interop.region_query(index, galactic)
	assert "GALACTIC" in str(excinfo.value).upper()

	# ... and must be allowed when the caller says the index is galactic too
	matches = cornish_interop.region_query(index, galactic, frame="Galactic")
	assert matches.ndim == 2

	# ... or when the check is explicitly waived
	waived = cornish_interop.region_query(index, galactic, skip_frame_check=True)
	assert waived.shape == matches.shape

@requires_cornish
def test_region_query_frame_none_means_the_default_not_no_check(index_with_points):
	'''
	frame=None must mean "the default frame (ICRS)", never "skip the check" — the two
	are spelled differently on purpose, so that omitting a frame cannot silently waive
	the comparison.
	'''
	index, _ = index_with_points
	galactic = ASTCircle(frame=ASTSkyFrame(system="Galactic"), center=[12.0, 5.0], radius=1.0)

	with pytest.raises(ValueError):
		cornish_interop.region_query(index, galactic, frame=None)

	icrs = ASTCircle(frame=ASTICRSFrame(), center=[12.0, 5.0], radius=1.0)
	explicit = cornish_interop.region_query(index, icrs, frame="ICRS")
	defaulted = cornish_interop.region_query(index, icrs, frame=None)
	assert np.array_equal(np.atleast_2d(explicit), np.atleast_2d(defaulted))

@requires_cornish
def test_region_query_rejects_a_non_string_frame(index_with_points):
	index, _ = index_with_points
	circle = ASTCircle(frame=ASTICRSFrame(), center=[12.0, 5.0], radius=1.0)
	with pytest.raises(TypeError):
		cornish_interop.region_query(index, circle, frame=1950)

@requires_cornish
def test_region_query_rejects_wrong_types(index_with_points):
	index, _ = index_with_points
	with pytest.raises(TypeError):
		cornish_interop.region_query(index, "not a region")
	with pytest.raises(TypeError):
		cornish_interop.region_query("not an index", ASTCircle(frame=ASTICRSFrame(), center=[0.0, 0.0], radius=1.0))

@requires_cornish
def test_region_query_return_key():
	'''
	The key column comes back as a recarray, in the same form the other region queries use.
	'''
	q = QLSC(depth=5)
	index = QLSCIndex(qlsc=q)
	index.add_points(ra=[12.0, 12.1, 200.0], dec=[5.0, 5.1, -40.0], keys=["a", "b", "c"])

	circle = ASTCircle(frame=ASTICRSFrame(), center=[12.0, 5.0], radius=1.0)
	matches = cornish_interop.region_query(index, circle, return_key=True)

	assert sorted(matches["key"].tolist()) == ["a", "b"]
	assert set(matches.dtype.names) == {"ra", "dec", "key"}

@requires_cornish
def test_region_query_single_match_is_not_squeezed():
	'''
	region_query is MatchArrayShape.STABLE, like every query added since radial_query:
	a single match is (1,2). radial_query supplies the candidates but its historical
	squeeze to (2,) must not leak through.
	'''
	q = QLSC(depth=5)
	index = QLSCIndex(qlsc=q)
	index.add_points(ra=[12.0, 200.0], dec=[5.0, -40.0])

	circle = ASTCircle(frame=ASTICRSFrame(), center=[12.0, 5.0], radius=0.5)
	matches = cornish_interop.region_query(index, circle)

	assert matches.shape == (1, 2)
	assert matches[0][0] == pytest.approx(12.0)
	assert matches[:, 0].shape == (1,) # the subscript that fails on a squeezed result

	# ... and the coarse step really is the squeezed one, so this is not a coincidence
	assert index.radial_query(12.0, 5.0, 0.5).shape == (2,)

@requires_cornish
def test_region_query_return_skycoord():
	'''
	Parity with the other queries after first-class astropy support landed.
	'''
	# astropy is a soft dependency of qlsc (cornish requires it, but this file's
	# always-run tests must not), so it is skipped for here rather than imported above
	SkyCoord = pytest.importorskip("astropy.coordinates").SkyCoord
	QTable = pytest.importorskip("astropy.table").QTable

	q = QLSC(depth=5)
	index = QLSCIndex(qlsc=q)
	index.add_points(ra=[12.0, 12.1, 200.0], dec=[5.0, 5.1, -40.0], keys=["a", "b", "c"])
	circle = ASTCircle(frame=ASTICRSFrame(), center=[12.0, 5.0], radius=1.0)

	coords = cornish_interop.region_query(index, circle, return_skycoord=True)
	assert isinstance(coords, SkyCoord)
	assert len(coords) == 2
	assert sorted(coords.ra.deg) == pytest.approx([12.0, 12.1])

	# with a flag as well, the coordinates come back in a QTable alongside it
	table = cornish_interop.region_query(index, circle, return_key=True, return_skycoord=True)
	assert isinstance(table, QTable)
	assert sorted(table["key"].tolist()) == ["a", "b"]
	assert isinstance(table["coord"], SkyCoord)

@requires_cornish
def test_region_query_with_no_matches():
	q = QLSC(depth=5)
	index = QLSCIndex(qlsc=q)
	index.add_points(ra=[200.0], dec=[-40.0])

	circle = ASTCircle(frame=ASTICRSFrame(), center=[12.0, 5.0], radius=0.5)
	matches = cornish_interop.region_query(index, circle)

	assert matches.shape[0] == 0

@requires_cornish
def test_region_query_on_a_moc(index_with_points):
	'''
	A MOC is a region like any other — including one built from qlsc's own pixels,
	which closes the loop: pixels out, region back in.
	'''
	index, _ = index_with_points
	q = index.qlsc
	ipix = sorted({int(q.ang2ipix(ra, 5.0)) for ra in (11.0, 11.5, 12.0, 12.5)})
	moc = cornish_interop.ipix2moc(q, ipix)

	matches = cornish_interop.region_query(index, moc)
	assert matches.shape[0] > 0

	# every match must be inside one of the originating pixels (up to the MOC's resolution)
	for ra, dec in matches:
		assert moc.pointInRegion([ra, dec]) is True
