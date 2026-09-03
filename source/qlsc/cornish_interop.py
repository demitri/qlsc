
'''
Interoperation between qlsc and `cornish <https://github.com/demitri/cornish>`_,
a Pythonic interface to the `Starlink AST <https://starlink.eao.hawaii.edu/starlink/AST>`_
library.

The two libraries divide the work cleanly:

* **AST defines regions.** Circles, boxes, polygons, ellipses, compound
  (boolean) regions, MOCs — all with exact containment, boolean algebra,
  serialization (STC-S, MOC, FITS) and plotting.
* **qlsc defines an index.** Which stored points *might* be in a region, found
  from z-order ipix ranges without touching the rest of the catalogue.

This module joins them with the same two-step shape every :class:`~qlsc.QLSCIndex`
region query already uses internally: a **coarse step that is guaranteed to be a
superset**, followed by an **exact per-candidate test**. Here the coarse step is
a cone search over the region's bounding circle and the exact test is AST's own
``pointInRegion``, so the result is exactly the set of indexed points that AST
considers inside the region — for any region AST can describe.

Where qlsc already defines the same shape AST does — a circle is a cone, a polygon
has great-circle edges — the whole query is handed to the matching built-in method
instead, so the exact test runs as one batched call into the C layer rather than one
``pointInRegion`` call per candidate. Which regions those are, and why an AST box is
not one of them, is documented on :py:func:`region_query`.

The two entry points most users want are methods on the index itself:
:py:func:`QLSCIndex.region_query` and :py:func:`QLSCIndex.image_query`, the latter
answering "which of my sources land on this image?" from a FITS header alone. Both
delegate to this module.

``cornish`` is an optional dependency; it is not required to use qlsc. Importing
this module without it installed succeeds, but calling any function in it raises
an :class:`ImportError` naming the missing package.

Coordinate frames
-----------------

qlsc stores bare (ra,dec) numbers and has no concept of a coordinate frame; an
AST region always carries one. Every function here therefore compares the
region's frame against the frame the caller declares the *index* to be in (the
``frame`` parameter) and refuses a mismatch rather than silently mixing
coordinate systems. Convert a region between frames with
``ASTRegion.regionWithMapping(frame=...)``.

``frame=None`` means "the default frame" (ICRS) everywhere in this module, never
"do not check" — waiving the comparison is a separate, explicit
``skip_frame_check=True``, so that forgetting to pass a frame and deliberately
taking responsibility for the frames matching cannot be spelled the same way.
'''

import logging
from typing import Iterable, Union

import numpy as np

try:
	import cornish
	from cornish import ASTRegion, ASTCircle, ASTPolygon, ASTMoc, ASTICRSFrame
	cornish_available = True
	cornish_import_error = None
except ImportError as e:
	cornish_available = False
	# Keep the original error: "cornish is not installed" and "cornish is installed but its
	# own starlink-pyast import failed" are different problems with the same symptom here,
	# and reporting the first when it is really the second sends the user the wrong way.
	cornish_import_error = e

from . import q3c
from .qlsc_base import QLSC, QLSCIndex, MatchArrayShape, _match_result, _as_stable_pairs
from .utilities import _as_float_scalar

logger = logging.getLogger("qlsc")

# The frame the index coordinates are assumed to be in when the caller does not say.
DEFAULT_FRAME = "ICRS"

#: The largest number of vertices a polygon may have and still be answered by
#: :py:func:`QLSCIndex.polygon_query`. Read from the extension module rather than
#: repeated here: the limit is q3c's own fixed buffer length, and two copies of a
#: number are two chances to disagree.
MAXIMUM_POLYGON_VERTICES = q3c.Q3C_MAX_N_POLY_VERTEX

#: The fewest vertices a polygon may have and still be a polygon. AST does not build
#: one with fewer, so this bound is latent; the dispatch checks it anyway rather than
#: relying on another library's constructor to hold it.
MINIMUM_POLYGON_VERTICES = 3

#: The largest bounding-circle radius, in degrees, for which an
#: :class:`~cornish.ASTPolygon` is answered with :py:func:`QLSCIndex.polygon_query`
#: rather than the generic path.
#:
#: Q3C projects a polygon onto a single cube face and refuses one that will not fit;
#: where that happens depends on the polygon's position (measured between 43 and 76
#: degrees across for a square at five positions - see
#: ``docs/examples/survey_footprint.rst``), and the API documentation quotes roughly 23
#: degrees across as the figure that is safe anywhere. A polygon inside a circle of this
#: radius is at most twice it across, comfortably inside that. Going over it costs
#: nothing but speed: the region takes the generic path, which has no size limit.
POLYGON_FAST_PATH_MAXIMUM_RADIUS = 10.0

def _require_cornish():
	'''
	Raise a descriptive :class:`ImportError` when the optional ``cornish`` package could
	not be imported. The original import error is chained, since a broken ``starlink-pyast``
	underneath an installed cornish looks identical from here otherwise.
	'''
	if not cornish_available:
		raise ImportError(
			"The 'cornish' package is required for qlsc.cornish_interop but could not be "
			f"imported ({cornish_import_error}). Install it with 'pip install cornish' (it in "
			"turn requires 'starlink-pyast'), or see https://github.com/demitri/cornish ."
			) from cornish_import_error

def _check_region(region, frame:Union[str,None], skip_frame_check:bool) -> None:
	'''
	Validate that an object is a sky-frame :class:`cornish.ASTRegion` in the expected frame.

	:param region: the object to check
	:param frame: name of the frame the index coordinates are in; None means the default (ICRS)
	:param skip_frame_check: if True, the frame comparison is waived (the type and sky-frame checks are not)
	'''
	if not isinstance(region, ASTRegion):
		raise TypeError(f"A cornish ASTRegion was expected; was given '{type(region).__name__}'.")

	region_frame = region.frame()
	if not region_frame.isSkyFrame:
		raise ValueError("The region must be defined in a sky frame; this region is in a "
						 f"'{region_frame.system}' frame, whose coordinates are not (ra,dec).")

	if skip_frame_check:
		return # the caller has explicitly taken responsibility for the frames matching

	if frame is None:
		frame = DEFAULT_FRAME

	if not isinstance(frame, str):
		raise TypeError(f"The 'frame' parameter must be a frame name (a string) or None for the "
						f"default ('{DEFAULT_FRAME}'); was given '{type(frame).__name__}'.")

	if region_frame.system.upper() != frame.upper():
		raise ValueError(
			f"The region is in the '{region_frame.system}' frame but the index coordinates were "
			f"declared to be in '{frame}'. qlsc stores bare coordinates and cannot convert them; "
			f"either convert the region first (region.regionWithMapping(frame=cornish.ASTSkyFrame(system='{frame}'))), "
			f"pass frame='{region_frame.system}' if the index really holds {region_frame.system} coordinates, "
			f"or pass skip_frame_check=True to waive this check.")

def ipix2region(qlsc:QLSC, ipix:int, frame=None) -> "ASTPolygon":
	'''
	Return the pixel with the given ipix number as a cornish :class:`~cornish.ASTPolygon`.

	The polygon has the pixel's four corners as its vertices and great-circle edges,
	which is exactly how both libraries define a polygon on the sphere — the two
	containment tests (:py:func:`QLSC.point_in_polygon` and AST's ``pointInRegion``)
	agree on the resulting shape.

	The one thing that does *not* carry over is vertex order: AST takes the interior
	of a polygon to be on the left of the path traced by the vertices, so the corners
	from :py:func:`QLSC.ipix2polygon` describe the *rest of the sky* if handed to AST
	unchanged. They are reversed here, and the result is verified by confirming that
	the pixel's own center lies inside the polygon — if that ever fails the winding
	convention of one of the two libraries has changed, and a RuntimeError is raised
	rather than a silently inverted region returned.

	:param qlsc: the :class:`~qlsc.QLSC` scheme the ipix value belongs to
	:param ipix: the ipix number of the pixel
	:param frame: the sky frame to define the polygon in; None means an ICRS frame
	:returns: a :class:`~cornish.ASTPolygon` covering the pixel
	'''
	polygon, _ = _ipix2region_with_center(qlsc, ipix, frame=frame)
	return polygon

def _ipix2region_with_center(qlsc:QLSC, ipix:int, frame=None):
	'''
	The body of :py:func:`ipix2region`, which also hands back the pixel center it had to
	compute for the winding check so that callers needing it do not compute it twice.

	:param qlsc: the :class:`~qlsc.QLSC` scheme the ipix value belongs to
	:param ipix: the ipix number of the pixel
	:param frame: the sky frame to define the polygon in; None means an ICRS frame
	:returns: a tuple of (:class:`~cornish.ASTPolygon`, the pixel's (ra,dec) center)
	'''
	_require_cornish()

	if not isinstance(qlsc, QLSC):
		raise TypeError(f"A qlsc.QLSC object was expected for 'qlsc'; was given '{type(qlsc).__name__}'.")

	if frame is None:
		frame = ASTICRSFrame()

	vertices = qlsc.ipix2polygon(ipix)          # (4,2) [ra,dec] pairs in degrees
	polygon = ASTPolygon(frame=frame, points=vertices[::-1]) # see the winding note above

	center = qlsc.ipix2ang_center(ipix)
	if not polygon.pointInRegion([center[0], center[1]]):
		raise RuntimeError(
			f"The polygon built for ipix {ipix} does not contain the pixel's own center "
			f"{center}: the vertex winding convention of qlsc or of AST has changed and "
			"this conversion is no longer correct. Please report this as a qlsc bug.")

	return polygon, center

def ipix2moc(qlsc:QLSC, ipix:Iterable[int], max_order:int=None, frame=None) -> "ASTMoc":
	'''
	Return a set of qlsc pixels as a cornish :class:`~cornish.ASTMoc` (an IVOA
	Multi-Order Coverage map).

	This is the bridge from a qlsc coverage — the distinct ipix values of a catalogue,
	the pixels touched by a query, a survey footprint — to the standard VO coverage
	object: a MOC can be intersected and unioned with other MOCs, asked for its area,
	tested for containment, plotted, and serialized to the IVOA string that Aladin,
	MOCpy, VizieR and pgsphere all read.

	A MOC is built from HEALPix cells, so it approximates the QLSC pixel boundaries to
	the resolution of its finest cells. ``max_order`` defaults to the qlsc scheme's
	``depth`` + 5 (capped at AST's maximum of 27), a margin chosen empirically: AST's
	region-to-MOC conversion returns an *empty* MOC, without complaint, when the order
	is too coarse for the region it is given, and the threshold at which that happens
	is not a simple function of the pixel size. Every pixel added here is therefore
	checked against the finished MOC, and a ValueError naming ``max_order`` is raised
	if the coverage came back empty — this function never returns a MOC that quietly
	omits pixels it was asked to include.

	Two AST limits are worth knowing:

	* ``max_order`` of 0 is refused: AST 9.3 **crashes the process** (SIGSEGV) when a
	  region is added to an order-0 MOC.
	* Very coarse schemes have pixels tens of degrees across, and AST refuses to map
	  some of them (particularly those containing a pole) onto HEALPix at all; it
	  raises ``Ast.NODEF``, which is left to propagate.

	The ipix values are validated as :py:func:`QLSC.ipix2ang` validates them: integers of
	any width, in range for the scheme, with a non-integral or out-of-range value raising
	``ValueError`` rather than being rounded to the pixel next door.

	:param qlsc: the :class:`~qlsc.QLSC` scheme the ipix values belong to
	:param ipix: an iterable of ipix numbers
	:param max_order: the HEALPix order of the finest cells in the MOC, 1..27
	:param frame: the sky frame to define the pixels in; None means an ICRS frame
	:returns: a :class:`~cornish.ASTMoc` covering the given pixels
	'''
	_require_cornish()

	if not isinstance(qlsc, QLSC):
		raise TypeError(f"A qlsc.QLSC object was expected for 'qlsc'; was given '{type(qlsc).__name__}'.")

	if max_order is None:
		max_order = min(qlsc.depth + 5, 27)
	elif max_order == 0:
		raise ValueError("A 'max_order' of 0 is refused: AST 9.3 crashes the process (SIGSEGV) "
						 "when a region is added to an order-0 MOC. Use 1 or greater.")

	if frame is None:
		frame = ASTICRSFrame()

	# The scheme's own validator, not int(): int() TRUNCATES, so an ipix value of 1.9
	# quietly became pixel 1, True became pixel 1, and the string "3" became pixel 3 -
	# a MOC covering somewhere the caller did not ask about, built without a word. Every
	# other entry point in qlsc refuses a non-integral or out-of-range pixel number
	# rather than rounding it to a neighbour, and this one now refuses it with the same
	# check and the same message. The whole iterable is validated before any region is
	# built, so a bad value halfway down a list is not found after an arbitrary amount
	# of work.
	pixels = qlsc._validated_ipix(list(ipix))

	moc = ASTMoc(max_order=max_order)
	centers = list()
	for value in np.atleast_1d(pixels):
		value = int(value)
		# the center comes back from the conversion, which needed it for its winding check
		region, center = _ipix2region_with_center(qlsc, value, frame=frame)
		moc.add(region)
		centers.append((value, center))

	if len(centers) == 0:
		raise ValueError("No ipix values were provided; an empty MOC has no meaningful coverage.")

	# AST silently produces an empty MOC when max_order is too coarse for the regions
	# being added, so confirm the coverage is real rather than trusting the call.
	for value, center in centers:
		if not moc.pointInRegion([center[0], center[1]]):
			raise ValueError(
				f"The MOC built at max_order={max_order} does not cover pixel {value} (center "
				f"{center}). AST's region-to-MOC conversion returns an empty or incomplete map "
				f"when the order is too coarse for the region; raise 'max_order' and try again.")

	return moc

#: Returned by :py:func:`_fast_path_matches` when the region has no fast path, so that
#: "there is no fast path for this region" cannot be confused with a query result. No
#: query returns this object, and nothing outside this module sees it.
_NO_FAST_PATH = object()

def _candidate_rows(candidates, return_key:bool):
	'''
	Yield the (ra, dec, key) of every row of a coarse-step result.

	The coarse step is :py:func:`QLSCIndex.radial_query`, whose plain array result is
	``MatchArrayShape.SQUEEZED`` - a single candidate arrives with shape (2,) rather than
	(1,2) - so the plain form is put through :py:func:`_as_stable_pairs` first. The rows
	are yielded rather than collected: the caller keeps only the ones the region contains,
	and a whole-sky coarse step over a large index can have a great many it does not.

	:param candidates: the return value of the coarse step, asked for as plain values or with ``return_key``
	:param return_key: whether the coarse step was asked for the "key" field
	:returns: a generator of (ra, dec, key) tuples; the key is None when it was not asked for
	'''
	if return_key:
		# a recarray with the fields ra,dec,key; atleast_1d because an empty or
		# single-row recarray can arrive as a 0-dimensional array
		for row in np.atleast_1d(candidates):
			yield row["ra"], row["dec"], row["key"]
	else:
		for row in _as_stable_pairs(candidates):
			yield row[0], row[1], None

def _fast_path_matches(index:QLSCIndex, region:"ASTRegion", disc_radius:float,
					   return_key:bool, return_skycoord:bool):
	'''
	Answer a region query with one of qlsc's own exact region queries, when the region
	is a shape qlsc defines identically to AST.

	Two shapes qualify, and the claim in each case is that the two libraries describe
	the *same set of points*, not a close approximation of it:

	* an :class:`~cornish.ASTCircle` is the locus of points within a fixed geodesic
	  distance of a centre, which is exactly what :py:func:`QLSCIndex.radial_query`
	  searches;
	* an :class:`~cornish.ASTPolygon` joins its vertices with great circles, which is
	  exactly how :py:func:`QLSC.point_in_polygon` — and therefore
	  :py:func:`QLSCIndex.polygon_query` — reads a polygon. That the two agree is not
	  assumed: it was measured over 690 probe points across five pixels with zero
	  disagreements when this bridge was built, and ``test_cornish_region_query.py``
	  re-measures the agreement over seeded catalogues on every run.

	An :class:`~cornish.ASTBox` deliberately does **not** qualify. AST's Box has sides
	parallel to the axes of its frame, so on a sky frame its edges follow lines of
	constant declination — small circles, not great circles — and the great-circle
	polygon through its four corners is neither a subset nor a superset of it. (Verified
	directly: the midpoint of the bottom edge of a box spanning dec 20..30 reads inside
	at dec=20.00 and outside at dec=19.98, i.e. the edge is the parallel of declination.)
	A box therefore takes the generic path, where AST itself decides membership.

	Anything else — MOCs, compound regions, region types cornish gains later — also
	takes the generic path, which is correct for every region AST can describe.

	The built-in query's result is returned as it comes, not rebuilt: both of them
	already honour ``return_key`` and ``return_skycoord`` and assemble the result with
	the same helper the generic path uses. The one difference is
	:py:func:`radial_query`'s historical squeeze of its *plain* array, which
	:py:func:`_as_stable_pairs` undoes; its keyed, SkyCoord and QTable results are not
	squeezed and need nothing done to them.

	:param index: the :class:`~qlsc.QLSCIndex` to search
	:param region: the region, already checked to be a bounded, non-negated sky region
	:param disc_radius: the radius in degrees of the region's bounding circle
	:param return_key: whether the caller asked for the "key" field
	:param return_skycoord: whether the caller asked for the coordinates as a SkyCoord
	:returns: the query result, or :py:data:`_NO_FAST_PATH` when there is no fast path for this region
	'''
	if isinstance(region, ASTCircle):
		centre = region.centre
		radius = _as_float_scalar(region.radius.to("deg").value, "circle radius")
		matches = index.radial_query(float(centre[0]), float(centre[1]), radius,
									 return_key=return_key, return_skycoord=return_skycoord)
		# radial_query is the one query whose plain array is squeezed
		if return_key or return_skycoord:
			return matches
		return _as_stable_pairs(matches)

	if isinstance(region, ASTPolygon):
		vertices = np.asarray(region.points, dtype=np.double)

		# Q3C's polygon routines take between 3 and Q3C_MAX_N_POLY_VERTEX vertices, and
		# the limits are read from the extension module rather than repeated here. AST
		# does not build a polygon with fewer than three, so the lower bound is a
		# latent case; it is checked anyway so that the dispatch is total rather than
		# total-by-coincidence.
		if not (MINIMUM_POLYGON_VERTICES <= len(vertices) <= MAXIMUM_POLYGON_VERTICES):
			logger.debug(f"The polygon has {len(vertices)} vertices, outside the "
						 f"{MINIMUM_POLYGON_VERTICES}..{MAXIMUM_POLYGON_VERTICES} Q3C's polygon "
						 "query accepts; the region query will use its generic path instead.")
			return _NO_FAST_PATH

		# Q3C projects a polygon onto a single cube face, so a large one is refused
		# outright (ValueError). The size at which that happens depends on where the
		# polygon sits - measured between 43 and 76 degrees across for a square at five
		# positions (see docs/examples/survey_footprint.rst) - and the API documentation
		# quotes ~23 degrees as the figure that is safe anywhere. The fast path is taken
		# only well inside that: a bounding circle of this radius holds a polygon no more
		# than twice it across. A larger polygon is not refused, it just takes the generic
		# path, which has no such limit.
		if disc_radius > POLYGON_FAST_PATH_MAXIMUM_RADIUS:
			logger.debug(f"The polygon's bounding circle has a radius of {disc_radius} deg, beyond the "
						 f"{POLYGON_FAST_PATH_MAXIMUM_RADIUS} deg for which Q3C's polygon query is known to "
						 "project at any position on the sky; the region query will use its generic path instead.")
			return _NO_FAST_PATH

		# The winding order is not reversed on the way in, unlike ipix2region on the way
		# out: AST reads the interior as the left of the vertex path, but Q3C's polygon
		# test is winding-agnostic and answers about the same (small) side either way.
		# test_polygon_query_is_winding_agnostic pins that.
		try:
			return index.polygon_query(polygon=vertices, return_key=return_key,
									   return_skycoord=return_skycoord)
		except q3c.RangeOverflowError:
			# The coverage step wanted more ipix ranges than it has room for. A cone has a
			# size to grow on a retry and a polygon does not, so polygon_query lets this
			# through - but the generic path answers the same question with a cone, which
			# does. Nothing is skipped: the query is re-run by the other route, and the
			# answer is the same set of points.
			logger.warning("The polygon query for this region overflowed Q3C's ipix range capacity; "
						   "the region query is falling back to its generic bounding-circle path, which "
						   "returns the same matches.")
			return _NO_FAST_PATH

	return _NO_FAST_PATH

def region_query(index:QLSCIndex, region:"ASTRegion", return_key:bool=False,
				 return_skycoord:bool=False, all_sky:bool=False,
				 frame:Union[str,None]=None, skip_frame_check:bool=False,
				 fast_path:bool=True) -> np.ndarray:
	'''
	Return the points of a :class:`~qlsc.QLSCIndex` that fall inside a cornish region.

	This is the function behind :py:func:`QLSCIndex.region_query`, which is how it is
	normally reached; it is kept as a module function so that the bridge can also be
	used without going through the index's own namespace.

	Any bounded AST region works — :class:`~cornish.ASTCircle`, :class:`~cornish.ASTBox`,
	:class:`~cornish.ASTPolygon`, :class:`~cornish.ASTMoc`, and the compound regions built
	from them with :class:`~cornish.ASTCompoundRegion` — because membership is decided by
	AST itself, not re-derived here.

	The query has the same two steps as :py:func:`QLSCIndex.radial_query` and
	:py:func:`QLSCIndex.polygon_query`:

	1. a **coarse step** that must return a superset of the answer: a cone search over
	   the region's bounding circle (``region.boundingCircle()``), which is where the
	   index does its work;
	2. an **exact step**: every candidate is tested with ``region.pointInRegion``.

	Because the coarse step only has to be a superset, an over-generous bounding circle
	costs time and never correctness. The reverse — a bounding circle that does *not*
	contain the whole region — silently loses points, so two cases are refused outright:

	* a **negated region** (``region.isNegated``), whose bounding circle AST reports as
	  the circle of the *un-negated* region, i.e. the exact complement of what bounds it;
	* any region whose bounding circle demonstrably fails to contain it, detected by
	  testing the point antipodal to the circle's center.

	Both are raised as a :class:`ValueError` recommending ``all_sky=True``, which replaces
	the coarse step with a scan of the whole index (a cone of radius 180 degrees). That is
	correct for any region at all, at the cost of testing every indexed point.

	**Fast paths.** Two region types describe a shape qlsc itself defines identically, and
	for those the whole query is handed to the matching built-in method — the exact test
	then runs as one batched call into the C layer instead of one ``pointInRegion`` call
	per candidate:

	================================  ==================================================
	region                            query
	================================  ==================================================
	:class:`~cornish.ASTCircle`       :py:func:`QLSCIndex.radial_query` — the same geodesic disc
	:class:`~cornish.ASTPolygon`      :py:func:`QLSCIndex.polygon_query` — the same great-circle edges
	================================  ==================================================

	A polygon takes its fast path only when it has at most
	``qlsc.q3c.Q3C_MAX_N_POLY_VERTEX`` vertices and its bounding circle is no larger than
	``POLYGON_FAST_PATH_MAXIMUM_RADIUS``, the size for which Q3C's single-face projection
	is documented to work anywhere on the sky. Every other region — a
	:class:`~cornish.ASTBox`, whose sky edges are parallels of declination rather than
	great circles, a MOC, a compound region — takes the generic path, and so does every
	region when ``all_sky`` is set. The matches are the same either way, and
	``fast_path=False`` forces the generic path for a caller who wants AST to be the sole
	arbiter of membership (the test suite uses it to check the two agree).

	One asymmetry is deliberate and worth knowing about. Q3C decides for itself whether a
	polygon will project onto a single cube face, and raises ``ValueError`` when it will
	not; the size guard above is set far enough inside the documented limit that this is
	not expected to happen, but if it ever does the error is **raised rather than
	answered by the generic path**, because a guard that has been proved wrong should say
	so rather than route around itself. Q3C's other refusal, a
	``qlsc.q3c.RangeOverflowError`` from the coverage step, *is* answered by the generic
	path — a cone has a size to grow and retry with where a polygon has none — with a
	warning logged and the same matches returned. (Telling the two apart without matching
	on the message text needs a dedicated exception class in the extension module, the
	way ``RangeOverflowError`` already has one; until then, the projection failure is
	simply not caught.)

	The result is assembled by the same helper the built-in region queries use, and like
	every query added since :py:func:`radial_query` it is ``MatchArrayShape.STABLE``: the
	plain array is always ``(n,2)``, including for a single match. ``radial_query``'s
	historical squeeze to ``(2,)`` is not inherited, even when it supplies the answer.

	:param index: the :class:`~qlsc.QLSCIndex` to search
	:param region: a cornish :class:`~cornish.ASTRegion` defined in a sky frame
	:param return_key: if set to True, include the key value as provided when added to the index
	:param return_skycoord: if set to True, return the coordinates as an :class:`~astropy.coordinates.SkyCoord`
	:param all_sky: if set to True, consider every point in the index instead of using the region's bounding circle
	:param frame: name of the frame the index coordinates are in, compared against the region's frame; None means the default (ICRS)
	:param skip_frame_check: if set to True, waive that comparison and take responsibility for the frames matching
	:param fast_path: if set to False, always use the generic bounding-circle path, whatever the region's type
	:returns: an array of matches with shape (n,2), a :class:`numpy.recarray` with the fields ``ra,dec,key``, a :class:`~astropy.coordinates.SkyCoord`, or an :class:`~astropy.table.QTable`
	'''
	_require_cornish()

	if not isinstance(index, QLSCIndex):
		raise TypeError(f"A qlsc.QLSCIndex object was expected for 'index'; was given '{type(index).__name__}'.")

	_check_region(region, frame, skip_frame_check)

	if all_sky:
		# a cone of radius 180 degrees centered anywhere covers the sphere
		center_ra, center_dec, radius = 0.0, 90.0, 180.0
	else:
		if region.isNegated:
			raise ValueError(
				"A negated region cannot be searched with a bounding-circle coarse step: AST "
				"returns the bounding circle of the un-negated region, which contains none of "
				"the negated region's area, so the search would silently return nothing. "
				"Pass all_sky=True to test every point in the index instead.")

		# ValueError from cornish when the region is empty or unbounded; NotA2DRegion when not 2D
		disc = region.boundingCircle()
		center_ra, center_dec = (float(v) for v in disc.centre)
		radius = _as_float_scalar(disc.radius.to("deg").value, "bounding circle radius")

		# Falsify the bound where it is cheap to do so: if the region contains the point
		# directly opposite the circle's center, the circle is not a bound. This catches
		# compound regions with a negated component, whose top-level isNegated is False,
		# and a polygon wound the other way round, which AST reads as the rest of the sky.
		if radius < 180.0 and region.pointInRegion([(center_ra + 180.0) % 360.0, -center_dec]):
			raise ValueError(
				f"The bounding circle AST reports for this region (center=({center_ra}, {center_dec}), "
				f"radius={radius} deg) does not contain the region: the region contains the point "
				"antipodal to that center. Using it as the coarse step would silently lose matches. "
				"Pass all_sky=True to test every point in the index instead.")

		if fast_path:
			# The bound has been checked above, so a fast path answers about the same region
			# the generic path would have.
			fast = _fast_path_matches(index, region, radius, return_key, return_skycoord)
			if fast is not _NO_FAST_PATH:
				return fast

	# the coarse step is always asked for as plain values: its rows are candidates to be
	# filtered here, not the result, and return_skycoord is applied to the matches below
	candidates = index.radial_query(center_ra, center_dec, radius, return_key=return_key)

	match_ra = list()
	match_dec = list()
	match_key = list()

	for row_ra, row_dec, row_key in _candidate_rows(candidates, return_key):
		if region.pointInRegion([float(row_ra), float(row_dec)]):
			match_ra.append(row_ra)
			match_dec.append(row_dec)
			if return_key:
				match_key.append(row_key)

	return _match_result(match_ra, match_dec, match_key, None, return_key, False,
						 MatchArrayShape.STABLE, return_skycoord=return_skycoord)

def image_query(index:QLSCIndex, fits_header=None, return_key:bool=False,
				return_skycoord:bool=False, all_sky:bool=False, frame:Union[str,None]=None,
				skip_frame_check:bool=False, maxerr=None, maxvert:int=None,
				fast_path:bool=True) -> np.ndarray:
	'''
	Return the points of a :class:`~qlsc.QLSCIndex` that fall on the field described by a
	FITS header — "which of my sources land on this image?".

	This is the function behind :py:func:`QLSCIndex.image_query`, which is how it is
	normally reached.

	The header's WCS is read by AST, which returns the sky region bounding the image's
	field (an :class:`~cornish.ASTPolygon` following the field outline, carrying the
	pixel-to-world frame set in its ``.wcs``); that region is then handed to
	:py:func:`region_query`, so everything that function documents — the coarse ipix step,
	the exact per-candidate test, the fast paths, the frame check — applies unchanged.
	qlsc does no WCS arithmetic of its own here.

	The region comes back in whatever frame the header declares, so a header with
	``RADESYS = 'FK5'`` searched against an index of ICRS coordinates raises rather than
	mixing the two; pass ``frame`` naming the index's own frame, or ``skip_frame_check``,
	exactly as :py:func:`region_query` describes. ``all_sky`` is forwarded for the same
	reason: :py:func:`region_query`'s two refusals both recommend it, and a caller who
	reaches one of them from here must be able to take the advice.

	.. warning::

		**A very wide field can crash the interpreter, and the fault is upstream.**
		``ASTRegion.fromFITSHeader`` segmentation-faults (SIGSEGV, exit 139) on a header
		describing a field 180 degrees across — reproduced on AST 9.3.0 /
		starlink-pyast 4.0.1 with cornish alone and qlsc not imported, using a 180x180
		``RA---CAR``/``DEC--CAR`` header at 1 degree per pixel. A 90-degree field is
		fine, and 180 degrees in right ascension by 90 in declination raises ``Ast.NODEF``
		rather than crashing. Nothing qlsc can do catches a segmentation fault; the note
		is here so the crash is not read as a qlsc defect. Ordinary images are nowhere
		near this.

	:param index: the :class:`~qlsc.QLSCIndex` to search
	:param fits_header: a FITS header — an astropy :class:`~astropy.io.fits.Header`, a fitsio header, an array of cards, a dict, or a string
	:param return_key: if set to True, include the key value as provided when added to the index
	:param return_skycoord: if set to True, return the coordinates as an :class:`~astropy.coordinates.SkyCoord`
	:param all_sky: if set to True, consider every point in the index instead of using the field's bounding circle
	:param frame: name of the frame the index coordinates are in, compared against the header's frame; None means the default (ICRS)
	:param skip_frame_check: if set to True, waive that comparison and take responsibility for the frames matching
	:param maxerr: maximum deviation of the field polygon from the true field outline, as an astropy :class:`~astropy.units.Quantity`; None uses cornish's default
	:param maxvert: maximum number of vertices in the field polygon; None uses cornish's default
	:param fast_path: if set to False, always use the generic bounding-circle path
	:returns: an array of matches with shape (n,2), a :class:`numpy.recarray` with the fields ``ra,dec,key``, a :class:`~astropy.coordinates.SkyCoord`, or an :class:`~astropy.table.QTable`
	'''
	_require_cornish()

	if fits_header is None:
		raise ValueError("A value for 'fits_header' must be specified.")

	# Passed only when given: cornish owns the defaults for both, and repeating them here
	# would be a second copy to fall out of step with the first.
	region_parameters = dict()
	if maxerr is not None:
		region_parameters["maxerr"] = maxerr
	if maxvert is not None:
		region_parameters["maxvert"] = maxvert

	region = ASTRegion.fromFITSHeader(fits_header, **region_parameters)

	return region_query(index, region, return_key=return_key, return_skycoord=return_skycoord,
						all_sky=all_sky, frame=frame, skip_frame_check=skip_frame_check,
						fast_path=fast_path)
