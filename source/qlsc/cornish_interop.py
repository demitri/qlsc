
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

from typing import Iterable, Union

import numpy as np

try:
	import cornish
	from cornish import ASTRegion, ASTPolygon, ASTMoc, ASTICRSFrame
	cornish_available = True
	cornish_import_error = None
except ImportError as e:
	cornish_available = False
	# Keep the original error: "cornish is not installed" and "cornish is installed but its
	# own starlink-pyast import failed" are different problems with the same symptom here,
	# and reporting the first when it is really the second sends the user the wrong way.
	cornish_import_error = e

from .qlsc_base import QLSC, QLSCIndex, MatchArrayShape, _match_result
from .utilities import _as_float_scalar

# The frame the index coordinates are assumed to be in when the caller does not say.
DEFAULT_FRAME = "ICRS"

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

	moc = ASTMoc(max_order=max_order)
	centers = list()
	for value in ipix:
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

def region_query(index:QLSCIndex, region:"ASTRegion", return_key:bool=False,
				 return_skycoord:bool=False, all_sky:bool=False,
				 frame:Union[str,None]=None, skip_frame_check:bool=False) -> np.ndarray:
	'''
	Return the points of a :class:`~qlsc.QLSCIndex` that fall inside a cornish region.

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

	The result is assembled by the same helper the built-in region queries use, and like
	every query added since :py:func:`radial_query` it is ``MatchArrayShape.STABLE``: the
	plain array is always ``(n,2)``, including for a single match. ``radial_query``'s
	historical squeeze to ``(2,)`` is not inherited, even though it supplies the candidates.

	:param index: the :class:`~qlsc.QLSCIndex` to search
	:param region: a cornish :class:`~cornish.ASTRegion` defined in a sky frame
	:param return_key: if set to True, include the key value as provided when added to the index
	:param return_skycoord: if set to True, return the coordinates as an :class:`~astropy.coordinates.SkyCoord`
	:param all_sky: if set to True, consider every point in the index instead of using the region's bounding circle
	:param frame: name of the frame the index coordinates are in, compared against the region's frame; None means the default (ICRS)
	:param skip_frame_check: if set to True, waive that comparison and take responsibility for the frames matching
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
		# compound regions with a negated component, whose top-level isNegated is False.
		if radius < 180.0 and region.pointInRegion([(center_ra + 180.0) % 360.0, -center_dec]):
			raise ValueError(
				f"The bounding circle AST reports for this region (center=({center_ra}, {center_dec}), "
				f"radius={radius} deg) does not contain the region: the region contains the point "
				"antipodal to that center. Using it as the coarse step would silently lose matches. "
				"Pass all_sky=True to test every point in the index instead.")

	# the coarse step is always asked for as plain values: its rows are candidates to be
	# filtered here, not the result, and return_skycoord is applied to the matches below
	candidates = index.radial_query(center_ra, center_dec, radius, return_key=return_key)

	match_ra = list()
	match_dec = list()
	match_key = list()

	if return_key:
		candidates = np.atleast_1d(candidates)
		rows = ((row["ra"], row["dec"], row["key"]) for row in candidates)
	else:
		# radial_query squeezes a single match down to shape (2,); restore the (n,2) form
		candidates = np.atleast_2d(candidates)
		rows = ((row[0], row[1], None) for row in candidates)

	for row_ra, row_dec, row_key in rows:
		if region.pointInRegion([float(row_ra), float(row_dec)]):
			match_ra.append(row_ra)
			match_dec.append(row_dec)
			if return_key:
				match_key.append(row_key)

	return _match_result(match_ra, match_dec, match_key, None, return_key, False,
						 MatchArrayShape.STABLE, return_skycoord=return_skycoord)
