API
===

.. module:: qlsc

The interfaces to QLSC are described here. The :doc:`user guide <user_guide/index>`
explains what they are for; this page is the reference.

.. contents::
   :local:
   :depth: 2


QLSC Class
----------

The pixellation scheme: conversions between sky coordinates, pixel numbers, and
positions on a cube face. See :doc:`user_guide/coordinates`.

Coordinates may be given as ``ra,dec`` in degrees, as
:class:`~astropy.units.Quantity` angles, or as a scalar
:class:`~astropy.coordinates.SkyCoord` passed as the single ``ra`` argument, which
is transformed to ICRS first.

.. autoclass:: QLSC
   :members:
   :inherited-members:


QLSC Index
----------

A searchable collection of sky positions. See :doc:`user_guide/indexing`.

.. autoclass:: QLSCIndex
   :members:
   :inherited-members:


Module Functions
----------------

.. autofunction:: qlsc.distance

.. autofunction:: qlsc.sindist

.. autofunction:: qlsc.xy2facenum


Cross-matching
--------------

.. autofunction:: qlsc.crossmatch


Geographic (longitude/latitude) Classes
---------------------------------------

The same scheme and index, in geographic vocabulary: latitude is declination, and
longitude is right ascension reported in [-180,180) rather than [0,360).

.. autoclass:: GeoQLSC
   :members:
   :inherited-members:

.. autoclass:: GeoQLSCIndex
   :members:
   :inherited-members:

.. autofunction:: qlsc.wrap_longitude


Starlink AST regions, via cornish
---------------------------------

.. module:: qlsc.cornish_interop

Bridges a QLSC scheme and index to `cornish <https://github.com/demitri/cornish>`_,
a Python interface to the Starlink AST library: a pixel becomes an ``ASTPolygon``,
a set of pixels becomes a MOC, and any ``ASTRegion`` becomes a query over an index.

.. note::

	This module needs `cornish <https://github.com/demitri/cornish>`_, which is an
	optional dependency and is not installed with qlsc. The module imports without
	it; calling any of the functions below raises ``ImportError`` naming the package
	and the reason the import failed.

.. autofunction:: qlsc.cornish_interop.ipix2region

.. autofunction:: qlsc.cornish_interop.ipix2moc

.. autofunction:: qlsc.cornish_interop.region_query


Point Generation
----------------

The method below is useful to generate a number of
evenly distributed points at deterministic locations
for testing purposes.

.. autofunction:: qlsc.generate.sunflower_points_on_sphere


Exceptions and constants
------------------------

Beyond the standard ``ValueError``, ``TypeError`` and ``KeyError`` raised for
invalid arguments, the compiled extension module ``qlsc.q3c`` provides:

``qlsc.q3c.RangeOverflowError``
	Raised when a region query asks the coverage step for more ipix ranges than
	it has room for. :py:func:`QLSCIndex.radial_query` and
	:py:func:`QLSCIndex.ellipse_query` recover from this on their own by widening
	the search region and retrying, so it does not normally reach the caller;
	:py:func:`QLSCIndex.polygon_query` has no size to widen and lets it through.
	See :doc:`examples/survey_footprint`.

``qlsc.q3c.ELLIPSE_MAXIMUM_SEMI_MAJOR_AXIS``
	The largest semi-major axis, in degrees, for which the elliptical region
	geometry is defined (90). :py:func:`QLSC.point_in_ellipse` and
	:py:func:`QLSCIndex.ellipse_query` raise ``ValueError`` above it rather than
	answering about a region that is not the one asked for. An ``axis_ratio`` of 1
	is a circle and is not subject to the limit.

A closed :py:class:`QLSCIndex` raises ``RuntimeError`` on any further use.
