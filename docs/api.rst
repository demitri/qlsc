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


Point Generation
----------------

The method below is useful to generate a number of
evenly distributed points at deterministic locations
for testing purposes.

.. autofunction:: qlsc.generate.sunflower_points_on_sphere


Exceptions
----------

Beyond the standard ``ValueError``, ``TypeError`` and ``KeyError`` raised for
invalid arguments, one exception comes from the embedded Q3C code:

``qlsc.q3c.RangeOverflowError``
	Raised when a region query asks the coverage step for more ipix ranges than
	it has room for. :py:func:`QLSCIndex.radial_query` and
	:py:func:`QLSCIndex.ellipse_query` recover from this on their own by widening
	the search region and retrying, so it does not normally reach the caller;
	:py:func:`QLSCIndex.polygon_query` has no size to widen and lets it through.
	See :doc:`examples/survey_footprint`.

A closed :py:class:`QLSCIndex` raises ``RuntimeError`` on any further use.
