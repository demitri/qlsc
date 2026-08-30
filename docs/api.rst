API
===

.. module:: qlsc

The interfaces to QLSC are described here.

QLSC Class
----------

.. autoclass:: QLSC
   :members:
   :inherited-members:
   
QLSC Index
----------

.. autoclass:: QLSCIndex
   :members:
   :inherited-members:

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

Cross-matching
--------------

.. autofunction:: qlsc.crossmatch

Point Generation
----------------

The method below is useful to generate a number of
evenly distributed points at deterministic locations
for testing purposes.

.. autofunction:: qlsc.generate.sunflower_points_on_sphere
