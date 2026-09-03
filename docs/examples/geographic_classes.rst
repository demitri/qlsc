Geographic Longitude and Latitude
====================================

.. contents::
   :local:
   :depth: 2

**The problem.** The quadrilateralized spherical cube was designed for the
Earth, and qlsc's astronomical vocabulary — right ascension in [0,360),
declination — is one reading of it. Sometimes the positions you actually have
are geographic: longitude in the usual [-180,180) convention, latitude. You
would rather not translate back and forth by hand every time you call the
library.

:py:class:`~qlsc.GeoQLSC` and :py:class:`~qlsc.GeoQLSCIndex` are the
geographic reading of :py:class:`~qlsc.QLSC` and :py:class:`~qlsc.QLSCIndex`.
**Nothing about the pixelation changes** — same ipix numbers, same geometry —
only the vocabulary, and, for latitude, the strictness.


Wrapping longitude
-------------------

:py:func:`~qlsc.wrap_longitude` is the standalone conversion, if all you need
is the mapping itself:

.. code-block:: python

	from qlsc import GeoQLSC, GeoQLSCIndex, wrap_longitude

	print(wrap_longitude(286.0))
	print(wrap_longitude(-74.0))
	print(wrap_longitude(180.0))
	print(wrap_longitude(540.0))

.. code-block:: text

	-74.0
	-74.0
	-180.0
	-180.0

The mapping is idempotent and periodic in 360 degrees: ``286`` and ``-74`` are
the same meridian, and exactly ``180`` comes back as ``-180`` — the same
meridian again, on the boundary of the half-open range, matching astropy's
``Angle.wrap_at(180 * u.deg)`` convention.


Indexing a few cities
-----------------------

.. code-block:: python

	scheme = GeoQLSC(depth=20)

	cities = {
		"New York": (-74.0060, 40.7128),
		"Paris":    (2.3522, 48.8566),
		"Tokyo":    (139.6917, 35.6895),
		"Sydney":   (151.2093, -33.8688),
	}

	index = GeoQLSCIndex(qlsc=scheme)
	for name, (lon, lat) in cities.items():
		index.add_point(lon=lon, lat=lat, key=name)

	print("cities indexed:", len(index))

.. code-block:: text

	cities indexed: 4

``GeoQLSCIndex`` needs a ``GeoQLSC`` behind it (not a plain ``QLSC``) so that
its geometry methods and its stored longitudes agree on the convention; the
constructor refuses the mismatch rather than reporting one convention while
storing another.


Either spelling, anywhere
----------------------------

Every method that takes a position accepts it as ``lon``/``lat`` or as
``ra``/``dec`` — never both — and both spellings reach the same pixel:

.. code-block:: python

	nyc_by_lon = scheme.lonlat2ipix(lon=-74.0060, lat=40.7128)
	nyc_by_ra  = scheme.ang2ipix(ra=285.9940, dec=40.7128)
	print("same pixel either spelling:", nyc_by_lon == nyc_by_ra)

	lon, lat = scheme.ipix2lonlat(nyc_by_lon)
	print(f"round-trip: lon={lon:.4f} lat={lat:.4f}")

.. code-block:: text

	same pixel either spelling: True
	round-trip: lon=-74.0060 lat=40.7128

Longitude is periodic on the way in — ``-74.0060`` and its ``285.9940``
equivalent are accepted anywhere a right ascension can be — but every method
of these classes reports longitude in [-180,180) on the way out, never
[0,360).

Querying works the same way. Note that the returned record's fields are still
named ``ra``/``dec``, not ``lon``/``lat`` — the geographic classes only add
the lon/lat spelling as keyword arguments; the values that come back keep
qlsc's field names, holding whichever convention this class reports:

.. code-block:: python

	matches = index.radial_query(lon=-74.0, lat=40.7, radius=1.0, return_key=True)
	print("field names are still ra,dec:", matches.dtype.names)
	for m in matches:
		print(f"  {m['key']:10s}  lon={m['ra']:9.4f}  lat={m['dec']:+8.4f}")

.. code-block:: text

	field names are still ra,dec: ('ra', 'dec', 'key')
	  New York    lon= -74.0060  lat=+40.7128

.. code-block:: python

	paris = index.radial_query(ra=2.3522, dec=48.8566, radius=0.01, return_key=True)
	print("found by ra/dec spelling too:", paris["key"][0])

.. code-block:: text

	found by ra/dec spelling too: Paris


Latitude is checked, not reinterpreted
-------------------------------------------

:py:class:`~qlsc.QLSC` reads a declination outside [-90,90] as an unnormalized
celestial coordinate — a point past the pole, folded to the far side and 180
degrees around in right ascension. A latitude has no such reading: there is no
place at latitude 100, so the geographic classes raise rather than silently
computing an answer for a point you did not mean:

.. code-block:: python

	try:
		GeoQLSC().point_in_polygon(lat=100.0, lon=0.0, polygon=[[0., -1.], [1., -1.], [1., 1.], [0., 1.]])
	except ValueError as e:
		print(type(e).__name__, e)

	try:
		index.add_point(lon=0., lat=-95.)
	except ValueError as e:
		print(type(e).__name__, e)

.. code-block:: text

	ValueError The 'lat' value(s) must be in the range [-90,90]; a latitude outside it is not a point on the far side of a pole, it is an error.
	ValueError The 'lat' value(s) must be in the range [-90,90]; a latitude outside it is not a point on the far side of a pole, it is an error.

The contrast with plain ``QLSC`` is exact — same out-of-range value, opposite
outcome:

.. code-block:: python

	from qlsc import QLSC

	celestial = QLSC(depth=6)
	p = celestial.ang2ipix(ra=100., dec=-95.)
	folded_ra, folded_dec = celestial.ipix2ang(p)
	print(f"QLSC folds dec=-95 to the far side of the pole: ra={folded_ra:.2f} dec={folded_dec:.2f}")
	print("GeoQLSC refuses the same value outright, as shown above.")

.. code-block:: text

	QLSC folds dec=-95 to the far side of the pole: ra=281.31 dec=-85.44
	GeoQLSC refuses the same value outright, as shown above.

This is the one way the geographic classes are stricter than the celestial
ones, and it holds at every entry point that takes a latitude, under either
spelling.


How longitude is stored
--------------------------

A point is stored with its longitude already wrapped into [-180,180), so two
insertions of the "same" meridian under different spellings collapse to one
row, the way any duplicate coordinate does (see
:doc:`../user_guide/indexing`'s discussion of coordinate identity):

.. code-block:: python

	wide = GeoQLSCIndex(qlsc=GeoQLSC(depth=20))
	wide.add_point(lon=286., lat=10., key="a")   # 286 == -74, unwrapped on the way in
	wide.add_point(lon=-74., lat=10., key="b")   # same point, already wrapped
	print("same coordinate is one row regardless of which spelling was used:", len(wide))

	stored = wide.radial_query(lon=-74., lat=10., radius=0.001, return_key=True)
	print("stored longitude is wrapped to [-180,180):", stored["ra"][0])
	print("only the first key survives (see the indexing user guide):", stored["key"][0])

.. code-block:: text

	same coordinate is one row regardless of which spelling was used: 1
	stored longitude is wrapped to [-180,180): -74.0
	only the first key survives (see the indexing user guide): a

An index file is interchangeable between :py:class:`~qlsc.GeoQLSCIndex` and
plain :py:class:`~qlsc.QLSCIndex`: the numbers on disk are the same either
way, and only the range they are read back in differs. A file written by
:py:class:`~qlsc.QLSCIndex` may hold longitudes in [0,360) rather than
[-180,180); those come back as stored.


The one exception: exported coverage
------------------------------------------

:py:func:`~qlsc.QLSCIndex.coverage` is inherited unchanged, and it declares
the ICRS frame — which has its own convention, [0,360), regardless of which
class produced the document:

.. code-block:: python

	cov = index.coverage(depth=2)
	ra_values = [v[0] for p in cov["polygons"] for v in p["vertices"]]
	print("min:", round(min(ra_values), 3), "max:", round(max(ra_values), 3))
	print("frame:", cov["frame"], "(always [0,360) regardless of GeoQLSCIndex's own convention)")

.. code-block:: text

	min: 0.0 max: 296.565
	frame: icrs (always [0,360) regardless of GeoQLSCIndex's own convention)

A sky-coverage document is an interchange artifact for other software, so its
numbers follow the frame's convention rather than this object's usual
reporting preference — the one place a :py:class:`~qlsc.GeoQLSCIndex` reports
right ascension rather than longitude.


See also
--------

* :py:class:`~qlsc.GeoQLSC` and :py:class:`~qlsc.GeoQLSCIndex` — the full
  class documentation, including every method's lon/lat signature.
* :py:func:`~qlsc.wrap_longitude` — the standalone conversion.
* :doc:`coverage_export` — what :py:func:`~qlsc.QLSCIndex.coverage` produces,
  including the [0,360) exception noted above.
* :doc:`../user_guide/indexing` — coordinate identity, keys, and the query
  methods these classes inherit unchanged in every way but vocabulary.
