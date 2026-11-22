Indexing and Searching
======================

.. contents::
   :local:
   :depth: 2

:py:class:`~qlsc.QLSCIndex` holds a collection of sky positions and answers
spatial questions about them: what is inside this cone, this polygon, this
ellipse. It is the part of the package that replaces "install PostgreSQL and the
Q3C extension".


Creating an index
-----------------

With no arguments you get an index in memory at ``depth=30``:

.. code-block:: python

	from qlsc import QLSC, QLSCIndex

	index = QLSCIndex()

Give it a file path and it persists. If the file exists and is a valid index it
is opened; otherwise it is created:

.. code-block:: python

	index = QLSCIndex(filepath="my_catalog.qlscidx",
	                  name="demo", description="a small demo index")

The path may be a string, a :class:`pathlib.Path`, or a ``file://`` URL. The
extension is up to you; ``.qlscidx`` is a useful convention. Pass a
:py:class:`~qlsc.QLSC` to use a different resolution:

.. code-block:: python

	index = QLSCIndex(qlsc=QLSC(depth=12), filepath="coarse.qlscidx")

An index supports the context manager protocol, and can be closed explicitly:

.. code-block:: python

	import os
	import numpy as np
	from qlsc import QLSC, QLSCIndex

	path = "demo.qlscidx"

	with QLSCIndex(filepath=path, name="demo", description="a small demo index") as index:
		index.add_points(ra=np.array([1., 2., 3.]), dec=np.array([1., 2., 3.]))
		print("inside the context manager:", len(index))

	print("file size on disk:", os.path.getsize(path), "bytes")

	reopened = QLSCIndex(filepath=path)
	print("reopened:", len(reopened), reopened.name, "/", reopened.description)
	reopened.close()
	try:
		len(reopened)
	except Exception as e:
		print("using a closed index:", type(e).__name__, e)

.. code-block:: text

	inside the context manager: 3
	file size on disk: 24576 bytes
	reopened: 3 demo / a small demo index
	using a closed index: RuntimeError This QLSCIndex has been closed and can no longer be used.

An in-memory index vanishes when the object goes out of scope; a file-backed one
does not, and closing it is the tidy thing to do.

.. note::

	The storage is an SQLite database, and connections are opened and closed
	around each operation so that several threads can use a file-backed index at
	once. The in-memory case is the exception: its one connection must stay open
	for the lifetime of the object, because the database *is* the connection.
	The schema is deliberately not part of the public interface and may change.


Adding points
-------------

:py:func:`~qlsc.QLSCIndex.add_point` takes one position;
:py:func:`~qlsc.QLSCIndex.add_points` takes many and is much faster than a loop.
Coordinates are degrees. Points may be given as separate ``ra`` and ``dec``
iterables or as an ``(n,2)`` ``points`` array, and either form may carry ``keys``:

.. code-block:: python

	index.add_point(56.59625, 78.00866, key="NGC 1234")
	index.add_points(ra=ra_array, dec=dec_array)
	index.add_points(points=points_array, keys=["a", "b", "c"])

A **coordinate is the identity of a point**. Adding the same coordinate twice
does not create a second row:

.. code-block:: python

	index = QLSCIndex()
	index.add_point(10., 20.)
	index.add_point(10., 20.)
	index.add_point(10., 20.0000001)
	print("points:", len(index))

.. code-block:: text

	points: 2

That is exact equality, not a match radius — the third point differs in the
seventh decimal place and is a separate row. It is worth knowing when loading a
real catalog: the Yale Bright Star Catalogue, for instance, contains fourteen
pairs of components listed at identical rounded positions, and an index built
from it holds fourteen fewer rows than the file has lines (see
:doc:`../examples/bright_star_mask`).

**Keys are labels, not identities.** They need not be unique, several points may
share one, and a point may have none:

.. code-block:: python

	index = QLSCIndex()
	index.add_points(ra=[10., 10., 11.], dec=[20., 20., 21.], keys=["a", "b", "c"])
	print("points:", len(index))
	print("keys:", sorted(index.keys()))

	index2 = QLSCIndex()
	index2.add_points(ra=[10., 11.], dec=[20., 21.], keys=["same", "same"])
	print("two points, one key:", len(index2), sorted(index2.keys()), "same" in index2)

.. code-block:: text

	points: 2
	keys: ['a', 'c']
	two points, one key: 2 ['same', 'same'] True

Note the first result: ``(10,20)`` was added twice with different keys, and the
second insert was discarded along with its key. The same happens if a coordinate
is added first without a key and then with one:

.. code-block:: python

	index = QLSCIndex()
	index.add_point(10., 20.)
	index.add_point(10., 20., key="named")
	print("points:", len(index), "keys:", sorted(index.keys()))

.. code-block:: text

	points: 1 keys: []

If your catalog has identifiers you care about, attach them on the first insert.


The dictionary-like interface
-----------------------------

An index answers ``len()``, ``in``, iteration, and
:py:func:`~qlsc.QLSCIndex.keys`, so a lot of code reads naturally:

.. code-block:: python

	index = QLSCIndex()
	index.add_points(ra=[1., 2., 3.], dec=[1., 2., 3.])
	index.add_point(4., 4., key="k")
	print("len:", len(index), "keys:", sorted(index.keys()), "bool:", bool(index))

	empty = QLSCIndex()
	print("empty index is falsy:", bool(empty), len(empty))

.. code-block:: text

	len: 4 keys: ['k'] bool: True
	empty index is falsy: False 0

Where it differs from a dictionary, and it matters:

* ``len(index)`` counts **points**, while ``index.keys()`` yields **keys**.
  Unkeyed points are not represented in ``keys()``, and a key shared by several
  points is yielded once per point, so the two counts can differ in either
  direction. ``index.number_of_points`` is the same value as ``len(index)``.
* ``keys()`` returns a generator, not a list, and it holds a database connection
  open until it is exhausted or closed.
* ``key in index`` asks whether *any* point carries that key. Keys must be
  strings — the underlying column has TEXT affinity, so a non-string would be
  compared as its text representation, and the class refuses rather than let that
  happen silently.

.. code-block:: python

	index = QLSCIndex()
	index.add_point(4., 4., key="k")
	try:
		42 in index
	except TypeError as e:
		print(type(e).__name__, e)

.. code-block:: text

	TypeError Index keys are strings (the underlying column has TEXT affinity); the key is '42' of type 'int'.

:py:func:`~qlsc.QLSCIndex.remove_point` deletes by key or by exact coordinate,
returns how many rows went, and raises rather than quietly doing nothing:

.. code-block:: python

	index = QLSCIndex()
	index.add_points(ra=[1., 2., 3.], dec=[1., 2., 3.], keys=["x", "y", "y"])
	print("removed by key 'y':", index.remove_point(key="y"), "-> len", len(index))
	print("removed by coordinate:", index.remove_point(ra=1., dec=1.), "-> len", len(index))
	try:
		index.remove_point(key="gone")
	except KeyError as e:
		print("removing something absent:", type(e).__name__, e)

.. code-block:: text

	removed by key 'y': 2 -> len 1
	removed by coordinate: 1 -> len 0
	removing something absent: KeyError "No point matching the key 'gone' was found in this index."

Removing by coordinate needs the coordinate that was stored — there is no
tolerance. Run a small :py:func:`~qlsc.QLSCIndex.radial_query` first and remove
the position it reports.


Querying
--------

Three region queries share a common shape:

+--------------------------------------------------+-----------------------------------------------+
| Method                                           | Region                                        |
+==================================================+===============================================+
| :py:func:`~qlsc.QLSCIndex.radial_query`          | a cone: centre and radius                     |
+--------------------------------------------------+-----------------------------------------------+
| :py:func:`~qlsc.QLSCIndex.polygon_query`         | a spherical polygon, convex or concave        |
+--------------------------------------------------+-----------------------------------------------+
| :py:func:`~qlsc.QLSCIndex.ellipse_query`         | an ellipse: centre, semi-major axis,          |
|                                                  | axis ratio, position angle                    |
+--------------------------------------------------+-----------------------------------------------+

All three take degrees, or :class:`~astropy.units.Quantity` angles if astropy is
installed, and all three accept ``return_key``. The cone and the ellipse also
accept ``return_distance``, the angular separation from the centre in degrees;
a polygon has no centre to measure from, so it has no such option.

**The return type depends on those flags.** With both left ``False`` you get a
plain ``(n,2)`` array of coordinates. With either set you get a
:class:`numpy.recarray` with named fields, in the order ``ra``, ``dec``,
``key``, ``distance``. Using a fresh index of fifty thousand points:

.. code-block:: python

	import numpy as np
	from qlsc import QLSCIndex

	rng = np.random.default_rng(seed=20200422)
	ra  = rng.uniform(0., 360., 50_000)
	dec = np.degrees(np.arcsin(rng.uniform(-1., 1., 50_000)))

	index = QLSCIndex()
	index.add_points(ra=ra, dec=dec, keys=[f"S{i:05d}" for i in range(50_000)])

	r1 = index.radial_query(ra=100., dec=10., radius=2.)
	r2 = index.radial_query(ra=100., dec=10., radius=2., return_key=True)
	r3 = index.radial_query(ra=100., dec=10., radius=2., return_distance=True)
	r4 = index.radial_query(ra=100., dec=10., radius=2., return_key=True, return_distance=True)

	for name, r in [("neither", r1), ("return_key", r2),
	                ("return_distance", r3), ("both", r4)]:
		print(f"{name:16s} type={type(r).__name__:10s} shape={r.shape} fields={r.dtype.names}")
	print(r4[0])

.. code-block:: text

	neither          type=ndarray    shape=(18, 2) fields=None
	return_key       type=recarray   shape=(18,) fields=('ra', 'dec', 'key')
	return_distance  type=recarray   shape=(18,) fields=('ra', 'dec', 'distance')
	both             type=recarray   shape=(18,) fields=('ra', 'dec', 'key', 'distance')
	(99.23592352159642, 8.642793354892902, 'S29764', 1.5525712338207096)

The record form is the more convenient one: it sorts by field
(``matches.sort(order="distance")``) and indexes by name
(``matches["key"]``).

.. warning::

	:py:func:`~qlsc.QLSCIndex.radial_query`'s plain ``(n,2)`` form is squeezed, so
	a **single match comes back with shape (2,)**, not ``(1,2)``:

	.. code-block:: python

		one = QLSCIndex()
		one.add_point(5., 5.)
		print(one.radial_query(ra=5., dec=5., radius=1.).shape)

	.. code-block:: text

		(2,)

	Code that iterates over rows will silently iterate over the two coordinates
	instead. Ask for ``return_distance=True`` — the record form always has shape
	``(n,)`` — or use ``np.atleast_2d`` on the result.

An empty result is an empty array, not ``None``, so ``len(result)`` is always
meaningful.


.. _how-a-region-query-works:

Querying many positions at once
-------------------------------

:py:func:`~qlsc.QLSCIndex.radial_query` also takes **arrays** of centres — a list
of targets to match against the index, which is what a catalog
cross-identification needs. The result is then always a
:class:`numpy.recarray` whose first field is ``target_index``, the position in the
target list that each match belongs to; the matches are one flat table, grouped by
target in the order the targets were given.

A target that matches nothing contributes no rows, so ``return_counts`` exists to
report the zeroes that the table cannot:

.. code-block:: python

	targets_ra  = np.array([100., 200., 300.])
	targets_dec = np.array([10., -20., 30.])

	matches, counts = index.radial_query(ra=targets_ra, dec=targets_dec, radius=2.0,
	                                     return_counts=True, return_distance=True)
	print("fields:", matches.dtype.names)
	print("counts per target:", counts)
	print("rows:", len(matches))
	for row in matches[:4]:
		print(f"  target {row['target_index']}  {row['ra']:9.5f} {row['dec']:+9.5f}  {row['distance']:.4f} deg")

.. code-block:: text

	fields: ('target_index', 'ra', 'dec', 'distance')
	counts per target: [18 14 18]
	rows: 50
	  target 0   99.23592  +8.64279  1.5526 deg
	  target 0   98.54506  +8.78621  1.8798 deg
	  target 0   98.39569  +9.51080  1.6551 deg
	  target 0   98.32358  +9.67745  1.6830 deg

``radius`` may be one value for every target or one value per target. The batch
form is about the bookkeeping, not about speed — it runs one query per target,
because the pixel ranges covering a cone depend on where the cone is and cannot be
shared. The method's own documentation gives the measured figure.

Note the interaction with qlsc's usual scalar-like rule: a size-1 array is a
scalar, so a target list with a single entry would be read as a scalar query. Set
``return_counts`` to ask for the batch contract regardless of shape.

For matching two *indices* against each other, rather than a target list against
one index, use :py:func:`~qlsc.crossmatch`.


Looking a point up by key
-------------------------

:py:func:`~qlsc.QLSCIndex.find_by_key` returns the coordinates carrying a key,
without any spatial query at all — the ``key`` column is indexed:

.. code-block:: python

	print(index.find_by_key("S29764"))

.. code-block:: text

	[(99.23592352, 8.64279335)]

Keys are not unique, so this returns every point carrying the key; an absent key
returns an empty array rather than raising.


Working with astropy
--------------------

Wherever a coordinate goes, a scalar :class:`~astropy.coordinates.SkyCoord` may be
given instead, as the ``ra`` argument with ``dec`` left unset; it is transformed to
ICRS first. An **array** SkyCoord is a target list, and so selects the batch form
above exactly as arrays of degrees do.

``return_skycoord`` asks for the matches back as an ICRS
:class:`~astropy.coordinates.SkyCoord` — always of length n, never squeezed.
Combined with any other flag the result is an :class:`~astropy.table.QTable` with a
``coord`` column plus ``target_index``, ``key`` and/or ``distance``, because a
SkyCoord cannot be a field of a recarray. Asking for either without astropy
installed raises ``ImportError``.


Exporting the coverage
----------------------

:py:func:`~qlsc.QLSCIndex.coverage` summarises which part of the sky an index
touches, as a sky-coverage v1 document: the occupied pixels at a chosen depth, each
as a polygon with the number of points in it.

.. code-block:: python

	cov = index.coverage(depth=3)
	print("keys:   ", sorted(cov))
	print("type:   ", cov["type"], "| version:", cov["version"], "| frame:", cov["frame"])
	print("polygons:", len(cov["polygons"]), "| first:", sorted(cov["polygons"][0]))
	print("meta:   ", {k: cov["meta"][k] for k in sorted(cov["meta"])})
	print("value:  ", cov["value"])

.. code-block:: text

	keys:    ['edges', 'frame', 'meta', 'polygons', 'type', 'units', 'value', 'version']
	type:    sky-coverage | version: 1 | frame: icrs
	polygons: 384 | first: ['id', 'value', 'vertices']
	meta:    {'depth': 3, 'producer': 'qlsc 1.2.0', 'scheme': 'qlsc'}
	value:   {'name': 'count', 'units': 'points'}

The document is plain JSON-compatible Python, so ``json.dump`` writes it straight
out. :py:func:`~qlsc.QLSC.coverage` is the same export for a bare set of pixels and
values, with no index behind it.

.. note::

	Cross-matching, coverage export, and the geographic classes are new in this
	release. Beyond this page and the :doc:`API reference <../api>`, they now have
	worked examples: :doc:`../examples/cross_match`, :doc:`../examples/coverage_export`,
	and :doc:`../examples/geographic_classes`.


How a region query works
------------------------

It is worth knowing what happens between the call and the answer, because it
explains several of the behaviours above.

1. The region is handed to the embedded Q3C code, which returns the ipix ranges
   that cover it — some fully inside, some straddling the boundary. **This is
   always computed at depth 30**, whatever the index's own depth, because the
   upstream code is only exercised at that resolution and produces degenerate
   ranges at small ``nside``.
2. Those ranges are shifted down to the index's depth, rounding *outward* so a
   partially covered pixel is still a candidate, then sorted and merged into
   disjoint intervals.
3. The intervals become a SQL ``WHERE`` clause of the form
   ``(ipix>=lo AND ipix<hi) OR ...``, which the database answers using its index
   on the ipix column.
4. Every candidate row returned is then tested **exactly** — against the cone
   radius with :py:func:`~qlsc.sindist`, or with the same containment test
   :py:func:`~qlsc.QLSC.point_in_polygon` and
   :py:func:`~qlsc.QLSC.point_in_ellipse` use.

Step 4 is why the answer never depends on the depth, why concave polygons are
safe (the coverage step is allowed to be generous; the exact test removes the
extras), and why a coarse index is slower rather than wrong.

Two things can go wrong in step 1, and both are handled:

* For a very small or zero-size region, Q3C sometimes produces no ranges at all.
  The cone and ellipse queries grow the region and retry — the region searched
  always contains the one you asked for, so step 4 still gives the right answer —
  and log a warning through the ``qlsc`` logger when they had to.
* At some positions and sizes Q3C runs out of room for the ranges it wants to
  emit and raises ``qlsc.q3c.RangeOverflowError``. The cone and ellipse queries
  handle this the same way, by widening. A polygon has no size to widen, so this
  one can reach you; see :doc:`../examples/survey_footprint`.

A cheaper, coarser cut
----------------------

Q3C also has a much cheaper calculation for the same shape of question, the one
its PostgreSQL extension uses to seed a positional join:
:py:func:`~qlsc.QLSC.nearby_ipix_ranges`. It skips step 1's quadtree refinement
entirely and returns **at most four** ranges, each a single square of the
quadtree at the level just coarse enough to hold the cone's bounding box on a
cube face:

.. code-block:: python

	import numpy as np
	from qlsc import QLSC, q3c

	q = QLSC(depth=20)
	nearby = q.nearby_ipix_ranges(83.822083, -5.391111, 0.5)
	print("nearby ranges:", nearby.tolist())
	print("pixels covered:", int(sum(e - s + 1 for s, e in nearby)))

	fulls, partials = q3c.radial_query(q._hprm, 83.822083, -5.391111, 0.5)
	refined = int(sum(e - s + 1 for s, e in np.vstack([fulls, partials])))
	print("radial_query ranges:", len(fulls) + len(partials), " pixels covered:", refined)
	print("whole sky past 35 deg:", q.nearby_ipix_ranges(83.822083, -5.391111, 40.).tolist(), "of", q.nbins, "pixels")

.. code-block:: text

	nearby ranges: [[2469606195200, 2469874630655], [2470143066112, 2470411501567]]
	pixels covered: 536870912
	radial_query ranges: 24  pixels covered: 100663320
	whole sky past 35 deg: [[0, 6597069766655]] of 6597069766656 pixels

Two ranges instead of twenty-four, covering five times as much sky. That is the
trade: it is a first cheap cut at "which pixels should I even look at?", not an
answer. Past ``qlsc.q3c.Q3C_MAXRAD`` (35 degrees) it stops doing geometry
altogether and hands back the whole sky.

.. warning::

	These ranges are not documented as a *guaranteed* cover of the cone. Where the
	cone spills onto two or three neighbouring cube faces, Q3C emits one square
	per neighbouring face, grown until it contains that face's whole bounding
	box. That construction argues coverage, and the package's tests
	measure it over a grid of centres, radii and depths — but neither is a proof.
	Anything that must not miss a point belongs in
	:py:func:`~qlsc.QLSCIndex.radial_query`, which is built on the refined ranges
	of step 1 and finishes with the exact test of step 4.


If the range computation ever produces nothing usable, the query falls back to
scanning every row rather than silently returning an empty result, and says so
through the logger:

.. code-block:: python

	import logging
	logging.basicConfig(level=logging.WARNING)

That is worth turning on the first time you run an unusual query.
