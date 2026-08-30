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
	file size on disk: 20480 bytes
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

	try:
		42 in index
	except TypeError as e:
		print(type(e).__name__, e)

.. code-block:: text

	TypeError Index keys are strings (the underlying column has TEXT affinity); was given '42' of type 'int'.

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
``key``, ``distance``:

.. code-block:: python

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

	The plain ``(n,2)`` form is squeezed, so a **single match comes back with
	shape (2,)**, not ``(1,2)``:

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

If the range computation ever produces nothing usable, the query falls back to
scanning every row rather than silently returning an empty result, and says so
through the logger:

.. code-block:: python

	import logging
	logging.basicConfig(level=logging.WARNING)

That is worth turning on the first time you run an unusual query.
