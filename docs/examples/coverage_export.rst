Exporting Sky Coverage
======================

.. contents::
   :local:
   :depth: 2

**The problem.** You have built up an index — a survey footprint, a set of
detections, anything with sky positions in it — and now something *else* needs
to know what part of the sky it covers: a plotting library, a proposal tool, a
teammate's script in another language. You want a portable description of
"these polygons, each with this many points in it," not qlsc's own SQLite
schema.

:py:func:`~qlsc.QLSCIndex.coverage` answers exactly that, as a plain
JSON-compatible Python dictionary in the `sky-coverage v1
<https://github.com/demitri/astplot-js>`_ format that `astplot-js
<https://github.com/demitri/astplot-js>`_'s ``CoverageLayer`` renders as a
filled-polygon coverage map or choropleth.


A patchy survey
----------------

Three fields of different density, so the coverage map has something to show:

.. code-block:: python

	import json
	import numpy as np
	from qlsc import QLSC, QLSCIndex

	rng = np.random.default_rng(seed=20260830)

	index = QLSCIndex()
	for center_ra, center_dec, n, spread in [(40., -20., 20_000, 1.5),
	                                          (200., 10., 8_000, 0.8),
	                                          (210., 12., 3_000, 0.5)]:
		ra  = rng.normal(center_ra, spread, n) % 360.
		dec = np.clip(rng.normal(center_dec, spread, n), -89.9, 89.9)
		index.add_points(ra=ra, dec=dec)

	print("indexed points:", len(index))

.. code-block:: text

	indexed points: 31000


Exporting the coverage
------------------------

.. code-block:: python

	cov = index.coverage(depth=4)
	print("keys:", sorted(cov))
	print("type:", cov["type"], "| version:", cov["version"], "| frame:", cov["frame"],
	      "| edges:", cov["edges"], "| units:", cov["units"])
	print("occupied pixels:", len(cov["polygons"]))
	print("meta:", {k: cov["meta"][k] for k in sorted(cov["meta"])})
	print("value descriptor:", cov["value"])

.. code-block:: text

	keys: ['edges', 'frame', 'meta', 'polygons', 'type', 'units', 'value', 'version']
	type: sky-coverage | version: 1 | frame: icrs | edges: geodesic | units: deg
	occupied pixels: 17
	meta: {'depth': 4, 'producer': 'qlsc 1.1.0', 'scheme': 'qlsc'}
	value descriptor: {'name': 'count', 'units': 'points'}

Only pixels holding at least one point appear — 17 out of the 1,536 pixels at
depth 4 for the whole sky, which is what a "patchy" footprint should look
like. ``depth`` picks the resolution of the map: coarser than the index's own
depth (the default is 5, or the index's depth if that is coarser), since the
index cannot recover pixel detail finer than it stored.

Each polygon carries its ipix number, its four corners, and the point count
that fell in it, and the counts add up to the whole index exactly, because no
point is dropped or double-counted:

.. code-block:: python

	first = cov["polygons"][0]
	print("first polygon keys:", sorted(first))
	print("id:", first["id"], "| value:", first["value"], "| vertices:", len(first["vertices"]))
	print("a vertex:", first["vertices"][0])

	total_points = sum(p["value"] for p in cov["polygons"])
	print("sum of per-pixel counts:", total_points, "== indexed points:", len(index))

.. code-block:: text

	first polygon keys: ['id', 'value', 'vertices']
	id: 347 | value: 36.0 | vertices: 4
	a vertex: [32.005383208083494, -27.923542300530904]
	sum of per-pixel counts: 31000.0 == indexed points: 31000


It really is plain JSON
--------------------------

The whole point of the format is that nothing beyond the standard library is
needed to write it out or read it back:

.. code-block:: python

	path = "coverage.json"
	with open(path, "w") as f:
		json.dump(cov, f)

	with open(path) as f:
		reloaded = json.load(f)

	print("round-trips through plain json:", reloaded == cov)
	print("polygons read back:", len(reloaded["polygons"]))

.. code-block:: text

	round-trips through plain json: True
	polygons read back: 17


A depth cannot outrun the index
-----------------------------------

An index only ever stores the pixel each point falls in at its *own* depth, so
asking for a finer coverage depth than that raises rather than silently
falling back to the index's depth:

.. code-block:: python

	try:
		index.coverage(depth=31)
	except ValueError as e:
		print(type(e).__name__, e)

.. code-block:: text

	ValueError The coverage depth (31) is finer than this index's depth (30); the index does not record where within its own pixels its points fall, so a finer pixel cannot be recovered.

The same happens for a coarser index built at a shallower depth on purpose:

.. code-block:: python

	coarse_index = QLSCIndex(qlsc=QLSC(depth=6))
	coarse_index.add_points(ra=[10., 11., 200.], dec=[5., 5., -30.])
	try:
		coarse_index.coverage(depth=7)
	except ValueError as e:
		print(type(e).__name__, e)

.. code-block:: text

	ValueError The coverage depth (7) is finer than this index's depth (6); the index does not record where within its own pixels its points fall, so a finer pixel cannot be recovered.


Pinning the color scale, and adding your own metadata
---------------------------------------------------------

``value_domain`` lets several documents share one color scale — useful for a
sequence of maps a viewer will flip between — and ``meta`` records whatever a
producer wants to keep attached, short of the three keys qlsc writes itself:

.. code-block:: python

	custom = index.coverage(depth=3, value_domain=[0, 50_000],
	                        meta={"survey": "demo-2026a"})
	print("value:", custom["value"])
	print("meta:", {k: custom["meta"][k] for k in sorted(custom["meta"])})

	try:
		index.coverage(depth=3, meta={"scheme": "not-qlsc"})
	except ValueError as e:
		print(type(e).__name__, e)

.. code-block:: text

	value: {'name': 'count', 'units': 'points', 'domain': [0.0, 50000.0]}
	meta: {'depth': 3, 'producer': 'qlsc 1.1.0', 'scheme': 'qlsc', 'survey': 'demo-2026a'}
	ValueError These 'meta' keys are written by qlsc itself and cannot be replaced: ['scheme'].

A document claiming to have been produced by something other than the qlsc
version that actually produced it would be worse than no provenance at all, so
the collision is refused rather than letting one value quietly win.


Exporting a bare set of pixels, no index required
------------------------------------------------------

:py:func:`~qlsc.QLSC.coverage` is the same export with no index behind it —
useful when the pixels and their values already exist, computed some other
way, and only need to become a document:

.. code-block:: python

	q = QLSC(depth=5)
	some_pixels = [1000, 1001, 1002, 6000]
	some_values = [3.5, 1.2, 9.9, 0.4]
	bare = q.coverage(ipix=some_pixels, values=some_values,
	                  value_name="brightness", value_units="mag")
	print("polygons:", len(bare["polygons"]), "| depth in meta:", bare["meta"]["depth"])
	print("value:", bare["value"])
	print([(p["id"], p["value"]) for p in bare["polygons"]])

.. code-block:: text

	polygons: 4 | depth in meta: 5
	value: {'name': 'brightness', 'units': 'mag'}
	[(1000, 3.5), (1001, 1.2), (1002, 9.9), (6000, 0.4)]

The ipix numbers need not be unique or sorted, and every one becomes its own
polygon entry in the order given — nothing is merged, which is why the same
pixel listed twice would produce two overlapping polygons rather than an
error.


See also
--------

* :py:func:`~qlsc.QLSCIndex.coverage` and :py:func:`~qlsc.QLSC.coverage` — the
  full method documentation.
* `astplot-js <https://github.com/demitri/astplot-js>`_ — the ``CoverageLayer``
  consumer this format was designed for.
* :doc:`geographic_classes` — a :py:class:`~qlsc.GeoQLSCIndex`'s coverage
  documents still declare right ascension in [0,360), the one place that class
  does not report its own [-180,180) convention.
