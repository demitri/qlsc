Cross-Matching Two Catalogs
============================

.. contents::
   :local:
   :depth: 2

**The problem.** You have run your own detection pipeline over a patch of sky,
and a reference catalog already covers the same field. Before doing anything
else with your detections, you need to know which of them are already known —
and, just as importantly, which ones are not.

This is a spatial join of two catalogs: for every point in one, find the
points in the other within some radius. :py:func:`~qlsc.crossmatch` does this
between two :py:class:`~qlsc.QLSCIndex` objects.


The two catalogs
-----------------

A reference catalog, and a set of detections that mostly — but not entirely —
reproduce it, with realistic astrometric noise:

.. code-block:: python

	import numpy as np
	from qlsc import QLSCIndex, crossmatch

	rng = np.random.default_rng(seed=20260222)

	# a reference catalog, spread over one region of sky
	n_ref = 6_000
	ref_ra  = rng.uniform(150., 152., n_ref)
	ref_dec = rng.uniform(-11., -9., n_ref)
	ref_keys = [f"REF-{i:05d}" for i in range(n_ref)]

	# our own detections: most are the reference sources with 0.3" of astrometric
	# noise, plus a few thousand sources the reference catalog does not have
	n_shared = 4_000
	shared_ra  = ref_ra[:n_shared]  + rng.normal(0., 0.3/3600., n_shared)
	shared_dec = ref_dec[:n_shared] + rng.normal(0., 0.3/3600., n_shared)

	n_new = 1_500
	new_ra  = rng.uniform(150., 152., n_new)
	new_dec = rng.uniform(-11., -9., n_new)

	det_ra  = np.concatenate([shared_ra, new_ra])
	det_dec = np.concatenate([shared_dec, new_dec])
	det_keys = [f"DET-{i:05d}" for i in range(n_shared)] + \
	           [f"DET-{i:05d}" for i in range(n_shared, n_shared + n_new)]

	print("reference catalog:", n_ref, "| our detections:", len(det_ra))

.. code-block:: text

	reference catalog: 6000 | our detections: 5500

Note that ``det_ra``/``det_dec`` hold every detection whether or not it has a
counterpart — ``crossmatch`` does not need the two catalogs pre-separated into
"known" and "new"; that separation is exactly what it computes.


Building the indices
---------------------

.. check-outputs: setup -- builds the indices the rest of the page queries

.. code-block:: python

	reference = QLSCIndex()
	reference.add_points(ra=ref_ra, dec=ref_dec, keys=ref_keys)

	detections = QLSCIndex()
	detections.add_points(ra=det_ra, dec=det_dec, keys=det_keys)


Finding the matches
--------------------

.. code-block:: python

	pairs = crossmatch(detections, reference, radius=1.0/3600., nearest_only=True)
	print("matched pairs:", len(pairs))
	print("fields:", pairs.dtype.names)
	for row in pairs[:5]:
		print(f"  {row.key_a}  {row.key_b}  {row.distance*3600:.3f}\"")

.. code-block:: text

	matched pairs: 3989
	fields: ('ra_a', 'dec_a', 'key_a', 'ra_b', 'dec_b', 'key_b', 'distance')
	  DET-02299  REF-02299  0.258"
	  DET-03867  REF-03867  0.521"
	  DET-00446  REF-00446  0.339"
	  DET-00763  REF-00763  0.475"
	  DET-00392  REF-00392  0.457"

``index_a`` (``detections`` here) is the catalog the result is grouped by:
each row pairs one detection with a counterpart from ``index_b``
(``reference``). ``nearest_only=True`` keeps only the closest counterpart of
each detection, which is what "does this detection already have an entry in
the reference catalog" needs.

The 3,989 matches are short of the 4,000 built-in overlap by eleven — at 0.3″
of per-axis noise a handful of shared pairs land more than 1″ apart by chance,
so they fall outside this particular radius. That is a real effect of the
match radius, not a defect in the join.


The other side of the join
---------------------------

A point of ``index_a`` with no counterpart is simply absent from the result,
which is what makes the count of *unmatched* detections easy to get by
subtraction:

.. code-block:: python

	unmatched = len(detections) - len(pairs)
	print("detections with no reference counterpart within 1\":", unmatched)
	print("expected around n_new:", n_new)

.. code-block:: text

	detections with no reference counterpart within 1": 1511
	expected around n_new: 1500

The eleven shared pairs that missed the radius above show up here too, on the
other side of the same effect: they count as "unmatched" along with the 1,500
genuinely new detections.


nearest_only and duplicate counterparts
-----------------------------------------

Without ``nearest_only``, a detection sitting in a locally crowded part of the
reference catalog can appear more than once — once per reference point within
the radius:

.. code-block:: python

	all_pairs = crossmatch(detections, reference, radius=1.0/3600., nearest_only=False)
	print("pairs without nearest_only:", len(all_pairs))
	print("pairs with nearest_only:  ", len(pairs))

	counts = {}
	for key in all_pairs["key_a"]:
		counts[key] = counts.get(key, 0) + 1
	multi = sorted(k for k, n in counts.items() if n > 1)
	print("detections with more than one reference counterpart:", len(multi))

.. code-block:: text

	pairs without nearest_only: 3991
	pairs with nearest_only:   3989
	detections with more than one reference counterpart: 2

Two detections happen to have two reference sources within 1″ at this
density; ``nearest_only`` collapses each of those pairs to its single closest
counterpart, which is why the ``nearest_only`` count is *not* simply
"3991 minus the extra rows" — it removes one row from each of the two crowded
detections, so the count drops by exactly two.


Cross-matching a catalog with itself
--------------------------------------

Every point of a catalog is trivially within any positive radius of itself,
which matters if you ever cross-match an index against itself — for a
duplicate check, say. With ``nearest_only``, that self-pair is always the
match, because nothing can be closer than zero:

.. code-block:: python

	same = crossmatch(reference, reference, radius=1.0/3600., nearest_only=True)
	print("self cross-match pairs:", len(same), "vs points:", len(reference))
	print("all at distance 0:", bool(np.all(same["distance"] == 0.)))

.. code-block:: text

	self cross-match pairs: 6000 vs points: 6000
	all at distance 0: True

A duplicate-detection search over one catalog therefore wants
``nearest_only=False`` and a filter that discards the zero-distance self-pairs
from the result.


The join does not depend on which side you iterate
------------------------------------------------------

``crossmatch`` always iterates over whichever index holds fewer points, for
speed, but the *set of matched pairs* — with ``nearest_only`` left off — does
not depend on which catalog was passed as ``index_a``:

.. code-block:: python

	forward = crossmatch(detections, reference, radius=1.0/3600., nearest_only=False)
	backward = crossmatch(reference, detections, radius=1.0/3600., nearest_only=False)
	print("forward pairs:", len(forward), "| backward pairs:", len(backward))

	forward_set = set(zip(forward["key_a"], forward["key_b"]))
	backward_set = set(zip(backward["key_b"], backward["key_a"]))
	print("same set of pairs either direction:", forward_set == backward_set)

.. code-block:: text

	forward pairs: 3991 | backward pairs: 3991
	same set of pairs either direction: True

This is the full join's symmetry, and it is worth keeping separate from
``nearest_only``: reducing "every pair within the radius" down to "the closest
counterpart of each ``index_a`` point" is a per-``index_a`` operation, so
*that* reduction, unlike the underlying join, does depend on which catalog you
call ``index_a``.


Edge cases
----------

An empty possible-match set returns an empty table rather than raising or
returning ``None``, and ``return_keys=False`` drops the key fields when the
identifiers are not needed:

.. code-block:: python

	empty_a = QLSCIndex()
	empty_a.add_point(0., 0.)
	none_nearby = crossmatch(empty_a, reference, radius=1.0/3600.)
	print("no possible match:", len(none_nearby))

	no_key_pairs = crossmatch(detections, reference, radius=1.0/3600.,
	                          nearest_only=True, return_keys=False)
	print("fields with return_keys=False:", no_key_pairs.dtype.names)

.. code-block:: text

	no possible match: 0
	fields with return_keys=False: ('ra_a', 'dec_a', 'ra_b', 'dec_b', 'distance')


See also
--------

* :py:func:`~qlsc.crossmatch` — the full method documentation, including its
  measured throughput and the memory cost of the smaller catalog.
* :doc:`../user_guide/indexing` — :py:func:`~qlsc.QLSCIndex.radial_query`'s
  batch form, which ``crossmatch`` is built on, for matching a fixed target
  list against one index rather than two indices against each other.
* :doc:`cone_search` — the single-cone version of the same underlying search.
