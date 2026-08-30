Working at the Pixel Level
==========================

.. contents::
   :local:
   :depth: 2

**The problem.** Sometimes the region you want is not a cone or a polygon but "the
pixels around this one" — a mask grown outward from a bad column, a set of tiles
to queue for follow-up around a detection, the neighbourhood of a pixel that came
back anomalous in a coverage map. Working in pixels rather than in angles means
integer set operations instead of geometry, which is both faster and easier to
reason about.

This example uses :py:func:`~qlsc.QLSC.ipix_neighbors`,
:py:func:`~qlsc.QLSC.ipix_up` and :py:func:`~qlsc.QLSC.ipix_down`, and ends with
the one trap in doing analysis this way: pixels are not all the same size.


The neighbours of a pixel
-------------------------

.. code-block:: python

	import numpy as np
	from qlsc import QLSC, QLSCIndex, distance

	q = QLSC(depth=6)
	ipix = q.ang2ipix(83.822083, -5.391111)
	neighbors = q.ipix_neighbors(ipix)

	print(f"depth 6, ipix {ipix} has {len(neighbors)} neighbours: {neighbors}")
	ra_c, dec_c = q.ipix2ang_center(ipix)
	for n in neighbors:
		r, d = q.ipix2ang_center(n)
		print(f"  {n:6d}  centre ({r:8.4f},{d:+8.4f})  {distance(ra_c, dec_c, r, d):.4f} deg away")

.. code-block:: text

	depth 6, ipix 9200 has 8 neighbours: [9167 9178 9179 9189 9191 9201 9202 9203]
	    9167  centre ( 81.9953, -7.9277)  2.4551 deg away
	    9178  centre ( 83.7581, -7.9579)  1.7527 deg away
	    9179  centre ( 85.5328, -7.9807)  2.5008 deg away
	    9189  centre ( 81.9953, -6.1816)  1.7527 deg away
	    9191  centre ( 81.9953, -4.4238)  2.5008 deg away
	    9201  centre ( 85.5328, -6.2231)  1.7644 deg away
	    9202  centre ( 83.7581, -4.4408)  1.7644 deg away
	    9203  centre ( 85.5328, -4.4536)  2.4880 deg away

Eight neighbours — four sharing an edge at about 1.76°, four sharing only a
corner at about 2.49°, which is the √2 you would expect from a square grid. Note
that the pixel numbers are *not* consecutive: 9167 and 9178 are in a different
branch of the quadtree even though they are physically adjacent. Z-order keeps
neighbours numerically close, but not contiguous.

Two special cases are geometry, not bugs:

.. code-block:: python

	corner = q.ang2ipix(45., 35.26438968)     # a vertex of the cube
	print(f"cube-corner pixel {corner}: {len(q.ipix_neighbors(corner))} neighbours")

	q0 = QLSC(depth=0)
	print(f"depth 0, face 1: {len(q0.ipix_neighbors(1))} neighbours -> {q0.ipix_neighbors(1)}")

.. code-block:: text

	cube-corner pixel 10922: 7 neighbours
	depth 0, face 1: 4 neighbours -> [0 2 4 5]

Three faces meet at a cube vertex, so only three pixels meet there — the diagonal
step across the vertex reaches no eighth pixel, and a corner pixel has seven
neighbours. At ``depth=0`` each face is a single pixel with four neighbours: the
four adjacent faces. Face 1's opposite, face 3, is correctly not among them.

The returned array is sorted and de-duplicated. It is *not* ordered by compass
direction, and could not be: two faces meeting at an edge have their coordinate
axes rotated relative to each other, so there is no consistent "north" to order
by. If you need direction, use :py:func:`~qlsc.QLSC.ipix2xy` on the results and
compare positions on the face.


Growing a region
----------------

Repeatedly taking the neighbours of everything in a set grows it outward one ring
at a time. This is breadth-first search on the pixel grid:

.. code-block:: python

	q8 = QLSC(depth=8)
	centre_ra, centre_dec = 83.822083, -5.391111
	seed = int(q8.ang2ipix(centre_ra, centre_dec))

	region = {seed}
	frontier = {seed}
	mean_area = 4 * np.pi / q8.nbins
	pixel_side = np.sqrt(mean_area) * 180. / np.pi
	print(f"  mean pixel side {pixel_side:.4f} deg")
	print("  ring   pixels   area/deg^2   max radius/deg")

	for ring in range(1, 6):
		new = set()
		for p in frontier:
			new.update(int(x) for x in q8.ipix_neighbors(p))
		frontier = new - region
		region |= frontier
		ra_r, dec_r = q8.ipix2ang_center(np.array(sorted(region)))
		radius = max(distance(centre_ra, centre_dec, a, b) for a, b in zip(ra_r, dec_r))
		area = len(region) * mean_area * (180./np.pi)**2
		print(f"  {ring:4d}   {len(region):6d}   {area:10.4f}   {radius:12.4f}")

.. code-block:: text

	  mean pixel side 0.3239 deg
	  ring   pixels   area/deg^2   max radius/deg
	     1        9       0.9442         0.8475
	     2       25       2.6228         1.4741
	     3       49       5.1407         2.1004
	     4       81       8.4978         2.7263
	     5      121      12.6943         3.3517

The pixel counts are exactly ``(2*ring + 1)**2``: the region is a square patch
of the cube face, growing by one row and column on every side each time. That is
the shape to expect as long as the growth stays on one face; crossing a face
boundary bends it.

Region growing is worth reaching for when the shape you want is defined by the
data rather than by geometry — flood-filling outward from a seed until some
per-pixel condition stops being true, for instance, which a cone search cannot
express at all.


Using a pixel set as a selection
--------------------------------

A set of ipix values is a selection. Compute the ipix of every catalog row at the
same depth and test membership; no query needed:

.. code-block:: python

	rng = np.random.default_rng(seed=20200422)
	n = 2_000_000
	ra  = rng.uniform(0., 360., n)
	dec = np.degrees(np.arcsin(rng.uniform(-1., 1., n)))

	# rebuild the region, stopping at ring 3 (49 pixels)
	region = {seed}
	frontier = {seed}
	for ring in range(3):
		new = set()
		for p in frontier:
			new.update(int(x) for x in q8.ipix_neighbors(p))
		frontier = new - region
		region |= frontier

	ipix = np.asarray(q8.ang2ipix(ra=ra, dec=dec))
	member = np.isin(ipix, np.fromiter(sorted(region), dtype=np.int64))
	print("region pixels:", len(region))
	print("catalogue objects in the region:", int(member.sum()))

.. code-block:: text

	region pixels: 49
	catalogue objects in the region: 442

This is the same idea a database uses when it stores the ipix alongside each row:
selection by integer comparison, at whatever resolution you chose.


.. _pixel-area-trap:

The trap: pixels are not the same size
--------------------------------------

Now check that answer against what the area predicts. The catalog has 48.5
objects per square degree, and 49 pixels at depth 8 with the mean pixel area of
0.1049 square degrees is 5.14 square degrees, so about 249 objects.

The measurement said 442.

The mean is the wrong number. As :ref:`explained in the user guide
<how-equal-are-the-pixels>`, Q3C's pixels vary in area by up to a factor of about
5.2, largest at the centre of a cube face and smallest at a cube vertex — and this
region happens to sit close to the centre of face 2:

.. code-block:: python

	print("seed pixel (face, x, y):", q8.ipix2xy(seed))

.. code-block:: text

	seed pixel (face, x, y): (2, -0.109375, -0.1015625)

:py:func:`~qlsc.QLSC.ipix_area` gives each pixel's true solid angle, so summing it
over the 49 pixels is a one-liner:

.. code-block:: python

	total = sum(q8.ipix_area(p) for p in sorted(region))
	deg2 = (180./np.pi)**2

	print(f"true region area:   {total*deg2:.4f} sq deg")
	print(f"mean-area estimate: {len(region)*mean_area*deg2:.4f} sq deg")
	print(f"ratio: {total/(len(region)*mean_area):.4f}")
	print(f"expected objects:   {total*deg2*n/41252.96:.0f}")

.. code-block:: text

	true region area:   9.5148 sq deg
	mean-area estimate: 5.1407 sq deg
	ratio: 1.8509
	expected objects:   461

461 predicted against 442 counted — well within the counting noise of a Poisson
number of that size, and a factor of 1.85 away from what the mean pixel area
would have told you.


The lesson is narrow but sharp: **counting objects in pixels is exact; converting
those counts to a surface density is not, unless you use the pixel's own area.**
Everything the index does — cone, polygon and ellipse queries — is unaffected,
because those filter their candidates exactly rather than treating a pixel as a
unit of area.


Moving between resolutions
--------------------------

A pixel's four children at the next depth, and its parent at the previous one, are
a bit shift away:

.. code-block:: python

	q5 = QLSC(depth=5)
	parent = q5.ang2ipix(83.822083, -5.391111)
	children = q5.ipix_down(parent)
	print(f"depth 5 pixel {parent} contains depth 6 pixels {children}")
	print(f"  and is itself inside depth 4 pixel {q5.ipix_up(parent)}")

	q6 = QLSC(depth=6)
	for c in children:
		r, d = q6.ipix2ang_center(c)
		print(f"    {c}: centre ({r:.4f},{d:+.4f})")

.. code-block:: text

	depth 5 pixel 2300 contains depth 6 pixels [9200, 9201, 9202, 9203]
	  and is itself inside depth 4 pixel 575
	    9200: centre (83.7581,-6.2052)
	    9201: centre (85.5328,-6.2231)
	    9202: centre (83.7581,-4.4408)
	    9203: centre (85.5328,-4.4536)

Note that ``ipix_up`` and ``ipix_down`` change the *value*, not the object: the
returned children are ipix values in the ``depth+1`` scheme, so you need a
``QLSC(depth=6)`` to interpret them, as above. This is the same shift the
:doc:`density map <density_map>` uses to roll fine counts up into coarse ones.


See also
--------

* :doc:`../user_guide/depth` — pixel sizes at every depth, and how to choose one.
* :ref:`how-equal-are-the-pixels` — why the areas vary, with the numbers.
* :doc:`density_map` — pixel counts used as a map.
