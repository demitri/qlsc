Coordinates and Conversions
===========================

.. contents::
   :local:
   :depth: 2

A position in this package has three representations, and :py:class:`~qlsc.QLSC`
converts between them:

1. **Sky coordinates** — right ascension and declination in degrees.
2. **The pixel number** — an integer ``ipix``, valid for one particular depth.
3. **Face coordinates** — a cube face number and an ``(x,y)`` position on that
   flat square face, both axes running from −1 to +1.

Angles are degrees throughout. If astropy is installed they may also be given as
:class:`~astropy.units.Quantity` values, and wherever an ``ra,dec`` pair is taken a
scalar :class:`~astropy.coordinates.SkyCoord` may be passed as ``ra`` alone
instead, with ``dec`` left unset; it is transformed to ICRS first.
:py:class:`~qlsc.GeoQLSC` and :py:class:`~qlsc.GeoQLSCIndex` are the same scheme and
index in geographic vocabulary — latitude for declination, and longitude reported
in [-180,180) rather than [0,360) — with :py:func:`~qlsc.wrap_longitude` to convert.

The third representation is the one to reach for when a calculation gets complicated. On a face
the geometry is a plane, not a sphere: pixel boundaries are straight lines at
regular intervals, "adjacent" means what you expect, and there is no wrap at
α = 360° and no pole. The source of
:py:func:`~qlsc.QLSC.ipix2polygon` is a short worked example of the idea.


Sky coordinates to pixel numbers
--------------------------------

:py:func:`~qlsc.QLSC.ang2ipix` takes scalars, a pair of arrays, or an ``(n,2)``
array of points, and returns the matching shape:

.. code-block:: python

	import numpy as np
	from qlsc import QLSC

	q = QLSC(depth=8)

	print("scalar:", q.ang2ipix(45., 60.), type(q.ang2ipix(45., 60.)).__name__)

	ipix = q.ang2ipix(ra=np.array([45., 0., 300.]), dec=np.array([60., 0., -20.]))
	print("array: ", ipix, ipix.dtype)

	print("points:", q.ang2ipix(points=np.array([[45., 60.], [0., 0.], [300., -20.]])))

.. code-block:: text

	scalar: 26010 int
	array:  [ 26010 114688 291017] int64
	points: [ 26010 114688 291017]

Arrays are passed straight through to the C layer, which requires
``dtype=double``; anything else is converted, which means a copy. For a catalog
of a few hundred million rows that is worth arranging in advance.

Right ascension outside ``[0,360]`` is wrapped without any copying. A declination
outside ``[-90,90]`` is *normalized*, not rejected — the point is reflected over
the pole and the right ascension shifted by 180°, which is the coordinate it
actually names:

.. code-block:: python

	print(q.ang2ipix(405., 60.) == q.ang2ipix(45., 60.))
	print(q.ang2ipix(0., 91.), q.ang2ipix(180., 89.))

.. code-block:: text

	True
	49160 49160

Non-finite values are refused rather than silently producing a garbage pixel:

.. code-block:: text

	ang2ipix(nan, 0) -> ValueError: The 'ra' and 'dec' values must be finite; was given ra=nan, dec=0.0.


Pixel numbers back to sky coordinates
-------------------------------------

:py:func:`~qlsc.QLSC.ipix2ang` — and its explicit synonym
:py:func:`~qlsc.QLSC.ipix2ang_center` — returns the **centre** of the pixel.
Both are vectorized in the same way as ``ang2ipix``: a scalar in gives two floats
out, an array in gives two arrays out.

.. code-block:: python

	ra, dec = q.ipix2ang(ipix)
	print(ra, dec)
	print(q.ipix2ang(int(ipix[0])))

.. code-block:: text

	[4.5000000e+01 2.2381050e-01 2.9986525e+02] [ 59.88418619   0.22380879 -19.92365872]
	(45.0, 59.88418619334755)

.. warning::

	Before qlsc v1.1, ``ipix2ang`` returned the pixel's "lower left" corner. It
	now returns the centre, matching the Q3C PostgreSQL extension from v2.0.1
	onward. This is a silent change of value, not an error, so code written
	against an older version needs checking. To get a corner, use
	:py:func:`~qlsc.QLSC.ipix2xy` or :py:func:`~qlsc.QLSC.ipix2polygon`.

The distinction is easy to see at a coarse depth:

.. code-block:: python

	q2 = QLSC(depth=2)
	print("ipix2ang_center(42):", q2.ipix2ang_center(42))
	print("ipix2xy(42):        ", q2.ipix2xy(42))
	for v in q2.ipix2polygon(ipix=42):
		print(f"   {v[0]:10.5f} {v[1]:+10.5f}")

.. code-block:: text

	ipix2ang_center(42): (53.13010235415598, 30.96375653207353)
	ipix2xy(42):         (2, -1.0, 0.5)
	     45.00000  +19.47122
	     63.43495  +24.09484
	     63.43495  +41.81031
	     45.00000  +35.26439

:py:func:`~qlsc.QLSC.ipix2polygon` returns the four corners in order, as
``ra,dec`` pairs; pass ``duplicate_endpoint=True`` if you want the first vertex
repeated at the end for drawing. Remember that the edges between those corners are
great circles, not lines of constant α or δ. It takes one ipix at a time, by
design: vectorizing it would still have to project the corners face by face.


The face plane
--------------

:py:func:`~qlsc.QLSC.ipix2xy`, :py:func:`~qlsc.QLSC.xy2ang` and
:py:func:`~qlsc.QLSC.xy2ipix` work in face coordinates that run from −1 to +1
across a face. ``ipix2xy`` gives the pixel's "lower left" corner; the pixel is
``2/nside`` wide in each direction.

.. code-block:: python

	print("xy2ang(0, 0, face 2):    ", q.xy2ang(facenum=2, x=0.0, y=0.0))
	print("xy2ang(1, 1, face 2):    ", q.xy2ang(facenum=2, x=1.0, y=1.0))
	print("xy2ang(0.5, 0.5, face 2):", q.xy2ang(facenum=2, x=0.5, y=0.5))
	print("xy2ipix(0.5, 0.5, face 2):", q.xy2ipix(facenum=2, x=0.5, y=0.5))
	print("ipix2xy of that:", q.ipix2xy(q.xy2ipix(facenum=2, x=0.5, y=0.5)))

.. code-block:: text

	xy2ang(0, 0, face 2):     (90.0, 0.0)
	xy2ang(1, 1, face 2):     (135.0, 35.264389682754654)
	xy2ang(0.5, 0.5, face 2): (116.56505117707799, 24.094842552110702)
	xy2ipix(0.5, 0.5, face 2): 192512
	ipix2xy of that: (2, 0.5, 0.5)

``(0,0)`` is the face centre, ``(±1,±1)`` are the cube's vertices, and the
transform is the tangent-plane one: ``x`` is the tangent of the angle from the
face centre along one axis. ``xy2ang`` also accepts a ``points`` array of
``(x,y)`` pairs for one face at a time.

.. caution::

	:py:func:`~qlsc.QLSC.ang2ipix_xy` also returns ``x`` and ``y``, but on a
	**half-scale**: its values run from −0.5 to +0.5 across a face, exactly half
	of what :py:func:`~qlsc.QLSC.ipix2xy` and :py:func:`~qlsc.QLSC.xy2ang` use.
	Feeding its output back into ``xy2ang`` gives a position half the angular
	distance from the face centre, with no error raised.

	.. code-block:: python

		q30 = QLSC(depth=30)
		for ra, dec in [(90., 0.), (134.999, 0.), (100., 10.)]:
			r = q30.ang2ipix_xy(ra=ra, dec=dec)
			print(f"  ang2ipix_xy({ra:8.3f},{dec:5.1f}) -> facenum={r['facenum']} "
			      f"x={r['x']:+.6f} y={r['y']:+.6f}")

		r = q30.ang2ipix_xy(ra=100., dec=10.)
		facenum, x, y = q30.ipix2xy(r["ipix"])
		print(f"  ipix2xy of the same pixel  -> facenum={facenum} x={x:+.6f} y={y:+.6f}")
		print(f"  ratio: {x/r['x']:.6f}, {y/r['y']:.6f}")

	.. code-block:: text

		  ang2ipix_xy(  90.000,  0.0) -> facenum=2 x=+0.000000 y=+0.000000
		  ang2ipix_xy( 134.999,  0.0) -> facenum=2 x=+0.499983 y=+0.000000
		  ang2ipix_xy( 100.000, 10.0) -> facenum=2 x=+0.088163 y=+0.089524
		  ipix2xy of the same pixel  -> facenum=2 x=+0.176327 y=+0.179047
		  ratio: 2.000000, 2.000000

	Multiply ``ang2ipix_xy``'s ``x`` and ``y`` by two to put them in the same
	system as everything else, or use ``ipix2xy`` on the ``ipix`` it also returns.

The projection of a face stays meaningful past the edge of the face itself, which
is how a calculation can step off one face onto its neighbour.
:py:func:`~qlsc.xy2facenum` says which face a point beyond the edge really
belongs to. :py:func:`~qlsc.QLSC.xy2ang`, in contrast, restricts its inputs to
the face:

.. code-block:: text

	xy2facenum(1.4, 0.0, 1): 2
	xy2facenum(0.4, 0.0, 1): 1
	xy2ang past the edge: ValueError: Values of 'x' and 'y' must be in the range [-1,1]; was given (1.4,0.0).


The integer grid
----------------

Underneath the face plane is the grid the scheme is actually built from: an ipix
value is the bit interleaving of a face number and a pair of **integers** in
``[0, nside)``. :py:func:`~qlsc.QLSC.ipix2xiyi` and
:py:func:`~qlsc.QLSC.xiyi2ipix` are that layer, and they are exact inverses of
each other — no trigonometry, no tolerance:

.. code-block:: python

	facenum, xi, yi = q.ipix2xiyi(9200)
	print("ipix2xiyi(9200):  ", (facenum, xi, yi), "  nside =", q.nside)
	print("ipix2xy(9200):    ", q.ipix2xy(9200))
	print("xiyi2ipix back:   ", q.xiyi2ipix(facenum, xi, yi))
	print("one pixel along x:", q.xiyi2ipix(facenum, xi + 1, yi))

	ipix = q30.ang2ipix(83.822083, -5.391111)
	print("depth 30:", ipix, "->", q30.ipix2xiyi(ipix))
	print("round trip is exact:", q30.xiyi2ipix(*q30.ipix2xiyi(ipix)) == ipix)

.. code-block:: text

	ipix2xiyi(9200):   (0, 28, 92)   nside = 256
	ipix2xy(9200):     (0, -0.78125, -0.28125)
	xiyi2ipix back:    9200
	one pixel along x: 9201
	depth 30: 2589827523372428891 -> (2, 478757421, 485909747)
	round trip is exact: True

``x`` and ``xi`` are the same position written two ways —
``x == (xi / nside) * 2 - 1`` exactly — but the integers are what to reach for
when a calculation walks the grid. Stepping ``xi`` or ``yi`` by one is an exact
move to the next pixel at any depth; nudging ``x`` or ``y`` and going back
around through ``(ra,dec)`` is an approximation, and it is least reliable
precisely where a grid walk spends its time, on a pixel boundary. The two
methods take one pixel at a time, and refuse a grid coordinate outside
``[0, nside)`` — that is a position on a *neighbouring* face, whose axes are
rotated with respect to this one, not a wrapped position on this one. Crossing a
face is what :py:func:`~qlsc.QLSC.ipix_neighbors` is for.

.. note::

	:py:func:`~qlsc.QLSC.xy2ipix` uses this path itself, rather than projecting
	the position back onto the sphere and calling :py:func:`~qlsc.QLSC.ang2ipix`,
	which would put a boundary position on whichever side the last few bits of
	two trigonometric round trips land on. ``xy2ipix(ipix2xy(ipix))`` returns
	``ipix`` at every depth.


Moving between depths, and between neighbours
---------------------------------------------

:py:func:`~qlsc.QLSC.ipix_up` and :py:func:`~qlsc.QLSC.ipix_down` step one level
of the quadtree; they are ``ipix >> 2`` and the four values from ``ipix << 2``.
:py:func:`~qlsc.QLSC.ipix_neighbors` returns the pixels touching a given pixel —
sharing an edge or a corner:

.. code-block:: python

	q6 = QLSC(depth=6)
	ipix = q6.ang2ipix(83.822083, -5.391111)
	print(ipix, q6.ipix_neighbors(ipix))
	print("at a cube corner:", len(q6.ipix_neighbors(q6.ang2ipix(45., 35.26438968))))
	q0 = QLSC(depth=0)
	print("depth 0, face 1:", q0.ipix_neighbors(1))

.. code-block:: text

	9200 [9167 9178 9179 9189 9191 9201 9202 9203]
	at a cube corner: 7
	depth 0, face 1: [0 2 4 5]

The result is sorted and de-duplicated, and never contains the pixel itself.
Three cases are geometry, not bugs: a pixel on a **corner of the cube** has seven
neighbours, because only three pixels meet at a cube vertex; at ``depth=0`` a face
has four (the four adjacent faces, correctly not the opposite one); and the
ordering is numeric rather than compass-directional, because two faces meeting at
an edge have their axes rotated with respect to each other, so "north" is not
well defined across the seam. If direction matters, take
:py:func:`~qlsc.QLSC.ipix2xy` of the results and compare positions.

The method takes one pixel at a time. :doc:`../examples/pixel_regions` uses it to
grow a region outward from a detection.


Containment tests
-----------------

Two methods answer "is this point inside that region?" without any index at all,
and their answers do not depend on the object's depth:

* :py:func:`~qlsc.QLSC.point_in_polygon` — a spherical polygon given as an
  ``(n,2)`` array of ``ra,dec`` vertices, with great-circle edges, in either
  winding order. It is the same test :py:func:`~qlsc.QLSCIndex.polygon_query`
  applies to its candidates. See :doc:`../examples/survey_footprint`.
* :py:func:`~qlsc.QLSC.point_in_ellipse` — an ellipse given by its centre,
  semi-major axis, axis ratio and position angle, measured from north through
  east. It is the test behind :py:func:`~qlsc.QLSCIndex.ellipse_query`. See
  :doc:`../examples/galaxy_aperture`.


.. _distances:

Distances
---------

Two module-level functions compute separations, both taking degrees:

.. code-block:: python

	from qlsc import distance, sindist
	import numpy as np

	# printed to 12 significant digits: the last digit or two of a full-precision
	# value depends on the platform's C math library (sin, cos, asin round
	# differently on macOS and Linux)
	print(f"distance((0,0),(0,1)):     {distance(0., 0., 0., 1.):.12g}")
	print(f"distance((0,89),(180,89)): {distance(0., 89., 180., 89.):.12g}")

	d = distance(10., 20., 11., 21.)
	s = sindist(10., 20., 11., 21.)
	print(f"distance = {d:.12g} deg")
	print(f"sindist  = {s:.12g}")
	print(f"sin^2(d/2) = {np.sin(np.radians(d)/2)**2:.12g}")
	print(f"2*asin(sqrt(sindist)) in degrees = {np.degrees(2*np.arcsin(np.sqrt(s))):.12g}")

.. code-block:: text

	distance((0,0),(0,1)):     1
	distance((0,89),(180,89)): 2
	distance = 1.37015343906 deg
	sindist  = 0.000142959314662
	sin^2(d/2) = 0.000142959314662
	2*asin(sqrt(sindist)) in degrees = 1.37015343906

:py:func:`~qlsc.distance` returns the angular separation in degrees.
:py:func:`~qlsc.sindist` returns the *haversine* of the separation,
:math:`\sin^2(d/2)` — not the sine of the separation. It is the quantity a cone
filter actually needs, because comparing it against
:math:`\sin^2(r/2)` decides membership with no inverse trigonometry; convert with
:math:`d = 2\arcsin\sqrt{\mathrm{sindist}}` when you want degrees. This is what
:py:func:`~qlsc.QLSCIndex.radial_query`'s ``return_distance`` does for you.
