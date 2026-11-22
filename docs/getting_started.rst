Getting Started
===============

.. contents::
   :local:
   :depth: 2


Installing
----------

The package is on `PyPI <https://pypi.org/project/qlsc/>`_:

.. code-block:: shell

	pip install qlsc

or from a clone of the `repository <https://github.com/demitri/qlsc>`_:

.. code-block:: shell

	git clone https://github.com/demitri/qlsc.git
	cd qlsc/source
	pip install .

QLSC is a C extension: the spherical cube arithmetic is the Q3C C code, compiled
into the module ``qlsc.q3c`` and wrapped by the Python classes. Installing
therefore compiles that code, which needs a C compiler and the NumPy headers.
Both are handled automatically by pip's build isolation — the source distribution
declares its build requirements — so in practice ``pip install qlsc`` works on any
machine with a working compiler toolchain. Binary wheels are built for the common
platforms as well: where one matches yours, pip installs it and no compiler is
needed at all, and where none does, pip falls back to building from source.

NumPy is the only hard dependency. `astropy <https://www.astropy.org>`_ and
`pandas <https://pandas.pydata.org>`_ are optional: if astropy is present, angles
and radii may be given as :class:`~astropy.units.Quantity` objects anywhere the
API takes an angle in degrees.

Check the installation:

.. code-block:: python

	import qlsc
	print(qlsc.__version__)

.. code-block:: text

	1.2.0


The two classes
---------------

Everything in the package hangs off two objects.

:py:class:`~qlsc.QLSC` is the pixellation scheme itself: pure geometry, no data.
It converts between sky positions, pixel numbers, and positions on a cube face.

.. code-block:: python

	from qlsc import QLSC

	q = QLSC(depth=2)
	print(q)
	print(q.nbins)               # how many pixels this scheme has
	print(q.ang2ipix(45, 45))    # ra,dec -> pixel number
	print(q.ipix2ang_center(5))  # pixel number -> ra,dec of the pixel centre

.. code-block:: text

	<qlsc.QLSC at 0x10619ad50, depth=2>
	96
	5
	(45.0, 43.31385665828306)

:py:class:`~qlsc.QLSCIndex` is a searchable collection of points that uses a
scheme. You put coordinates in and ask spatial questions of them.

Both are covered in depth in the :doc:`user guide <user_guide/index>`; the rest of
this page is the short version.


Your first cone search
----------------------

Build an index. With no arguments it is created in memory at the finest
resolution the scheme supports (``depth=30``):

.. check-outputs: setup -- introduces the imports and the index the blocks below use

.. code-block:: python

	from qlsc import QLSCIndex
	from qlsc.generate import sunflower_points_on_sphere

	index = QLSCIndex()

Add points. :py:func:`~qlsc.generate.sunflower_points_on_sphere` produces any
number of roughly evenly spaced points on a sphere as an ``(n,2)`` array of
``ra,dec`` in degrees, which is one of the forms
:py:func:`~qlsc.QLSCIndex.add_points` accepts:

.. code-block:: python

	points = sunflower_points_on_sphere(n=1e6)
	print(points.shape)

	index.add_points(points=points)
	print(len(index))

.. code-block:: text

	(1000000, 2)
	1000000

Now ask what lies within half a degree of a position. Radii and coordinates are in
degrees:

.. code-block:: python

	matches = index.radial_query(ra=56., dec=-40., radius=0.5)
	print(type(matches).__name__, matches.shape)
	print(matches[:3])

.. code-block:: text

	ndarray (19, 2)
	[[ 55.97517129 -39.62402338]
	 [ 56.23910081 -39.71483435]
	 [ 56.50303028 -39.80576503]]

That is the whole search. One million points went in, nineteen came back, and
nothing scanned the other 999,981 rows: the query was turned into a small set of
pixel-number ranges first. See :doc:`examples/cone_search` for what that costs
in wall-clock time compared to the obvious NumPy loop.

Ask for the separations as well and the result becomes a
:class:`numpy.recarray` with named fields, which can be sorted by distance:

.. code-block:: python

	matches = index.radial_query(ra=56., dec=-40., radius=0.5, return_distance=True)
	matches.sort(order="distance")

	for row in matches[:5]:
		print(f"{row['ra']:10.5f} {row['dec']:+10.5f}   {row['distance']*3600:8.1f} arcsec")

.. code-block:: text

	  55.91287  -40.00953      242.7 arcsec
	  56.07598  -39.86202      539.2 arcsec
	  56.17680  -40.10085      607.6 arcsec
	  55.74975  -40.15735      892.2 arcsec
	  56.01368  -40.24888      896.7 arcsec

The distance costs nothing to ask for — the cone filter has already computed it
to decide whether each candidate is inside.

If astropy is installed, the radius can be a :class:`~astropy.units.Quantity`:

.. code-block:: python

	from astropy import units as u
	print(len(index.radial_query(ra=56., dec=-40., radius=30*u.arcmin)))

.. code-block:: text

	19


Keeping the index on disk
-------------------------

An in-memory index disappears when the program ends. Give the constructor a file
path and it persists, so a second script — or the same script run again — can open
it without rebuilding:

.. code-block:: python

	index = QLSCIndex(filepath="my_catalog.qlscidx")

	index.add_point(56.59625, 78.00866, key="NGC 1234")
	index.add_point(56.60000, 78.01000, key="NGC 1234b")
	print(len(index))
	index.close()

	# ... later, in another process ...
	index = QLSCIndex(filepath="my_catalog.qlscidx")
	print(len(index))
	print("NGC 1234" in index)

.. code-block:: text

	2
	2
	True

The file name and extension can be anything; ``.qlscidx`` is suggested so that
you remember what the file is. Any point can carry a ``key``, an arbitrary string
that comes back with the match when you ask for it:

.. code-block:: python

	hits = index.radial_query(ra=56.6, dec=78.0, radius=1./60.,
	                          return_key=True, return_distance=True)
	for row in hits:
		print(f"  {row['key']:>10}  {row['distance']*3600:6.2f} arcsec")

.. code-block:: text

	   NGC 1234b   36.00 arcsec
	    NGC 1234   31.30 arcsec

Keys make the index behave much like a dictionary — ``len(index)``,
``key in index``, ``index.keys()``, ``index.remove_point(key=...)``. See
:doc:`user_guide/indexing` for the details, including the ways in which it is
*not* a dictionary.


Choosing a resolution
---------------------

``QLSCIndex()`` with no scheme uses ``depth=30``, the finest available and the
value the Q3C PostgreSQL extension uses. A coarser scheme stores smaller pixel
numbers and produces fewer, wider ranges per query:

.. code-block:: python

	index = QLSCIndex(qlsc=QLSC(depth=10))
	index.add_points(points=sunflower_points_on_sphere(n=1e5))
	print(index.qlsc.depth, len(index))

.. code-block:: text

	10 100000

The answers a query returns do not depend on the depth — every region query is
filtered exactly, whatever the pixel size — so the depth is a performance and
storage choice, not a correctness one. :doc:`user_guide/depth` gives the pixel
sizes at every depth and the reasoning for picking one.


Where to next
-------------

* :doc:`user_guide/scheme` — what the pixels actually are, and why the numbering
  matters.
* :doc:`user_guide/indexing` — everything :py:class:`~qlsc.QLSCIndex` can do.
* :doc:`examples/index` — six worked problems with real output.
