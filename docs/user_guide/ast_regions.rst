.. check-outputs-page: needs-cornish -- every block needs the optional 'cornish' package, which is not installed with qlsc

Regions from Starlink AST
=========================

:py:func:`~qlsc.QLSCIndex.radial_query`, :py:func:`~qlsc.QLSCIndex.polygon_query` and
:py:func:`~qlsc.QLSCIndex.ellipse_query` offer the three shapes Q3C knows how to cover.
Real footprints are not always one of the three: a survey field with a chip gap, an
image outline read from a FITS header, "this pointing but not that one".

`cornish <https://github.com/demitri/cornish>`_ — a Pythonic interface to the Starlink
AST library — describes all of those, and :py:func:`~qlsc.QLSCIndex.region_query` puts
them in front of the index. It is an **optional** dependency: qlsc does not install it,
does not import it until one of these two methods is called, and says so plainly if it
is missing.

.. code-block:: bash

	pip install cornish

Any region AST can describe
---------------------------

.. code-block:: python

	import numpy as np
	from qlsc import QLSC, QLSCIndex
	from qlsc.generate import sunflower_points_on_sphere
	from cornish import ASTCircle, ASTBox, ASTICRSFrame, ASTCompoundRegion, RegionOperation

	index = QLSCIndex(qlsc=QLSC(depth=15))
	index.add_points(points=sunflower_points_on_sphere(n=1e5))

	frame = ASTICRSFrame()
	circle = ASTCircle(frame=frame, center=[45., 12.], radius=2.0)

	matches = index.region_query(circle)

	print(len(matches), "sources in the circle")
	print(np.round(matches[:3], 4))

.. code-block:: text

	28 sources in the circle
	[[43.813  10.8787]
	 [44.5039 10.6069]
	 [44.24   11.3189]]

The point of going through AST is not the circle — that one is
:py:func:`~qlsc.QLSCIndex.radial_query` in a hat — but everything that is not a circle.
Regions combine:

.. code-block:: python

	box = ASTBox.fromCorners(frame=frame, corners=([46., 13.], [50., 17.]))
	footprint = ASTCompoundRegion(regions=[box, circle], operation=RegionOperation.OR)

	print(len(index.region_query(circle)), "sources in the circle")
	print(len(index.region_query(box)), "sources in the box")
	print(len(index.region_query(footprint)), "sources in the union of the two")

.. code-block:: text

	28 sources in the circle
	37 sources in the box
	65 sources in the union of the two

The sources on an image
-----------------------

The question a catalogue is asked most often has its own method.
:py:func:`~qlsc.QLSCIndex.image_query` takes a FITS header, lets AST read its world
coordinate system, and queries the sky region bounding the field. qlsc does no WCS
arithmetic of its own.

.. code-block:: python

	from astropy.io import fits

	header = fits.Header()
	header["NAXIS"], header["NAXIS1"], header["NAXIS2"] = 2, 2048, 2048
	header["CTYPE1"], header["CTYPE2"] = "RA---TAN", "DEC--TAN"
	header["CRPIX1"] = header["CRPIX2"] = 1024.5
	header["CRVAL1"], header["CRVAL2"] = 45.0, 12.0
	header["CDELT1"], header["CDELT2"] = -0.0005, 0.0005
	header["RADESYS"] = "ICRS"

	on_image = index.image_query(header)

	print(len(on_image), "sources on the image")
	print(np.round(on_image, 4))

.. code-block:: text

	2 sources on the image
	[[44.6671 11.7598]
	 [45.0941 12.2015]]

With a real file it is two lines::

	with fits.open("image.fits") as hdulist:
		matches = index.image_query(hdulist[0].header)

How the query runs
------------------

Every region query in this class has the same two steps, and this one is no exception: a
**coarse step** over ipix ranges that is guaranteed to return a superset, then an
**exact test** on each candidate. Here the coarse step is a cone over the region's
bounding circle and the exact test is AST's own ``pointInRegion``, so the answer is
exactly the set of indexed points AST considers inside the region.

Two region types describe a shape qlsc already defines identically, and for those the
whole query is handed to the built-in method instead — the exact test then runs as one
batched call into the C layer rather than one call into AST per candidate:

======================  =========================================================
region                  query
======================  =========================================================
``ASTCircle``           :py:func:`~qlsc.QLSCIndex.radial_query`, the same geodesic disc
``ASTPolygon``          :py:func:`~qlsc.QLSCIndex.polygon_query`, the same great-circle edges
======================  =========================================================

An ``ASTBox`` is **not** one of them. Its sides are parallel to the axes of its frame,
so on the sky they are parallels of declination rather than great circles, and the
polygon through its four corners is neither inside it nor around it. A box takes the
general path, where AST decides membership.

Which path ran never changes the answer, and ``fast_path=False`` says so:

.. code-block:: python

	fast = index.region_query(footprint)
	generic = index.region_query(footprint, fast_path=False)

	print("same answer:", np.array_equal(fast, generic))

.. code-block:: text

	same answer: True

What is refused, and why
------------------------

A **negated region** — everything *outside* a circle, say — is refused rather than
searched. AST reports the bounding circle of a negated region as the circle of the
un-negated one, which contains none of the negated region's area, so using it as the
coarse step would return nothing at all and look like an honest answer. The same check
catches a compound region with a negated component, and a polygon wound the other way
round (AST fills a polygon to the *left* of its vertex path, so reversing the vertices
asks about the rest of the sky). Pass ``all_sky=True`` to answer those by testing every
point in the index.

qlsc stores bare numbers and has no idea what frame they are in; an AST region always
knows its own. The two are compared, and a mismatch is refused rather than silently
mixed — a header with ``RADESYS = 'FK5'`` searched against an index that was declared
ICRS raises. Say which frame the index is really in with ``frame="FK5"``, convert the
region first with ``regionWithMapping``, or take responsibility for the match with
``skip_frame_check=True``. ``frame=None`` means ICRS, never "do not check".

Going the other way
-------------------

:py:mod:`qlsc.cornish_interop` also runs the conversion in the other direction:
:py:func:`~qlsc.cornish_interop.ipix2region` turns one pixel into an ``ASTPolygon`` and
:py:func:`~qlsc.cornish_interop.ipix2moc` turns a set of pixels into an IVOA
Multi-Order Coverage map, which is what Aladin, MOCpy, VizieR and pgsphere read. See
:doc:`../examples/pixel_regions`.
