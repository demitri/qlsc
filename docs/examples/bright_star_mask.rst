Bright Star Avoidance With a Real Catalog
=========================================

.. check-outputs-page: needs-network -- every block depends on the Bright Star Catalogue fetched from the VizieR TAP service

.. contents::
   :local:
   :depth: 2

**The problem.** A bright star ruins the sky around it. It saturates the
detector, bleeds along columns, scatters light into a halo tens of arcminutes
across, and leaves diffraction spikes over everything nearby. Every imaging
survey therefore carries a bright-star mask, and every target list has to be
checked against it: *is this position too close to a naked-eye star to be worth
observing?*

The check is a cone search per target, with a radius that depends on the star's
brightness. Twenty thousand targets is twenty thousand cone searches, which is
where an index stops being a convenience.

Unlike the other examples on this site, this one uses **real data**: the Yale
Bright Star Catalogue (Hoffleit & Jaschek), VizieR catalog ``V/50``, fetched over
the network. Nothing beyond the Python standard library and NumPy is needed to
get it.


Fetching the catalog
--------------------

VizieR exposes a TAP service that answers ADQL over plain HTTP, so a small
function with ``urllib`` is enough. It caches to a local file so that rerunning
the example does not hit the server again:

.. code-block:: python

	import csv
	import os
	import urllib.parse
	import urllib.request

	import numpy as np
	from qlsc import QLSCIndex

	TAP = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
	ADQL = '''SELECT HR, RAJ2000, DEJ2000, Vmag
	          FROM "V/50/catalog"
	          WHERE Vmag IS NOT NULL'''

	def bright_star_catalogue(cache="bsc5.csv"):
		'''Fetch the Yale Bright Star Catalogue (V/50) from VizieR, with a local cache.'''
		if not os.path.exists(cache):
			form = urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL",
			                               "FORMAT": "csv", "QUERY": ADQL}).encode()
			with urllib.request.urlopen(TAP, data=form, timeout=120) as response:
				with open(cache, "wb") as f:
					f.write(response.read())
		with open(cache, newline="") as f:
			rows = [r for r in csv.DictReader(f) if r["Vmag"]]
		hr   = np.array([r["HR"] for r in rows])
		ra   = np.array([float(r["RAJ2000"]) for r in rows])
		dec  = np.array([float(r["DEJ2000"]) for r in rows])
		vmag = np.array([float(r["Vmag"]) for r in rows])
		return hr, ra, dec, vmag

	hr, ra, dec, vmag = bright_star_catalogue()
	print(f"stars: {len(hr)}   V range: {vmag.min():.2f} to {vmag.max():.2f}")

.. code-block:: text

	stars: 9096   V range: -1.46 to 7.96

The brightest entry, HR 2491 at V = −1.46, is Sirius, and the catalog runs down to about
seventh magnitude — essentially everything visible to the unaided eye.

.. note::

	If you have no network access, generate a synthetic catalog instead: draw
	positions uniformly on the sphere as in :doc:`cone_search` and magnitudes from
	a power-law number count. The rest of this page works unchanged; only the
	specific stars named at the end will differ.


Building the index
------------------

The Harvard Revised number becomes each point's key, which is what lets us get
back to the magnitude after a search:

.. code-block:: python

	index = QLSCIndex()
	index.add_points(ra=ra, dec=dec, keys=[f"HR {h}" for h in hr])
	print("indexed:", len(index))

	magnitude = {f"HR {h}": v for h, v in zip(hr, vmag)}

.. code-block:: text

	indexed: 9082

Fourteen stars went missing, and that is not a bug. The index treats a coordinate
as the identity of a point, and this catalog contains fourteen pairs of entries —
components of visual doubles, with distinct HR numbers — listed at *identical*
rounded positions. For example HR 887 and HR 888 both sit at (44.80292,
+21.34028). The second insert of each pair is discarded, along with its key.

Whether that matters depends on what you are doing. For a bright-star mask it
does not: the pair is one star's worth of glare either way. If you needed both
components you would perturb the duplicate positions, or key the index on
something you can look up separately. It is the sort of thing worth checking
against ``len()`` every time you load a real catalog.


The mask radius
---------------

The halo of a star scales with its brightness. A real survey calibrates this
against its own imaging; here is a simple illustrative relation — 30 arcminutes
at V = 0, halving every two magnitudes:

.. code-block:: python

	def halo_radius(v):
		'''An illustrative masking radius in degrees for a star of magnitude v.'''
		return 30.0 * 2.0 ** (-v / 2.0) / 60.0

	MAX_RADIUS = halo_radius(0.0)
	print(f"largest halo: {MAX_RADIUS*60:.1f} arcmin (V=0);  V=6 -> {halo_radius(6)*60:.2f} arcmin")

.. code-block:: text

	largest halo: 30.0 arcmin (V=0);  V=6 -> 3.75 arcmin

Because the radius varies from star to star, the search cannot be a single query
per target at "the" radius. The pattern that works is: **search once at the
largest radius any star could need, then filter the matches with each star's own
radius.** The keys make the second step a dictionary lookup.


Screening the target list
-------------------------

.. code-block:: python

	rng = np.random.default_rng(seed=20200422)
	n_targets = 20_000
	t_ra  = rng.uniform(0., 360., n_targets)
	t_dec = np.degrees(np.arcsin(rng.uniform(-1., 1., n_targets)))

	rejected = 0
	worst = None
	for i, (cra, cdec) in enumerate(zip(t_ra, t_dec)):
		near = index.radial_query(ra=cra, dec=cdec, radius=MAX_RADIUS,
		                          return_key=True, return_distance=True)
		if len(near) == 0:
			continue
		v = np.array([magnitude[k] for k in near["key"]])
		inside = near["distance"] <= halo_radius(v)
		if inside.any():
			rejected += 1
			offender = near[inside][np.argmin(near[inside]["distance"] / halo_radius(v[inside]))]
			if worst is None or offender["distance"] < worst[1]["distance"]:
				worst = (i, offender)

	print(f"targets rejected: {rejected} of {n_targets} ({100*rejected/n_targets:.2f}%)")
	i, o = worst
	print("closest offending pair:")
	print(f"  target {i}: ra={t_ra[i]:.5f} dec={t_dec[i]:+.5f}")
	print(f"  star {o['key']} (V={magnitude[o['key']]:.2f}) at {o['distance']*60:.3f} arcmin, "
	      f"halo {halo_radius(magnitude[o['key']])*60:.3f} arcmin")

.. code-block:: text

	targets rejected: 97 of 20000 (0.48%)
	closest offending pair:
	  target 5405: ra=17.00736 dec=+32.01001
	  star HR 327 (V=6.25) at 0.165 arcmin, halo 3.439 arcmin

Half a percent of the sky lost to bright stars, from a catalog of nine thousand.
That is the number a survey planner wants, and the twenty thousand cone searches
that produced it took 1.49 seconds — about 75 microseconds each, including the
magnitude lookup and the filtering.

Note the comparison ``near["distance"] <= halo_radius(v)``: it is elementwise
across the whole record array, so no second loop is needed. This is the payoff of
``return_distance`` — the separations are already there.


Looking at one field
--------------------

The same query answers the "what is around here?" question directly. Pointing it
at the Pleiades:

.. code-block:: python

	plei = index.radial_query(ra=56.75, dec=24.1167, radius=2.0,
	                          return_key=True, return_distance=True)
	plei.sort(order="distance")
	print(f"{len(plei)} BSC stars within 2 deg of (56.75, +24.117):")
	for row in plei[:10]:
		print(f"  {row['key']:>8}  V={magnitude[row['key']]:5.2f}  {row['distance']*60:6.2f}'")

.. code-block:: text

	16 BSC stars within 2 deg of (56.75, +24.117):
	   HR 1165  V= 2.87    6.68'
	   HR 1156  V= 4.18   13.68'
	   HR 1149  V= 3.87   22.01'
	   HR 1152  V= 6.43   27.90'
	   HR 1142  V= 3.70   29.09'
	   HR 1178  V= 3.63   29.85'
	   HR 1180  V= 5.09   29.96'
	   HR 1151  V= 5.76   30.23'
	   HR 1140  V= 5.46   31.79'
	   HR 1145  V= 4.30   32.29'

Sixteen naked-eye stars inside two degrees, and the brightest of them, HR 1165 at
V = 2.87, is Alcyone. The cluster is a good test case precisely because it is the
kind of place a mask has to work: a dozen overlapping halos in one field.


See also
--------

* :doc:`cone_search` — the same query explained in detail, and benchmarked.
* :doc:`../user_guide/indexing` — what keys are and are not.
* `VizieR <https://vizier.cds.unistra.fr>`_ — the catalog service; the Bright
  Star Catalogue is `V/50 <https://vizier.cds.unistra.fr/viz-bin/VizieR-3?-source=V/50>`_.
