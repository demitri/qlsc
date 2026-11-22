# Changelog

## Unreleased

### New

* Lat/lon mode, which the README has promised since 2020: `GeoQLSC` and `GeoQLSCIndex`
  are the scheme and the index in geographic vocabulary. The mapping is small, because
  most of it is the identity — **latitude is declination** (same range, same sign), and
  **longitude is right ascension on the way in** (right ascension is periodic and every
  entry point already reduces it, so `lon=-74` may be given anywhere an ra can). The
  only real difference is on the way out: these classes report longitude in
  **[-180,180)** where `QLSC` reports right ascension in [0,360), and they do so from
  *every* method that returns one, the inherited `ipix2ang`/`xy2ang`/`ipix2polygon`
  included. Both classes take the pair under either spelling (`lon`,`lat` or `ra`,`dec`,
  never both at once) and add `lonlat2ipix`/`ipix2lonlat`/`ipix2lonlat_center`.
  `GeoQLSCIndex` stores longitudes already wrapped, so a point added as `lon=286` and
  one added as `lon=-74` are the same point, and what a query returns is what was put
  in; a latitude outside [-90,90] raises rather than being read as a point on the far
  side of a pole. Nothing about the pixelation, the ipix numbering, or the C code
  changes. The public `qlsc.wrap_longitude()` performs the [0,360) → [-180,180) mapping
  on its own. There is one documented exception to the convention: `coverage()` writes
  its vertices in [0,360), because a sky-coverage document declares the ICRS frame and
  is read by other software, so its numbers follow that frame rather than the producing
  object's reporting preference — the document a `GeoQLSC` emits is byte-for-byte the
  one a `QLSC` emits.
* First-class astropy support. Every method that takes an ra,dec pair now also accepts
  an `astropy.coordinates.SkyCoord` (scalar or array) in place of the pair — passed as
  the *first* of the two parameters with the second left unset, e.g. `ang2ipix(coord)`,
  `index.radial_query(coord, radius=1*u.deg)`, `point_in_ellipse(coord, ra_cen=center)`
  — as well as through `points=` on `ang2ipix`/`add_points` and as the `polygon=`
  vertices of `point_in_polygon`/`polygon_query`. A bare coordinate frame object works
  the same way. **The coordinate is transformed to ICRS**, the frame every ra,dec value
  in this package is taken to be in: a galactic or FK5 coordinate is converted, not
  reinterpreted. Giving both a SkyCoord and the second angle raises rather than
  silently discarding one of them. Individual angles may also be given as
  `astropy.units.Quantity` values, as radii already could. A coordinate carrying a
  distance keeps only its direction; the index is a sphere. An **array** `SkyCoord`
  given to `radial_query` is a target list, so it selects the batch form exactly as
  arrays of degrees do.
* `return_skycoord=True` on `QLSC.ipix2ang`, `ipix2ang_center`, `xy2ang`,
  `ipix2polygon`, and on `QLSCIndex.radial_query`, `polygon_query`, `ellipse_query`,
  `find_by_key` returns an ICRS `SkyCoord` instead of degrees (scalar in, scalar out;
  query matches always come back as an array SkyCoord, never squeezed). Combined with
  any of `return_key`, `return_distance` or a batch query's `target_index` a query
  returns an `astropy.table.QTable` with a `coord` column holding the SkyCoord plus the
  requested columns, since a SkyCoord cannot be a field of a `numpy.recarray`; its
  columns are in the recarray's order, `target_index` first. No default return value
  changes, and asking for a SkyCoord without astropy installed raises a clear
  `ImportError` — astropy remains an optional dependency.
* `QLSC.ipix2ang`, `QLSC.ipix2ang_center`, and `QLSC.ipix2xy` accept an array-like of
  ipix values (NumPy array, list, or tuple) in addition to a single value, following
  the same scalar/array convention as `ang2ipix`: a scalar (or 0-d array) returns
  scalars, anything with a shape returns arrays of that shape. The loop runs in the
  C extension. Integer values of any width are accepted, unsigned ones included;
  non-integer input raises rather than being truncated to a different pixel, at the
  `qlsc.q3c` level as well as through `QLSC` (a bare "give me an int64 array"
  conversion refuses an existing float array but will happily *build* one from
  `[1.5, 2.5]`, so the dtype of the input is settled before anything is converted).
  `ipix2polygon` remains scalar-only (see its docstring for why).
* `QLSCIndex.radial_query(..., return_distance=True)` returns the angular separation
  of each match from the query center in degrees. The cone filter already computes
  sin²(d/2) for every candidate, so the value costs nothing. With either
  `return_distance` or `return_key` set the result is a `numpy.recarray` with the
  fields `ra,dec` plus `key` and/or `distance`; with both flags off the return value
  is unchanged (an array of shape (n,2)). Points matched as equivalent
  representations of the query center (ra=360 vs ra=0, any ra at a pole) report a
  separation of exactly 0 rather than a floating point residual.
* `QLSCIndex` gained a dict-like interface over its keys: `len(index)`,
  `key in index`, `for key in index` / `index.keys()` (a generator), and
  `index.remove_point(key=...)` or `index.remove_point(ra=..., dec=...)`, which
  returns the number of points removed and raises `KeyError` if there were none.
  Keys are not required to be unique, so removing by key removes every point
  carrying it, and `keys()` yields a shared key once per point; points added
  without a key are not represented in `keys()` at all (they can only be removed
  by coordinate). A non-string key raises `TypeError` rather than silently missing:
  the underlying column has SQLite's TEXT affinity.
* `QLSCIndex.find_by_key(key)` returns the coordinates of every point carrying a
  key, as a `numpy.recarray` with the fields `ra,dec` (always an array, even for a
  single match). An absent key raises `KeyError`, as `remove_point` does; use
  `key in index` for a check that does not raise.
* The `QLSCIndex` class docstring now states the contract of the mapping-style
  operations: the index is a collection of labelled points — a multimap — rather
  than a dictionary, so `len()` counts points, iteration yields one key per point
  (repeating a shared key, omitting keyless points), and lookup and removal by key
  act on every point carrying it.
* `QLSCIndex.polygon_query(polygon)` and
  `QLSCIndex.ellipse_query(ra, dec, semi_major_axis, axis_ratio, position_angle)`
  join `radial_query` as region searches over the index. Both take the coarse ipix
  ranges from the embedded Q3C code (`q3c_poly_query`, `q3c_ellipse_query`, newly
  exported from the C extension) and then filter every candidate exactly, so the
  result is precisely the set of indexed points inside the region. `polygon_query`
  takes the polygon in the same form as `QLSC.point_in_polygon` and supports
  `return_key`; `ellipse_query` supports `return_key` and `return_distance` (the
  separation from the center of the ellipse) and accepts `astropy.units.Quantity`
  angles. Both always return an `(n,2)` array, for one match as well as for many —
  they do not inherit `radial_query`'s squeezed `(2,)` single-match shape, which is
  kept for compatibility and queued for a 2.0 release.

  An `axis_ratio` of exactly 1 is a circle, and `ellipse_query` answers it by running
  `radial_query`, so the two agree by construction. This matters because Q3C's ellipse
  is a conic in the gnomonic (tangent plane) projection about its center — semi-axes
  `tan(a)` and `tan(a) * axis_ratio` — and only the near hemisphere projects. A
  genuinely elliptical region is therefore limited to a semi-major axis below 90
  degrees and raises `ValueError` above it, rather than answering about a region that
  is not the one asked for; a circle has no such limit and takes the full `[0,180]`
  range. Along the major axis the region reaches exactly `semi_major_axis` degrees; the
  semi-minor axis is `arctan(tan(a) * axis_ratio)`, close to `a * axis_ratio` for a
  small ellipse and noticeably wider for a large one.

  Concave polygons are supported. Q3C's coarse coverage step over-covers them —
  it calls a pixel fully covered when all four of its corners are inside the
  polygon, which a notch can escape (upstream issue #55) — but over-coverage only
  enlarges the candidate set, and the exact filter removes the extras. The polygon
  must be small enough to project onto a single cube face, roughly 23 degrees
  across; anything larger raises `ValueError` rather than returning a partial answer.
* `QLSC.point_in_ellipse(ra, dec, ra_cen, dec_cen, semi_major_axis, axis_ratio, position_angle)`
  tests whether a point lies inside an ellipse on the sphere, as a sibling of
  `point_in_polygon`. It is the exact test `ellipse_query` filters with, and wraps
  Q3C's own `q3c_in_ellipse`. Like Q3C's SQL interface it takes the axis ratio
  (semi-minor over semi-major, in `(0,1]`), not the eccentricity. The same 90 degree
  limit and the same circular exemption apply as for `ellipse_query`; a point exactly
  on the edge of an elliptical region counts as outside (Q3C's convention), while one
  exactly on the edge of a circle counts as inside (the cone search convention).
* `QLSC.ipix_neighbors(ipix)` returns the pixels touching a given pixel, as a sorted
  array of unique ipix values — the building block for region growing and
  neighborhood searches. Eight neighbors in the interior of a cube face; **seven**
  for the 24 pixels on a corner of the cube, where three faces meet and there is no
  eighth pixel across the vertex; four at `depth=0`, where a face is a single pixel.
* `QLSC.coverage(ipix, values=None, ...)` and `QLSCIndex.coverage(depth=None, ...)`
  export sky coverage as a "sky-coverage" version 1 document — the JSON interchange
  format specified by the [astplot-js](https://github.com/demitri/astplot-js) project,
  whose `CoverageLayer` renders it as a filled-polygon coverage map or choropleth in
  any projection. Both return a plain Python dictionary; writing it out is one
  `json.dumps` with the standard library encoder, no custom `default=` hook (NumPy
  input does not leak `np.int64`/`np.float64` into the document). Each pixel becomes
  one polygon carrying its four corners as `vertices`, its ipix number as `id`, and
  optionally a scalar as `value`. The document declares degrees, the ICRS frame and
  geodesic edges, and records the qlsc version, the scheme and the depth in `meta`;
  no `equinox` is written, since ICRS has none and the format makes declaring one an
  error. `QLSCIndex.coverage` counts the index's points per pixel with a single
  `GROUP BY` over the ipix column — one query, no per-point work in Python — and
  defaults to depth 5 (6,144 pixels for the whole sky), or the index's depth if that
  is coarser. A depth finer than the index's raises `ValueError`. A per-pixel value that
  `float()` would convert but that is not a number — a numeric string, a boolean — is
  refused rather than coerced, matching qlsc's strictness about ipix values.
  `examples/sky_coverage_density_map.py` writes such a document from the WISE sample
  catalog.
* `QLSCIndex.radial_query(ra=array, dec=array, radius=...)` runs a cone search for a
  whole list of targets — matching a target list against the catalog — with `radius`
  either one value for every target or one value per target. The result is a
  `numpy.recarray` whose first field is `target_index`: the position in the target list
  each match belongs to. A target that matches nothing contributes no rows;
  `return_counts=True` returns `(matches, counts)` instead, where `counts` has one
  `np.int64` entry per target and a zero for a target that found none. A scalar query
  is unchanged, including the rule that a size-1 array is a scalar value; a target list
  with exactly one entry is therefore expressed with `return_counts=True`. The batch
  form is about the bookkeeping rather than speed: it runs one query per target (the
  ipix ranges covering a cone depend on where the cone is), and measured about 9%
  faster than the equivalent Python loop. `ra` and `dec` now accept an
  `astropy.units.Quantity` as `radius` already did, scalar or array: `np.asarray` on a
  Quantity discards its units, so a target list in radians would otherwise have been
  read as degrees and quietly returned the wrong points.
* `qlsc.crossmatch(index_a, index_b, radius, nearest_only=False, return_keys=True)`
  spatially joins two catalogs locally — the Q3C headline feature, without PostgreSQL.
  It returns a `numpy.recarray` of the pairs within the radius, with the fields
  `ra_a, dec_a, key_a, ra_b, dec_b, key_b, distance` (the key fields unless
  `return_keys=False`; a point stored without a key reports `None`). `nearest_only`
  keeps only the closest counterpart of each `index_a` point. The pairs are grouped by
  the `index_a` point in ascending `(ra_a, dec_a)` order and sorted by separation within
  a group, so the ordering is determined by the data alone. The search iterates over
  whichever index holds fewer points and runs the same z-order range search
  `radial_query` does against the other, filtering every candidate exactly with
  `sindist`; the two indices do **not** have to share a depth. Both in-memory and
  file-backed indices work, in any combination. Measured at about 15,000 points of the
  smaller catalog per second on in-memory indices at depth 30 (50,000 against 200,000
  points at one arcsecond, 22,593 pairs, in 3.3 seconds).
* `qlsc.cornish_interop` queries an index with the regions of
  [cornish](https://github.com/demitri/cornish), the Pythonic interface to Starlink
  AST. `region_query(index, region)` takes any bounded AST region — circle, box,
  polygon, MOC, or a compound (boolean) region built from them — and returns the
  indexed points inside it, using the same two steps the built-in region queries use:
  the ipix ranges covering the region's bounding circle as a coarse superset, then
  AST's own `pointInRegion` as the exact filter, so membership is decided by AST and
  the index only narrows the candidates. `ipix2region(qlsc, ipix)` returns a pixel as
  an `ASTPolygon` (the two libraries agree exactly on what a spherical polygon is, but
  wind their vertices oppositely, which this handles), and `ipix2moc(qlsc, ipix_list)`
  returns a set of pixels as an IVOA MOC — the standard coverage object Aladin, MOCpy,
  VizieR and pgsphere read.

  `cornish` is an optional dependency, guarded the same way `astropy` and `pandas` are:
  qlsc neither requires nor installs it, the module imports without it, and calling
  anything in it then raises an `ImportError` naming the package and chaining the
  original import error.

* `QLSCIndex.region_query(region)` is that bridge as a method on the index, and it
  recognizes the shapes qlsc already defines: an `ASTCircle` is a geodesic disc, so the
  query is handed to `radial_query`, and an `ASTPolygon` joins its vertices with great
  circles exactly as Q3C does, so it is handed to `polygon_query`. The exact test then
  runs as one batched call into the C layer instead of one `pointInRegion` call per
  candidate. The answer is the same set of points either way; `fast_path=False` forces
  the general path. An `ASTBox` deliberately keeps the general path: its sides are
  parallel to the axes of its frame, so on the sky they are parallels of declination
  rather than great circles, and the polygon through its corners is neither inside it
  nor around it. A polygon takes its fast path only under `qlsc.q3c.Q3C_MAX_N_POLY_VERTEX`
  vertices and inside `cornish_interop.POLYGON_FAST_PATH_MAXIMUM_RADIUS`; a larger one
  is not refused, it takes the general path.
* `QLSCIndex.image_query(fits_header)` answers "which of my sources land on this
  image?" from a FITS header alone: AST reads the world coordinate system and returns
  the sky region bounding the field, the index prefilters by ipix ranges and tests each
  candidate exactly. qlsc does no WCS arithmetic of its own. The region carries the
  frame the header declares, so an `FK5` header searched against an index declared ICRS
  raises rather than mixing the two. One upstream caveat is documented on the method: a
  field 180 degrees across segmentation-faults inside AST's own header handling (a
  cornish/starlink-pyast defect, reproducible without qlsc), which nothing on this side
  can catch. Ordinary images are nowhere near it.
* `qlsc.q3c.Q3C_MAX_N_POLY_VERTEX` is published by the extension module: the largest
  number of vertices Q3C's polygon routines accept (100). It was already enforced in C;
  now the Python layer can read it rather than keeping a second copy of the number.

* `QLSC.ipix2xiyi()` and `QLSC.xiyi2ipix()` expose the integer grid the scheme is
  built from: an ipix value is the bit interleaving of a face number and a pair of
  integers in `[0,nside)`, and these two convert between them in both directions,
  exactly. The decode half already existed inside the C extension but only as an
  intermediate step of `ipix2xy`, which threw the integers away in favour of the
  rescaled floats; the encode half is the vendored `q3c_xiyi2ipix`, previously
  unbridged. The round trip reproduces every pixel of every depth, with no
  trigonometry and no tolerance, which makes stepping the grid — the natural way to
  walk from a pixel to its neighbour — exact where nudging `x`,`y` and going back
  around through `(ra,dec)` is not. A grid coordinate outside `[0,nside)` is refused
  rather than reduced modulo `nside`: it is a position on a *neighbouring* face,
  whose axes are rotated with respect to this one. Both are also available at the
  extension level as `qlsc.q3c.ipix2xiyi` / `qlsc.q3c.xiyi2ipix`.
* `QLSC.nearby_ipix_ranges()` (and the `lon`/`lat` spelling on `GeoQLSC`) exposes
  Q3C's cheap neighbourhood calculation, `q3c_get_nearby` — the one its PostgreSQL
  extension uses to seed a positional join, previously unbridged. It skips the
  quadtree refinement of a cone search and returns at most four ranges, each a single
  quadtree square just coarse enough to hold the cone's bounding box on a cube face:
  a first cheap cut at "which pixels should I even look at?", not an answer. Its
  contract is documented as measured rather than as hoped — Q3C emits one square per
  neighbouring cube face and nothing in its code bounds the rest of that face's
  bounding box to that square, so the ranges are **not** a guaranteed cover, and both
  the docstring and `test_nearby_ranges.py` say so plainly. A negative radius raises
  instead of returning the empty result Q3C answers it with, and the whole-sky
  sentinel Q3C uses past `Q3C_MAXRAD` is reported as the pixels it means. Also
  available as `qlsc.q3c.nearby_query`.

### Documentation

* The documentation is now a structured site rather than an API dump: a landing
  page, a "Getting Started" tutorial, a four-part user guide (the scheme, depth and
  resolution, coordinates and conversions, indexing and searching), and six worked
  examples with science framing. Every code block in the new pages was run and its
  real output pasted; the catalogs are seeded so the numbers reproduce.
* The worked examples cover a transient cone search (with a measured comparison
  against a brute-force NumPy scan), bright-star avoidance against the Yale Bright
  Star Catalogue fetched from VizieR, survey-footprint membership with
  `polygon_query` including concave footprints, elliptical apertures on an inclined
  galaxy disc, pixel-level region growing with `ipix_neighbors`, and sky density
  maps built from pixel counts.
* Documented, with measurements, that the pixels are **not** equal-area: Q3C uses the
  tangent-plane cube rather than the area-preserving Chan & O'Neill transform, so a
  pixel's solid angle is `dx dy (1 + x² + y²)^(-3/2)` and varies by up to a factor of
  3^(3/2) ≈ 5.196 between the centre of a cube face and a cube vertex. Region queries
  are unaffected — they filter candidates exactly — but converting pixel counts to
  surface densities is not: use `QLSC.ipix_area` for a pixel's own area and
  `4*pi/nbins` only where the mean is genuinely what is wanted. (The measurements
  that established this are what prompted the `ipix_area` fix and the README
  correction listed above; the documentation and the README now quote the same
  numbers and cross-reference rather than duplicate.)
* Documented that `QLSC.ang2ipix_xy` returns face coordinates on **half** the scale
  used by `ipix2xy`, `xy2ang` and `xy2ipix` (−0.5 to +0.5 across a face rather than
  −1 to +1), and that `qlsc.sindist` returns the haversine sin²(d/2), not the sine
  of the separation — the finding that prompted the docstring correction listed
  above.
* `docs/generate_figures.py` generates every figure in `docs/_static/figures` from
  seeded synthetic catalogs; none are hand-made. It needs matplotlib and astropy,
  neither of which is required to build the documentation.
* `docs/check_example_outputs.py` re-runs the documentation's example code and diffs
  it against the output the pages show, so the claim that every output block was
  really produced by the code above it is enforced rather than merely asserted. It
  covers 77 code/output pairs across twelve pages. Blocks that cannot run
  hermetically are marked in the page source with the reason (`run-only` for
  wall-clock timings, `setup` for a block whose output is not shown, and the
  page-level `needs-network` for the VizieR fetch, which runs only under
  `--network`). Its docstring states its blind spots — unpaired code blocks, numbers
  quoted in prose, and the figures, all of which stay a manual duty. Running it for
  the first time found seventeen output blocks whose leading whitespace had been
  lost in transcription; those are now the real output. `docs/README.md` documents
  how to build, check, and regenerate the documentation.
* The API reference covers this release's additions: `crossmatch`, `QLSC.coverage`
  and `QLSCIndex.coverage`, the batch form of `radial_query` (`target_index`,
  `return_counts`), `find_by_key`, `return_skycoord` and SkyCoord acceptance,
  `GeoQLSC`/`GeoQLSCIndex`/`wrap_longitude`, `qlsc.cornish_interop` (with a note
  that it needs the optional `cornish` package), and
  `qlsc.q3c.ELLIPSE_MAXIMUM_SEMI_MAJOR_AXIS`. The user guide covers batch queries,
  key lookup, astropy interoperability and the sky-coverage export with checked
  examples; dedicated worked-example pages for cross-matching, coverage and the
  geographic classes are still to come and are flagged as such in the examples index.
* Added `.readthedocs.yaml`, which Read the Docs now requires; without it the build
  has no configuration file to read.
* Fixed a malformed inline literal in the `QLSCIndex` docstring (``` ``file://` ```)
  that produced a docutils warning on every documentation build.

### Development

* `scripts/sync_q3c.py` mechanizes the vendored-Q3C upstream sync documented in
  `CLAUDE.md`: it fetches a q3c tag, checks that the vendored files are the previous
  release plus exactly the two documented local patches, decides for each patch
  whether the new release already contains the fix, re-applies the ones still needed
  by exact match (never fuzzed — an anchor that no longer matches stops the run),
  diffs every `common.h` declaration, struct and macro the C wrapper depends on,
  updates the `Q3C_VERSION` string, and runs the guard tests and the full suite.
  `--dry-run` reports without writing. This is a repository tool; it is not part of
  the installed package.
  Since the v2.0.6 sync it carries no local patches; its patch machinery is tested
  against a synthetic patch so that it stays ready for the next one.

### Breaking changes

* Omitting a required argument of `QLSCIndex.radial_query` or `QLSC.ang2ipix_at_depth`
  now raises `ValueError` where it used to raise `TypeError`. Their parameters were
  made formally optional so that a single `SkyCoord` can be given as the one positional
  coordinate argument (`radial_query(coord, radius=...)`), and Python can no longer
  report the missing argument itself; both methods check and name it instead. The
  arguments are exactly as required as they were, and the message is more specific —
  only the exception type differs, which matters solely to code catching `TypeError`
  around these two calls. `QLSC.face_number` already took its coordinate as optional
  and is unchanged.
* `QLSC.ipix_area`'s `depth` parameter now means something, and is validated. It was
  previously ignored along with `ipix`, so any value was accepted; it is now an
  absolute depth in `[0, qlsc.depth]` naming the ancestor pixel to measure (`depth=0`
  is a whole cube face), and it defaults to the scheme's own depth. A depth finer than
  the scheme's raises `ValueError`: the underlying Q3C function takes a levels-up
  offset rather than a depth, and past the scheme's own resolution it shifts the
  pixel's coordinates out of the cube face and returns an area larger than a face
  without complaint. Code passing a depth this scheme cannot resolve — for example
  `QLSC(depth=5).ipix_area(ipix, 30)` — was already getting a meaningless answer and
  must now pass a depth in range.
* An empty `QLSCIndex` is now falsy. `QLSCIndex` defines `__len__`, so
  `if index:` asks whether the index holds any points, where it was previously
  always true. Use `if index is not None:` for the old meaning.
* `QLSCIndex.add_point` and `QLSCIndex.add_points` require keys to be strings
  (or `None`, meaning "no key"), where any object was previously accepted and
  stored as its text representation. `add_point(..., key=5)` stored the string
  `'5'`, after which `5 in index` raised and only `'5' in index` found the point —
  the key that went in was not the key that came back out. Code passing integer
  or other non-string keys must convert them: `key=str(catalog_id)`. Existing
  index files are unaffected; their keys were always stored as text.

### Fixed

* The embedded Q3C code is updated from v2.0.5 to v2.0.6. Upstream fixed both bugs
  qlsc had been patching locally since 1.1.0 — the undersized buffer in `q3c_poly.c`
  reachable through `point_in_polygon` (q3c issue #50) and the integer overflow in
  the quadtree decomposition of tiny query regions (#51) — so both patches are
  dropped and the embedded files are now byte-for-byte upstream's. Pixel numbering
  is unchanged: `ang2ipix` and everything that encodes or decodes an ipix give the
  same answers, and existing index files remain valid. What changes:
  * `point_in_ellipse`, and so `ellipse_query`, now give the right answer for small
    ellipses. Upstream rewrote the containment test (the old one, derived with
    Maple, subtracted nearly equal numbers), and the old one was wrong more often
    than its release note suggests: over 4,000 random ellipses with semi-major axes
    from 1e-7 to 80 degrees and 200 points each, placed within twice the semi-major
    axis of the centre, 38,364 of the 800,000 answers changed, and on every one of
    those, checked against the tangent plane conic evaluated to 50 digits, the new
    answer was right and the old one wrong. Nearly all of the changes are for ellipses below
    1e-5 degrees (36 milliarcseconds), where the old test misclassified points well
    inside and outside the edge, not just on it; the largest ellipse with a changed
    answer had a semi-major axis of 1.9e-4 degrees (0.7 arcseconds), and none of
    the 2,193 sampled ellipses above 1e-3 degrees changed any answer. The
    measurement is `scripts/q3c_206_ellipse_differential.py`.
  * `nearby_ipix_ranges` closes the gap its warning described: where a cone spills
    onto two or three neighbouring cube faces, the one square Q3C emits for each
    neighbouring face is now grown until it contains that face's whole bounding
    box, rather than sitting at one corner of it. Upstream found thin ellipses
    crossing a cube corner that the old square missed. The ranges can be coarser
    in those cases; the method is still documented as a measured, not proven, cover.
  * Region queries near cube corners that produced more ipix ranges than Q3C's
    output arrays hold are now retried by Q3C itself with fewer, coarser ranges.
    qlsc's own handling of that case (`RangeOverflowError`, and the index widening
    its search region) remains as a guard that should no longer fire.
* `QLSC.xy2ipix()` returns the pixel containing the position it is given, at every
  depth. It computed its answer by projecting the position back onto the sphere with
  `xy2ang` and handing that to `ang2ipix` — two trips through trigonometry, whose last
  bits decided the answer for a position sitting on a pixel boundary. The corners
  `ipix2xy` returns are exactly such positions, so the two methods were not inverses:
  at `depth=30`, `xy2ipix(ipix2xy(ipix))` came back with a *neighbouring* pixel for
  1206 of a 3000-pixel sample, 40% of the scheme. It now scales the position to the
  face's integer grid and truncates — the same arithmetic Q3C applies once it has
  projected a sky position — and converts that with the new `xiyi2ipix`. Positions in
  the interior of a pixel were never in doubt and are unchanged; a test pins that they
  still agree with the route through the sky.
* `GeoQLSC.ang2ipix`, `lonlat2ipix`, `point_in_polygon` and `point_in_ellipse` refuse a
  latitude outside [-90,90] instead of quietly answering for somewhere else.
  `GeoQLSCIndex`, `face_number` and `ang2ipix_xy` already refused one; these delegated to
  the astronomical normalization, so `lonlat2ipix(lon=10, lat=100)` returned the pixel of
  the point at (190, 80) — a real place on the far side of the Earth — with no error. The
  check now lives at the single seam every geographic entry point crosses, and
  `test_geo_latitude_gate.py` derives the list of methods to check from their signatures,
  so a geographic method added later that takes a latitude fails the gate until it is
  covered. `QLSC` is unchanged: an unnormalized declination is a legitimate celestial
  coordinate and is still read as the point past the pole.
* `QLSC.coverage` and `QLSCIndex.coverage` validate a caller's `meta` instead of copying
  it in unexamined, restoring the guarantee the module advertises: a dictionary
  `json.dumps` serializes with the standard library encoder alone, `allow_nan=False`
  included. The values a producer naturally records there were exactly the ones that
  broke it — a `np.int64` count out of `np.unique` is not serializable at all, and a NaN
  is not valid JSON — and both failed at `json.dumps`, a call after the document was
  handed back, naming a type rather than the entry that carried it. NumPy scalars and
  arrays are now converted (their Python equivalents are exact), non-finite numbers,
  non-string keys, circular structures and values with no JSON form raise, and the
  message names the offending entry. The caller's own dictionary is left untouched.
* The documentation's two density figures showed densities about 2.55 times too high.
  Both divided every pixel's count by `ipix_area(0, depth)` — pixel 0 sits at a cube
  vertex, which makes it the *smallest* pixel in the scheme, not a representative one —
  while the page beside them said the mean area had been used. The pixel-scale figure had
  the same fault, plotting the smallest pixel's side against a page that tabulates the
  mean. The density figures now divide each count by that pixel's own area (a true
  density), the scale figure plots the mean side, the prose agrees with both, and
  `test_figure_pixel_areas.py` fails if a figure takes one pixel's area as the area of a
  pixel again.
* `QLSC.point_in_ellipse` now gives the same answer for a point however its coordinates
  are written. `ra=0` and `ra=360` are the same meridian and at a pole the right ascension
  is arbitrary, but neither of the two containment tests is exact on those: the circular
  path excluded the center of a zero-radius circle when the center was written the other
  way (`point_in_ellipse(360, 20, 0, 20, 0, axis_ratio=1)` was False where
  `radial_query` found the point, contradicting the docstring's claim that the two apply
  the identical test), and the elliptical path put a point at a pole a hair off the
  center. `QLSCIndex.ellipse_query` had the complementary defect: it matched any
  representation of its center outright, so it returned points its own documented filter
  called outside. Both now rewrite such a point as the center and run the ordinary test,
  which keeps Q3C's containment rules — including that a zero-size ellipse holds nothing
  at a pole — while making the answer independent of the representation.
* `qlsc.cornish_interop.ipix2moc` validates its pixel numbers instead of truncating them.
  It called `int()` on each value, so `ipix2moc(q, [1.9])` built a MOC covering pixel 1,
  `True` covered pixel 1, and the string `"3"` covered pixel 3 — a map of somewhere the
  caller never asked about, produced without a word. It now applies the same check every
  other entry point applies, over the whole iterable before any region is built.
* `QLSC.ang2ipix_at_depth` now reads its coordinate the way `ang2ipix` does. It called the
  C function directly, so a declination past the pole was truncated to the pole — a
  different pixel from the one `ang2ipix` names for the same coordinate — and a non-finite
  value reached the undefined float-to-int conversion of upstream q3c issue #52 that every
  other entry point guards against.
* Documented a trap that has always been there: **`QLSC.ang2ipix_xy` returns face
  coordinates on `[-0.5,0.5]`, while `ipix2xy`, `xy2ang` and `xy2ipix` use `[-1,1]`** —
  exactly a factor of two. Feeding `ang2ipix_xy`'s `x,y` straight into `xy2ang` returns
  the wrong coordinate and raises nothing, because a half-scale value is a perfectly
  legal full-scale one; it just names a different point. `ang2ipix_xy` now carries a
  warning with the doubling recipe, the other three point at it, and
  `test_face_xy_scale.py` pins the factor. The values themselves are unchanged: making
  the scales agree would break existing code with no error to warn it, so it is queued
  for a 2.0 release alongside the `radial_query` return-shape cleanup.
* `qlsc.sindist` (and `qlsc.q3c.sindist`) was documented as "the sine of the angular
  distance between two points on a sphere". It returns the **haversine**, `sin²(d/2)`,
  which is what Q3C computes and what every caller inside the package already assumed.
  The difference matters: `sin(d)` peaks at 90 degrees and falls back to zero at 180,
  so code that took the old description at face value and used the value as a distance
  ordering was wrong for every separation past a quadrant. Only the documentation
  changed; the returned values are unchanged. The docstring now also gives the
  conversion back to degrees, and `distance` now states that its result is in degrees.
* `QLSC.ipix_area` ignored both of its arguments and returned the mean pixel area,
  `4π / nbins`, with a `NOTE` in its docstring saying so and claiming the pixels "are
  very close" in size. They are not: the largest pixel of a cube face, at its centre,
  is `3√3 ≈ 5.196` times the smallest, at a cube vertex, and the spread is already a
  factor of 2.5 at `depth=2`. Against the mean the error runs from 1.91× too small to
  2.72× too large. The method now returns the area of the pixel actually asked for,
  from Q3C's own `q3c_pixarea`.
* `QLSC.ipix2face` returned the wrong face for the last pixel of every face at depths
  above about 26, and a `float` at every depth. It divided the ipix value by
  `math.pow(nside,2)`, and past 2^53 a double cannot hold the operands exactly: at
  depth 30 the very last pixel of the scheme was reported as being on **face 6**, which
  does not exist. The division is now integer arithmetic, and the return value is an
  `int` as the signature always said. Non-integer and out-of-range pixel numbers now
  raise `ValueError` instead of being truncated or clamped onto some other face.
* `QLSC.face_number(ipix=...)` raised `TypeError`. `ra` and `dec` were required
  positional parameters, so the documented pixel-number form could not be called at
  all; the only way in was the undocumented `face_number(None, None, ipix=n)`, and
  even that failed for `ipix=0` because the branch was guarded by `if ipix:`, which is
  false for the first pixel of face 0. `ra` and `dec` now default to `None`, every
  test is against `None` rather than truthiness, and the pixel form answers from
  `ipix2face` — exact integer arithmetic — rather than by converting the pixel to a
  coordinate and back.
* `qlsc.q3c.ellipse_query` accepted a semi-major axis of 90 degrees or more, where the
  companion predicate `qlsc.q3c.in_ellipse` refused it. Q3C's ellipse is a conic in the
  tangent plane at its center and means nothing past the near hemisphere, so a direct
  caller of the C function got coarse ranges for a region that was not the one asked
  for — safe today only because upstream returns the whole sphere for any region of
  `Q3C_MAXRAD` (35) degrees or more. Both functions now enforce the same limit. The
  `Python` layer was never affected: it validates first, and answers a circle of any
  radius with a cone search. `qlsc.q3c` now also exports `Q3C_MAXRAD` and
  `ELLIPSE_MAXIMUM_SEMI_MAJOR_AXIS`, and the whole-sphere fallback the safety argument
  rests on is pinned by `test_q3c_maxrad_fallback.py` so an upstream sync cannot remove
  it unnoticed.
* `QLSC.ipix2polygon(ipix, duplicate_endpoint=True)` did nothing. The branch built
  the closed five-row polygon and the line after it stripped the extra row back off
  unconditionally, so both forms returned the same four rows and the parameter was
  silently ignored. It now returns the closed form, shape `(5,2)`, with the first
  vertex repeated as the last. Either form may still be passed to `point_in_polygon`
  and `polygon_query`: a repeated vertex makes a zero-length edge, which the crossing
  test never counts.
* `QLSC.ipix2ang` rejected NumPy integers (`isinstance(ipix, int)` is `False` for
  `np.int64`), so the ipix values returned by `ang2ipix` for array input could not be
  fed back to it. Integer values of any width are now accepted.
* `QLSC.ipix2xy` did not bounds check its ipix value and silently returned the
  coordinates of no pixel at all for out-of-range input; it now raises `ValueError`,
  as `ipix2ang` already did. The check is also enforced in the C extension.
* `QLSC.face_number` returned the wrong cube face for any negative right ascension,
  and for some values a number that is not a face at all — `face_number(-315, 20)`
  returned `-2`. The vendored `q3c_get_facenum` is documented "strictly 0<=ra<360"
  and computes the face as `(char)fmod((ra + 45) / 90, 4) + 1`, which truncates
  toward zero for a negative argument; the C wrapper did not reduce ra first, so
  `face_number` and `ang2ipix_xy` disagreed about the same coordinate.
  `face_number` also no longer rejects ra outside [-360,360] — right ascension is
  periodic and is now accepted at any value, as `ang2ipix` accepts it. The range
  check had only hidden the fact that everything in [-360,0) was answered wrongly.
* `QLSC.point_in_polygon` and `QLSCIndex.polygon_query` chose the wrong cube face
  for a polygon whose **first vertex** had a wrapped right ascension — the natural
  way to write a polygon straddling ra=0. Both reach `q3c_get_facenum` through
  `q3c_get_facenum_poly`, which reads that vertex. A square at ra=-140 was rejected
  as "too large" while the identical square at ra=220 worked. The vertex values are
  now reduced at the C boundary, into a copy: the caller's array is not modified.
* Every entry point that hands a right ascension to the embedded Q3C code was
  audited for this; the audit is recorded as a block comment at the top of
  `qlsc_c_module.c` and pinned by `source/tests/test_ra_wrapping.py`, which asserts
  that every public API answers identically for `ra`, `ra+360k` and `ra-360k`.
  `qlsc.q3c.facenum` and the two polygon functions also now reject non-finite
  coordinates, as their siblings already did (upstream q3c issue #52).
* Lookups and removals by key scanned the whole index: the schema indexed `ipix`
  but not `key`, so `key in index`, `find_by_key` and `remove_point(key=...)` were
  proportional to the size of the catalog. New index files carry an index on `key`,
  and an existing file gains one when it is opened — the change is additive, so an
  older qlsc reading the same file simply ignores it and the schema version does not
  move. A file that cannot be written to (read-only, or locked) keeps the index it
  has and logs a warning saying so; the results are identical either way.
* `QLSCIndex.add_points` silently dropped points when given fewer keys than
  coordinates: the length check was the chained comparison
  `len(ra) != len(dec) != len(keys)`, which Python evaluates as
  `(len(ra) != len(dec)) and (len(dec) != len(keys))` — false whenever the
  coordinates agree, whatever the key count is. Three coordinates and two keys
  stored two points. It now raises `ValueError`, and the rows are built with a
  guard that refuses to combine sequences of unequal length, so the invariant is
  enforced where it is relied on and not only where the arguments are checked.
* `QLSCIndex.add_points` accepted `points` together with `ra`,`dec` and silently
  discarded the `ra`,`dec` pair; the check (`all([x is True for x in ...])`) could
  not be satisfied by an array. It now raises `ValueError`.
* `QLSCIndex.radial_query` built its record array through `np.core.records`, which
  NumPy 1.25 deprecated and NumPy 2 removed; it now uses the public `np.rec`
  spelling, which has been available since long before the declared NumPy floor.
* Removed three build warnings: the C extension discarded the `const` qualifier on
  the `q3c_prm` bit tables in `ipix_to_xy`, and setuptools warned that the
  C source directories `qlsc/c_code` and `qlsc/c_code/q3c` were absent from the
  `packages` configuration. They are data files of the `qlsc` package, not
  importable packages, and are now declared as such; the sdist and wheel contain
  exactly the same files as before.

### Packaging

* Releases now ship binary wheels, so installing no longer requires a working C
  compiler and the NumPy headers. They cover Linux (x86_64 and aarch64,
  manylinux and musllinux), macOS (arm64 and x86_64) and Windows (64-bit), for
  CPython 3.10 through 3.14. The wheels are built and tested by GitHub Actions
  (`.github/workflows/wheels.yml`, using cibuildwheel) when a version tag is
  pushed, and uploaded to PyPI from that run as part of the manual release; the
  full test suite is run against every wheel before it is published. The source
  distribution is built by the same workflow and remains the install path for
  every platform and interpreter outside that list.
* The extension now compiles against an explicit NumPy C API level
  (`NPY_TARGET_VERSION`) rather than whichever level the building NumPy happens
  to default to. Without it, a wheel built with a current NumPy would refuse to
  import under NumPy older than 1.25 while the package still declared
  `numpy>=1.18`; CI now installs the built wheel against the oldest NumPy
  available for the oldest supported interpreter and runs the tests.

## 1.1.0 (2026-08-30)

### Breaking changes

* **`QLSC.ipix2ang` now returns the *center* of the pixel** rather than the "lower left"
  corner, matching the contract change made in the Q3C PostgreSQL plugin v2.0.1.
  `ipix2ang_center` is now a synonym; pixel corners remain available via `ipix2xy`
  and `ipix2polygon`. No error is raised — the returned values are simply different —
  so code depending on corner semantics must switch to `ipix2xy`/`xy2ang`.
* `QLSCIndex.radial_query` now validates its radius to the range [0,180] degrees
  (`ValueError` otherwise), and the radius comparison is inclusive: a radius of 180°
  includes an exactly antipodal point, and a radius of 0 returns exactly coincident
  points (including equivalent coordinate representations, e.g. ra=360 vs ra=0).
* Operations on a closed `QLSCIndex` raise `RuntimeError` instead of silently
  operating on a fresh, empty database.
* The unused `qlsc.q3c.radial_query_it` function (never reachable through the
  Python API, and not thread-safe due to internal static caches) was removed
  from the C module.

### New

* `QLSC.point_in_polygon(ra, dec, polygon)`: spherical point-in-polygon tests
  (polygons up to 100 vertices, must fit within a single face projection).
* `QLSCIndex.close()` and context-manager support (`with QLSCIndex(...) as index:`).
* `qlsc.q3c.version()` returns the version of the embedded Q3C code.
* `qlsc.q3c.RangeOverflowError` (a `RuntimeError` subclass) is raised on internal
  Q3C range-capacity overflows instead of a plain `RuntimeError`.

### Fixed

* Embedded Q3C code updated from v2.0.0 to v2.0.5, picking up the upstream fix for
  cone searches with radii just below 35° that silently missed matches in some areas
  of the sky (and could crash), plus input-hardening fixes. The embedded copy now
  compiles via upstream's `Q3C_STANDALONE` (one local patch: an undersized buffer in
  `q3c_poly.c` reachable through `point_in_polygon`; see CLAUDE.md).
* The module-level `qlsc.distance`, `qlsc.sindist`, and `qlsc.xy2facenum` wrappers
  raised `NameError` when called; `xy2facenum` was additionally never exported by the
  C module. All three now work, and `xy2facenum` validates its face number.
* Radial queries on indices below depth 30 could silently return zero matches
  (degenerate Q3C ranges at small nside); queries now always search at depth 30 and
  map ranges down to the index depth.
* Radial queries with zero or very small radii could silently miss indexed points at
  some sky positions; the candidate search now widens its radius as needed.
* `QLSC.ang2ipix` now converts non-double input arrays as documented instead of
  raising `TypeError`.
* `QLSCIndex.add_points` stored un-normalized coordinates when any `dec` was outside
  [-90,90] (the normalized values were computed and discarded), and the `ra=`/`dec=`
  argument form raised a spurious `ValueError` in that case.
* Opening an existing index file created at a different depth than the provided
  `QLSC` raised `NameError`; opening a corrupt or locked index file now raises a
  clear error instead of being misreported as "not a QLSC index".
* `QLSCIndex.number_of_points` now uses `COUNT(*)` (previously `max(rowid)`, which
  returned `None` for an empty index).
* A per-call memory leak in the C radial-query wrapper (two leaked array references
  per call), a matching leak in array-form `ang2ipix`, and several connection leaks
  on error paths.
* `QLSC.ang2ipix` returned a bare `int` for length-1 array input (breaking
  `add_points` with a single point); array input now always returns an array.
* C-level functions now raise an exception instead of crashing the interpreter
  when handed an invalid `hprm` capsule, and `point_in_polygon` rejects non-finite
  coordinates instead of returning `False`.
* `ipix2xy` computed the face number in floating point, placing ipix values near
  face boundaries on the wrong face at depths ≥ 26.
* A signed integer overflow in the embedded Q3C quadtree decomposition for tiny
  query regions (radii below ~1e-7°; upstream q3c issue #51) produced undefined
  behavior and garbage coordinates; patched locally (UBSan-verified).
* Non-finite (NaN/Inf) coordinates are now rejected with `ValueError` at every
  public boundary (upstream q3c issue #52): NaN slipped past range checks like
  `abs(dec) > 90` into undefined float-to-int conversions, and `ang2ipix`
  silently returned `None` for NaN scalars.
* Removed broken vestigial result caches in the C `ang2ipix`/`ang2ipix_xy`
  wrappers: a missing `return` meant they never served cached values, but they
  leaked one Python integer per repeated-coordinate call.

## 1.0.6

Baseline release ([PyPI](https://pypi.org/project/qlsc/)).
