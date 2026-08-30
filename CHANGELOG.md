# Changelog

## Unreleased

### New

* `QLSC.ipix2ang`, `QLSC.ipix2ang_center`, and `QLSC.ipix2xy` accept an array-like of
  ipix values (NumPy array, list, or tuple) in addition to a single value, following
  the same scalar/array convention as `ang2ipix`: a scalar (or 0-d array) returns
  scalars, anything with a shape returns arrays of that shape. The loop runs in the
  C extension. `ipix2polygon` remains scalar-only (see its docstring for why).
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
* `QLSCIndex.polygon_query(polygon)` and
  `QLSCIndex.ellipse_query(ra, dec, semi_major_axis, axis_ratio, position_angle)`
  join `radial_query` as region searches over the index. Both take the coarse ipix
  ranges from the embedded Q3C code (`q3c_poly_query`, `q3c_ellipse_query`, newly
  exported from the C extension) and then filter every candidate exactly, so the
  result is precisely the set of indexed points inside the region. `polygon_query`
  takes the polygon in the same form as `QLSC.point_in_polygon` and supports
  `return_key`; `ellipse_query` supports `return_key` and `return_distance` (the
  separation from the center of the ellipse), accepts `astropy.units.Quantity`
  angles, and with `axis_ratio=1` returns exactly what `radial_query` does for the
  same radius.

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
  (semi-minor over semi-major, in `(0,1]`), not the eccentricity.
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
  is coarser. A depth finer than the index's raises `ValueError`. `examples/sky_coverage_density_map.py`
  writes such a document from the WISE sample catalog.

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

### Breaking changes

* An empty `QLSCIndex` is now falsy. `QLSCIndex` defines `__len__`, so
  `if index:` asks whether the index holds any points, where it was previously
  always true. Use `if index is not None:` for the old meaning.

### Fixed

* `QLSC.ipix2ang` rejected NumPy integers (`isinstance(ipix, int)` is `False` for
  `np.int64`), so the ipix values returned by `ang2ipix` for array input could not be
  fed back to it. Integer values of any width are now accepted; floats and other
  non-integer types still raise `ValueError` rather than being truncated.
* `QLSC.ipix2xy` did not bounds check its ipix value and silently returned the
  coordinates of no pixel at all for out-of-range input; it now raises `ValueError`,
  as `ipix2ang` already did. The check is also enforced in the C extension.
* `QLSCIndex.radial_query` built its record array through `np.core.records`, which
  NumPy 1.25 deprecated and NumPy 2 removed; it now uses the public `np.rec`
  spelling, which has been available since long before the declared NumPy floor.
* Removed three build warnings: the C extension discarded the `const` qualifier on
  the `q3c_prm` bit tables in `ipix_to_xy`, and setuptools warned that the
  C source directories `qlsc/c_code` and `qlsc/c_code/q3c` were absent from the
  `packages` configuration. They are data files of the `qlsc` package, not
  importable packages, and are now declared as such; the sdist and wheel contain
  exactly the same files as before.

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
