# Changelog

## Unreleased

### New

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
