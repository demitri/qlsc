# Changelog

## Unreleased

### Packaging

* The extension now compiles against an explicit NumPy C API level
  (`NPY_TARGET_VERSION`) rather than whichever level the NumPy doing the build
  happens to default to. NumPy's default tracks the headers that build the
  extension (NumPy 2.5 defaults to the 1.25 API level), and a compiled
  extension refuses to import under a NumPy older than the level it was built
  against - so a binary built with a current NumPy would have required
  numpy>=1.25 at runtime while the package declared `numpy>=1.18`.

## 1.1.0 (unreleased)

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
