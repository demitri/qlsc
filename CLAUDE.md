# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

QLSC is a Python package (published on PyPI as `qlsc`) implementing the quadrilateralized spherical cube — a scheme for segmenting a sphere into approximately equal-area pixels, aimed at astronomical spatial indexing and cone searches. It embeds C code vendored from Sergey Koposov's Q3C PostgreSQL extension (copyright Koposov) and provides Q3C's catalog-indexing functionality without requiring PostgreSQL. Docs are on Read the Docs (https://qlsc.readthedocs.io); the README is the primary conceptual introduction to the scheme itself.

## Layout Quirk

The Python package does not live at the repo root: everything installable is under `source/` (`source/setup.py`, `source/qlsc/`, `source/tests/`). The repo root holds the README, docs, examples, figures (LaTeX/graphics sources for the README diagrams), and `c_test/` (an Xcode project for exercising the C code directly).

## Commands

```bash
# Build & install from source (compiles the C extension; NumPy headers required)
cd source && python setup.py install

# Run the test suite — see the gotcha below first
cd source/tests && python -m pytest

# Run a single test file / test
python -m pytest test_ang2ipix.py
python -m pytest test_distance.py -k test_distance

# Build the docs (Sphinx; deps in docs/requirements.txt, incl. sphinx_autodoc_typehints)
cd docs && make html
```

**Test gotcha:** the tests do `import qlsc`, and since the source tree contains no built `.so`, they run against the *installed* package, not the working tree. After changing anything under `source/qlsc/`, reinstall before running tests or you will test stale code. Test files prefixed `xtest_` (rather than `test_`) are deliberately disabled.

**Releasing:** the version lives only in `source/qlsc/version.py` (kept out of `__init__.py` so `setup.py` can read it without importing dependencies). The commit convention is a dedicated "Bumping version to X.Y.Z." commit. The procedure (established for v1.1.0): update `CHANGELOG.md`; `cd source && rm -rf dist && python -m build` (builds sdist + a wheel for this machine only; everyone else builds the sdist, which works under pip's default build isolation via `pyproject.toml`); `twine check --strict dist/qlsc-X.Y.Z*`; `twine upload --repository qlsc dist/qlsc-X.Y.Z*` (the `[qlsc]` section in `~/.pypirc`; twine ≥6 requires sections be listed under its `[distutils] index-servers` header); `git tag vX.Y.Z && git push origin master vX.Y.Z`; `gh release create vX.Y.Z` with the version's CHANGELOG section as notes. The PyPI project page renders the repo-root README (read by `setup.py` at build time, Markdown, absolute image URLs) — it is baked at upload and immutable per release; the v1.1.0 page predates this and shows plain text, which self-corrects at the next upload.

## Architecture

Three layers, from bottom up:

1. **Vendored Q3C C code** — `source/qlsc/c_code/q3c/` (`q3cube.c`, `q3c_poly.c`, headers). This is third-party code from the Q3C project; treat it as an embedded dependency — the four files are upstream Q3C v2.0.5, compiled without PostgreSQL via upstream's own `Q3C_STANDALONE` macro (defined in `setup.py`, which also carries the `Q3C_VERSION` string to keep in step), with exactly TWO local divergences, both for reachable bugs reported upstream (each guarded by a test in `source/tests/test_source_invariants.py`; on the next upstream sync, keep each patch unless the upstream release contains the fix):
   - `q3c_poly.c`: `q3c_check_sphere_point_in_poly` declares `points[8]` where upstream has `points[4]` — a stack buffer overflow for polygons whose projected bounding box crosses onto 3-4 extra faces (upstream issue #50), reachable via `QLSC.point_in_polygon`. The companion caller-side fix (5-row, not 3-row, projection arrays) lives in `qlsc_c_module.c`, which is qlsc's own code.
   - `q3cube.c`: the quadtree refinement depth is computed as `log2(nside/n0)` (three sites) where upstream uses the raw ratio — a signed integer overflow (`2 * 2^30`) for tiny query regions producing garbage coordinates (upstream issue #51), reachable via `q3c.radial_query` with radii below ~1e-7 deg.

   Known upstream issues deliberately NOT patched because their code is unreachable from qlsc (re-check on each sync): #53/#54 live in the un-vendored PostgreSQL layer (`q3c.c`); from #55, `q3c_poly_cover_check`'s concave-polygon over-coverage and `q3c_fast_get_xy_minmax`'s uninitialized `Q3C_POLYGON` branch sit on the unexported poly-query path, and `Q3C_BOX_INTERSECT` in `common.h` is broken but unused. Issue #52 (NaN/Inf reaching undefined float→int conversions) is fixed at qlsc's own boundaries (Python validation plus wrapper checks), not in the vendored files. Never edit these files for style or convenience, but a bug that actually impacts qlsc (reachable from the wrapper) DOES get patched locally and immediately — correctness outranks byte-for-byte purity. Document any such divergence here and in the commit, and reconcile it on the next upstream sync. A bug that is dormant in qlsc (unreachable code) is documented only, keeping the files byte-identical. To sync a newer upstream release: copy the four files from https://github.com/segasai/q3c verbatim, check whether any function signature used by `qlsc_c_module.c` changed, and re-apply the `points[8]` patch above unless upstream has fixed it. All the actual spherical-cube math (ipix encoding/decoding, radial query bin selection) lives here. The central runtime object is the `q3c_prm` struct ("hprm"), which holds precomputed bit-interleaving tables for a given `nside`. Note: as of qlsc v1.1, `QLSC.ipix2ang` returns the pixel *center*, matching upstream's contract change in Q3C v2.0.1 (`ipix2ang_center` is a synonym); pixel corners are available via `ipix2xy` and `ipix2polygon`.

2. **CPython extension** — `source/qlsc/c_code/qlsc_c_module.c`, compiled by `setup.py` into the extension module `qlsc.q3c`. It wraps the Q3C functions and passes the `q3c_prm` struct around as an opaque PyCapsule (`init_q3c()` creates it; every geometry call takes it back). Functions accept scalars or NumPy arrays (arrays must be `dtype=double`; other dtypes get copied).

3. **Python API** — `source/qlsc/qlsc_base.py` defines the two public classes:
   - **`QLSC(depth)`** — pure geometry: conversions among (ra,dec), ipix, and per-face (x,y) coordinates (`ang2ipix`, `ipix2ang_center`, `ipix2polygon`, `ipix_up`/`ipix_down`, …). `depth` ∈ [0,30]; each depth level subdivides every bin into four; depth 30 is the Q3C/PostgreSQL default. The face-plane (x,y) coordinate system (origin at face center, axes −1..+1) is the recommended working space for nontrivial calculations.
   - **`QLSCIndex(qlsc, filepath)`** — a cone-search index backed by SQLite (in-memory by default, or a file; table `qlsc_ipix`). Points go in via `add_point`/`add_points` (with optional string keys); `radial_query()` builds SQL from ipix ranges produced by the C-level `q3c.radial_query` (always computed at depth 30 and mapped down to the index depth) and then exact-filters with `sindist`. The SQLite backing is explicitly documented as an implementation detail. Connections are opened/closed per operation (deliberate, for multithreading) — except the in-memory case, where `memory_db_connection` must stay open for the life of the object or the database vanishes.

   Supporting modules: `qlsc_functions.py` (module-level `distance`/`sindist`/`xy2facenum` wrappers), `utilities.py` (private, package-internal math helpers — docstring warns they are not general purpose), `generate.py` (`sunflower_points_on_sphere`, evenly distributed test points used heavily by the tests).

Ipix numbering follows a z-order curve, so numerically close ipix values are spatially close — this is what makes range-based SQL queries efficient and is the core invariant of the whole scheme.

`astropy` and `pandas` are optional soft dependencies (guarded imports in `qlsc_base.py`): angles/radii may be passed as astropy `Quantity`, and the availability flags gate those paths. The only hard dependency is NumPy.

## Conventions

- Python source is indented with tabs, not spaces. Match it.
- Docstrings are reStructuredText/Sphinx style (`:param:`, `:returns:`) — the API docs are generated from them.
