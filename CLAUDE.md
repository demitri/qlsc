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

**Releasing:** the version lives only in `source/qlsc/version.py` (kept out of `__init__.py` so `setup.py` can read it without importing dependencies). The commit convention is a dedicated "Bumping version to X.Y.Z." commit.

## Architecture

Three layers, from bottom up:

1. **Vendored Q3C C code** — `source/qlsc/c_code/q3c/` (`q3cube.c`, `q3c_poly.c`, headers). This is third-party code from the Q3C project; treat it as an embedded dependency and never edit it locally — the four files are byte-for-byte upstream Q3C v2.0.5, compiled without PostgreSQL via upstream's own `Q3C_STANDALONE` macro (defined in `setup.py`, which also carries the `Q3C_VERSION` string to keep in step). To sync a newer upstream release: copy the four files from https://github.com/segasai/q3c verbatim, then check whether any function signature used by `qlsc_c_module.c` changed. Known upstream issue as of v2.0.5, deliberately not patched locally because the code is unreachable from qlsc (the wrapper exports no polygon functions): `q3c_check_sphere_point_in_poly` in `q3c_poly.c` declares `points[4]` where large polygons need `points[8]` — a stack overflow in q3c proper; confirm it's fixed (or still dormant here) on the next sync. All the actual spherical-cube math (ipix encoding/decoding, radial query bin selection) lives here. The central runtime object is the `q3c_prm` struct ("hprm"), which holds precomputed bit-interleaving tables for a given `nside`. Note: as of qlsc v1.1, `QLSC.ipix2ang` returns the pixel *center*, matching upstream's contract change in Q3C v2.0.1 (`ipix2ang_center` is a synonym); pixel corners are available via `ipix2xy` and `ipix2polygon`.

2. **CPython extension** — `source/qlsc/c_code/qlsc_c_module.c`, compiled by `setup.py` into the extension module `qlsc.q3c`. It wraps the Q3C functions and passes the `q3c_prm` struct around as an opaque PyCapsule (`init_q3c()` creates it; every geometry call takes it back). Functions accept scalars or NumPy arrays (arrays must be `dtype=double`; other dtypes get copied).

3. **Python API** — `source/qlsc/qlsc_base.py` defines the two public classes:
   - **`QLSC(depth)`** — pure geometry: conversions among (ra,dec), ipix, and per-face (x,y) coordinates (`ang2ipix`, `ipix2ang_center`, `ipix2polygon`, `ipix_up`/`ipix_down`, …). `depth` ∈ [0,30]; each depth level subdivides every bin into four; depth 30 is the Q3C/PostgreSQL default. The face-plane (x,y) coordinate system (origin at face center, axes −1..+1) is the recommended working space for nontrivial calculations.
   - **`QLSCIndex(qlsc, filepath)`** — a cone-search index backed by SQLite (in-memory by default, or a file; table `qlsc_ipix`). Points go in via `add_point`/`add_points` (with optional string keys); `radial_query()` builds SQL from ipix ranges produced by the C-level `radial_query_it` and then exact-filters with `sindist`. The SQLite backing is explicitly documented as an implementation detail. Connections are opened/closed per operation (deliberate, for multithreading) — except the in-memory case, where `memory_db_connection` must stay open for the life of the object or the database vanishes.

   Supporting modules: `qlsc_functions.py` (module-level `distance`/`sindist`/`xy2facenum` wrappers), `utilities.py` (private, package-internal math helpers — docstring warns they are not general purpose), `generate.py` (`sunflower_points_on_sphere`, evenly distributed test points used heavily by the tests).

Ipix numbering follows a z-order curve, so numerically close ipix values are spatially close — this is what makes range-based SQL queries efficient and is the core invariant of the whole scheme.

`astropy` and `pandas` are optional soft dependencies (guarded imports in `qlsc_base.py`): angles/radii may be passed as astropy `Quantity`, and the availability flags gate those paths. The only hard dependency is NumPy.

## Conventions

- Python source is indented with tabs, not spaces. Match it.
- Docstrings are reStructuredText/Sphinx style (`:param:`, `:returns:`) — the API docs are generated from them.
