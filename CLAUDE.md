# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

QLSC is a Python package (published on PyPI as `qlsc`) implementing the quadrilateralized spherical cube — a scheme for segmenting a sphere into approximately equal-area pixels, aimed at astronomical spatial indexing and cone searches. It embeds C code vendored from Sergey Koposov's Q3C PostgreSQL extension (copyright Koposov) and provides Q3C's catalog-indexing functionality without requiring PostgreSQL. Docs are on Read the Docs (https://qlsc.readthedocs.io); the README is the primary conceptual introduction to the scheme itself.

## Layout

- Everything installable is under `source/` (`source/setup.py`, `source/qlsc/`, `source/tests/`), not at the repo root.
- The repo root holds the README, `docs/`, examples, figures (LaTeX/graphics sources for the README diagrams), and `c_test/` (an Xcode project for exercising the C code directly).
- `scripts/` holds repository tools that are not part of the package (`sync_q3c.py`, evidence scripts).
- `AI/` is gitignored: private working notes, not versioned. Nothing a clone or a reviewer needs may live only there.

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

# Check that every documented example still prints what its page says; then the strict build
cd docs && python check_example_outputs.py && python -m sphinx -W --keep-going -b html . _build/html
```

**Test gotcha:** the tests do `import qlsc`, and since the source tree contains no built `.so`, they run against the *installed* package, not the working tree. After changing anything under `source/qlsc/`, reinstall before running tests or you will test stale code. Test files prefixed `xtest_` (rather than `test_`) are deliberately disabled.

**Full-coverage runs** (the global rule puts test suites on a server) use r740's `~/venvs/qlsc`, built from `/usr/local/anaconda/bin/python` with numpy, astropy, pandas, matplotlib, mpmath, pytest, pytest-xdist, build and twine. Without astropy, matplotlib or cornish, tests skip by name. Install cornish from `$GH/cornish/source`, not PyPI: PyPI's cornish 1.1 fails to import (`No module named 'cornish.channel'`). Copy the working tree there (rsync excluding `.git`, `AI/` and the build directory), `pip install --force-reinstall --no-deps --no-build-isolation ./source`, then `pytest -n 24` from `source/tests`. A clean run has 0 skipped.

**Documentation outputs** are checked by `docs/check_example_outputs.py` on Linux CI, so a printed value must not depend on the platform: print floats to 12 significant digits, never a full-precision repr. Outputs that print the qlsc version change with every version bump; regenerate them by running the blocks.

**Python floor:** `python_requires='>=3.7'`. The tests import `scripts/sync_q3c.py`, so everything in the tree must stay within the floor; `test_python_version_floor.py` gates syntax and stdlib modules for the whole tree, and also requires every text read or write to pass `encoding=` (Windows CI decodes cp1252 otherwise).

## Releasing

The version lives only in `source/qlsc/version.py` (kept out of `__init__.py` so `setup.py` can read it without importing dependencies). The commit convention is a dedicated "Bumping version to X.Y.Z." commit. Binary wheels are built in CI, so the release is **tag-first**: the tag is what produces the files that get uploaded, and nothing is uploaded before it exists.

1. Update `CHANGELOG.md` (the version heading carries the tag date) and `source/qlsc/version.py`, commit, `git push origin master`. Before tagging, a `gh workflow run wheels.yml --ref master` exercises every platform without a tag.
2. `git tag vX.Y.Z && git push origin vX.Y.Z` — the tag triggers `.github/workflows/wheels.yml`, which builds and tests wheels for every supported platform and builds the sdist.
3. Wait for that run to go green, including its `verify` job: find it with `gh run list --workflow wheels.yml --limit 5` (the tag run is the one whose event is `push`, at the tag's SHA) and follow it with `gh run watch <run-id>`. `verify` is the gate — it runs `twine check --strict` over the whole set, checks that every expected interpreter/platform wheel is present (`.github/scripts/verify_artifacts.py`), that the tag matches `version.py`, and that the wheels import under the oldest supported NumPy.
4. Download that run's artifacts into a *fresh* staging directory and collect the files into a *separate fresh* release directory (`gh run download` writes one subdirectory per artifact, so the files have to be gathered; the directories must not be left among them, and a reused directory would smuggle a stale local build into the upload):

   ```bash
   cd <repo root>
   rm -rf release-staging release-dist && mkdir release-dist
   gh run download <run-id> -D release-staging
   find release-staging -type f \( -name '*.whl' -o -name '*.tar.gz' \) -exec mv -n {} release-dist/ \;
   ```

5. Verify locally, then upload — always naming the two globs rather than `release-dist/*`, so anything that is not a wheel or an sdist can never reach twine:

   ```bash
   python .github/scripts/verify_artifacts.py release-dist   # needs "packaging"
   twine check --strict release-dist/*.whl release-dist/*.tar.gz
   twine upload --repository qlsc release-dist/*.whl release-dist/*.tar.gz
   ```

   (the `[qlsc]` section in `~/.pypirc`; twine ≥6 requires sections be listed under its `[distutils] index-servers` header). Wheels and sdist go up together, in one upload. `rm -rf release-staging release-dist` afterwards.
6. `gh release create vX.Y.Z` with the version's CHANGELOG section as notes.

**If the tag build fails**, nothing has been published — the upload is step 5 and the tag only triggers a build — so the tag itself is the only thing to clean up, and the version number is still free:

```bash
git tag -d vX.Y.Z && git push origin :refs/tags/vX.Y.Z   # then fix, and re-tag the new commit
```

Re-tagging the same version is fine *as long as nothing was uploaded*. Once `twine upload` has succeeded, even partially, that version is spent: PyPI refuses a second upload of a filename it already has, and deleting a release to re-upload is not an option (it breaks anyone who pinned it). So if a failure appears after an upload has begun, do not re-tag — finish the interrupted upload with `twine upload --skip-existing …` if the artifacts are good, or bump to the next patch version and start again if they are not.

Do **not** upload a locally built `python -m build` result: it produces a wheel for one machine and one interpreter only, and uploading it would claim the release while leaving every other platform on the sdist. A local `cd source && rm -rf dist && python -m build` is useful for checking that the build works before tagging — just don't upload it. Publishing is deliberately not automated in the workflow (no trusted publishing, no API token in the repo); the workflow marks where such a job would go. The PyPI project page renders the repo-root README (read by `setup.py` at build time, Markdown, absolute image URLs); it is baked at upload and immutable per release.

## Architecture

Three layers, from bottom up.

### 1. Vendored Q3C C code

`source/qlsc/c_code/q3c/` holds four files (`q3cube.c`, `q3c_poly.c`, `common.h`, `my_bits.h`): upstream Q3C v2.0.6 plus the local patches below, compiled without PostgreSQL via upstream's own `Q3C_STANDALONE` macro. `setup.py` defines that macro and carries the `Q3C_VERSION` string, which must match the vendored release. All the spherical-cube math (ipix encoding/decoding, region decomposition) lives here. The central runtime object is the `q3c_prm` struct ("hprm"), which holds precomputed bit-interleaving tables for a given `nside`.

**Local patches.** Exactly one:

- `q3cube.c`, `q3c_get_nearby`: the loop for a cone spilling onto two or three neighbouring faces (`for (;; n0 /= 2)`) clamps the *lower* cell index to the face as well as the upper. Upstream clamps only the upper one. A box touching the far edge of a neighbouring face (xmin or ymin at 0.5, or within double rounding of it) then never fits one cell, `n0` reaches 0 and the loop never ends, hanging `QLSC.nearby_ipix_ranges` with the GIL held. Which inputs trigger it depends on the platform's libm. Upstream issue: segasai/q3c#60. Guard: `test_nearby_ranges_terminate_on_a_far_face_edge` (`test_nearby_ranges.py`; each cone runs in a subprocess with a timeout). Drop the patch at the sync that brings upstream's fix. The other two sites in `q3c_get_nearby` (main face, `multi_flag == 1`) also leave the lower index unclamped; there the effect would be an out-of-range square rather than a hang, no input reaching it has been found, and they are not patched.

**Rules for the vendored files:**

- Never edit them for style or convenience.
- A bug reachable from the wrapper is patched locally, immediately: register it in `PATCHES` in `scripts/sync_q3c.py`, give it a runtime guard test, list it above, and report it upstream. `test_registered_patches_are_documented_and_guarded` checks that each `PATCHES` entry's guard test is named in this file and exists.
- A bug unreachable from the wrapper is documented, not patched, so the files stay byte-identical to upstream.
- The caller-side buffers in `qlsc_c_module.c` for polygon projection are `[5][Q3C_MAX_N_POLY_VERTEX]` (the main face plus up to four neighbours); keep them at 5 rows.

**Known upstream behaviour, deliberately not patched:** `q3c_poly_cover_check` can mark a coarse pixel "fully covered" although a notch of a concave polygon escapes it. This is reachable via `QLSCIndex.polygon_query` but cannot change an answer: it only enlarges the candidate set, and every candidate goes through `point_in_polygon`. The direction (over-coverage, never under-coverage) is asserted in `test_concave_polygon_coarse_coverage_is_generous_but_filtered` (`source/tests/test_index_polygon_query.py`).

**Preconditions qlsc enforces itself:**

- NaN/Inf inputs are rejected at qlsc's own boundaries (Python validation plus wrapper checks); the vendored C does not reject them on qlsc's paths.
- Several vendored functions documented "strictly 0<=ra<360" do not enforce it (`q3c_get_facenum` and everything reaching it, including both polygon paths via `q3c_get_facenum_poly`); the wrapper unwraps ra before the call. **`qlsc_c_module.c` carries a per-function ra-precondition audit in a block comment at the top of the file — re-run it whenever a Q3C function is newly exported or the vendored code is synced.** `source/tests/test_ra_wrapping.py` is its guard.

**Syncing to a newer upstream release:** `python3 scripts/sync_q3c.py vX.Y.Z --dry-run`, then without `--dry-run`. The script fetches the tag, confirms the vendored files are the release `Q3C_VERSION` names plus exactly the registered patches, re-applies each patch by exact match (stopping, never fuzzing, when an anchor does not match), diffs the `common.h` interface the wrapper uses, updates `Q3C_VERSION`, reinstalls and runs the suite, and prints the human work it cannot do. When it stops on a patch, read the new upstream code: if upstream fixed the bug (often with different text, which exact matching cannot recognize), drop the patch from `PATCHES`, its entry above and the reference to it in its guard test's comments. Also re-triage the open upstream issues for newly reachable bugs, re-run the ra audit, and update this file and the CHANGELOG. Tests: `source/tests/test_sync_q3c.py`.

**Pixel centres:** `QLSC.ipix2ang` returns the pixel *center*, matching Q3C's contract (`ipix2ang_center` is a synonym); pixel corners come from `ipix2xy` and `ipix2polygon`.

### 2. CPython extension

`source/qlsc/c_code/qlsc_c_module.c`, compiled by `setup.py` into the extension module `qlsc.q3c`. It wraps the Q3C functions and passes the `q3c_prm` struct around as an opaque PyCapsule (`init_q3c()` creates it; every geometry call takes it back). Functions accept scalars or NumPy arrays (arrays must be `dtype=double`; other dtypes get copied).

### 3. Python API

`source/qlsc/qlsc_base.py` defines the two public classes:

- **`QLSC(depth)`** — pure geometry: conversions among (ra,dec), ipix, and per-face (x,y) coordinates (`ang2ipix`, `ipix2ang_center`, `ipix2polygon`, `ipix_up`/`ipix_down`, …), plus the two scheme-independent containment tests `point_in_polygon` and `point_in_ellipse` that the index's region queries use as their exact filters. `depth` ∈ [0,30]; each depth level subdivides every bin into four; depth 30 is the Q3C/PostgreSQL default. The face-plane (x,y) coordinate system (origin at face center, axes −1..+1) is the recommended working space for nontrivial calculations.
- **`QLSCIndex(qlsc, filepath)`** — a region-search index backed by SQLite (in-memory by default, or a file; table `qlsc_ipix`). Points go in via `add_point`/`add_points` (with optional string keys). The three region queries — `radial_query()`, `polygon_query()`, `ellipse_query()` — share one shape: coarse ipix ranges from the C layer (`q3c.radial_query`/`poly_query`/`ellipse_query`, always computed at depth 30 and mapped down to the index depth by `_index_depth_ranges`), SQL over those ranges, then an exact per-candidate filter in Python (`sindist`, `QLSC.point_in_polygon`, `QLSC.point_in_ellipse`). The coarse ranges are always a superset — that is the invariant everything rests on, and it is why an over-covering coarse step is harmless. For the cone and the ellipse the search region is grown (floor 1e-6 deg, ×10 on empty output or a range overflow) because q3c produces no or incomplete ranges for degenerate regions; a polygon has no such knob, so an empty range set falls back to a full scan with a logged warning. The polygon exact filter runs as one batched C call (`q3c.check_sphere_points_in_poly`, q3c's `invocation=1` path) so the polygon is projected onto the cube faces once rather than once per candidate. The SQLite backing is explicitly documented as an implementation detail. Connections are opened/closed per operation (deliberate, for multithreading) — except the in-memory case, where `memory_db_connection` must stay open for the life of the object or the database vanishes.

Two contracts to preserve:

1. **`radial_query` alone squeezes its plain result**, so one match is `(2,)` not `(1,2)`. Every other query returns `(n,2)` always; the choice is carried explicitly by the `MatchArrayShape` enum. Removing the squeeze is a breaking change reserved for 2.0 (see `AI/FEATURE_IDEAS.md`).
2. **Q3C's ellipse is a conic in the gnomonic projection about its center** (semi-axes `tan(a)` and `tan(a)·axis_ratio`, verified numerically in `test_ellipse_geometry.py`), so it only means anything inside the near hemisphere: `axis_ratio < 1` is capped at `ELLIPSE_MAXIMUM_SEMI_MAJOR_AXIS` = 90 deg with a loud `ValueError`, while `axis_ratio == 1` is a circle that `ellipse_query` answers by *calling* `radial_query` — never by q3c's conic, which rejects everything past 90 deg and is not monotonic in the semi-major axis beyond it.

Supporting modules: `qlsc_functions.py` (module-level `distance`/`sindist`/`xy2facenum` wrappers), `utilities.py` (private, package-internal math helpers — docstring warns they are not general purpose), `generate.py` (`sunflower_points_on_sphere`, evenly distributed test points used heavily by the tests).

Ipix numbering follows a z-order curve, so numerically close ipix values are spatially close — this is what makes range-based SQL queries efficient and is the core invariant of the whole scheme.

`astropy` and `pandas` are optional soft dependencies (guarded imports in `qlsc_base.py`): angles/radii may be passed as astropy `Quantity`, and the availability flags gate those paths. The only hard dependency is NumPy.

## Conventions

- Python source is indented with tabs, not spaces. Match it.
- Docstrings are reStructuredText/Sphinx style (`:param:`, `:returns:`) — the API docs are generated from them.
- Comments, docstrings and this file state what is true now and what to do. History (a bug's earlier behaviour, the release that fixed it, who found it, when) belongs in commit messages and the CHANGELOG.
