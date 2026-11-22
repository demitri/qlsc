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

**Releasing:** the version lives only in `source/qlsc/version.py` (kept out of `__init__.py` so `setup.py` can read it without importing dependencies). The commit convention is a dedicated "Bumping version to X.Y.Z." commit. Since binary wheels are built in CI, the release is **tag-first**: the tag is what produces the files that get uploaded, so nothing is uploaded before it exists. The procedure:

1. Update `CHANGELOG.md` and `source/qlsc/version.py`, commit, `git push origin master`.
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

Do **not** upload a locally built `python -m build` result: it produces a wheel for one machine and one interpreter only, and uploading it would claim the release while leaving every other platform on the sdist. A local `cd source && rm -rf dist && python -m build` remains useful for checking that the build works before tagging — just don't upload it. Publishing is deliberately not automated in the workflow (no trusted publishing, no API token in the repo); the workflow marks where such a job would go. The PyPI project page renders the repo-root README (read by `setup.py` at build time, Markdown, absolute image URLs) — it is baked at upload and immutable per release; the v1.1.0 page predates this and shows plain text, which self-corrects at the next upload.

## Architecture

Three layers, from bottom up:

1. **Vendored Q3C C code** — `source/qlsc/c_code/q3c/` (`q3cube.c`, `q3c_poly.c`, headers). This is third-party code from the Q3C project; treat it as an embedded dependency — the four files are upstream Q3C v2.0.6, compiled without PostgreSQL via upstream's own `Q3C_STANDALONE` macro (defined in `setup.py`, which also carries the `Q3C_VERSION` string to keep in step), with exactly ONE local divergence, for a reachable bug (on the next upstream sync, keep it unless the upstream release contains the fix):
   - `q3cube.c`, `q3c_get_nearby`: the loop v2.0.6 added for a cone spilling onto two or three neighbouring faces (`for (;; n0 /= 2)`) clamps only the *upper* cell index to the face. A box touching the far edge of a neighbouring face (xmin or ymin at 0.5, or within double rounding of it) puts the lower index at `n0`, it never equals the clamped upper index, `n0` halves to 0, and the loop spins forever holding the GIL — an uninterruptible hang of `QLSC.nearby_ipix_ranges` for a cone whose bounding box lands, to within double rounding, on the far edge of a neighbouring face (found by the opus review of the 1.2.0 release range). Such cones are constructible — bisecting the radius to the threshold where a cone's box crosses a second face edge hits one about 1.5–2% of the time — but random inputs essentially never do, and which inputs trigger it depends on the platform's libm. The patch clamps the lower index the same way. Guard: `test_nearby_ranges_terminate_on_a_far_face_edge` (`test_nearby_ranges.py`, subprocess with timeout, inputs from both platforms). **Not yet reported upstream**; filing an issue on segasai/q3c is the user's decision. The older sites in `q3c_get_nearby` (main face, `multi_flag == 1`) have the same unclamped lower index; there it would give an out-of-range square rather than a hang, and probing face edges with zero and tiny radii found no case, so it is documented here, not patched.

   qlsc carried two local patches against v2.0.5 — `points[8]` in `q3c_poly.c` (upstream issue #50, a stack buffer overflow reachable via `QLSC.point_in_polygon`) and a `log2` quadtree refinement depth in `q3cube.c` (#51, a signed overflow for tiny query regions) — and v2.0.6 fixed both upstream with different text (`points[8]` with a comment; `q3c_res_depth()`, which gives the same `1 + log2(nside/n0)` qlsc's loop did), so the sync dropped them by hand. Their runtime regressions remain: `test_point_in_polygon.py` (#50) and `test_radial_query_center_coverage.py` (#51). The caller-side half of #50 (5-row, not 3-row, projection arrays) lives in `qlsc_c_module.c`, qlsc's own code; upstream's `q3c.c` uses 5 rows too since its #50 fix (segasai/q3c 5d97e25).

   Upstream issues #50–#55 (all filed from qlsc's audits) are closed as of v2.0.6; #53/#54 were in the un-vendored PostgreSQL layer (`q3c.c`), and from #55 upstream removed the unused `Q3C_BOX_INTERSECT` and gave `q3c_fast_get_xy_minmax`'s empty `Q3C_POLYGON` case a `default: break` (unreachable either way). One #55 item stays as it was, deliberately: `q3c_poly_cover_check`'s concave-polygon over-coverage (upstream added a comment saying full squares need the `q3c_in_poly()` recheck, and changed nothing else). It *is* reachable, via `QLSCIndex.polygon_query`, and is deliberately left unpatched because it cannot change an answer — it only marks coarse pixels "fully covered" that a notch of a concave polygon partly escapes, which enlarges the candidate set, and every candidate is then put through `point_in_polygon` before it is returned. It is over-coverage, never under-coverage; the direction is asserted, with the over-coverage itself demonstrated, in `test_concave_polygon_coarse_coverage_is_generous_but_filtered` (`source/tests/test_index_polygon_query.py`), which is the guard to re-run if this ever looks doubtful. Issue #52 (NaN/Inf reaching undefined float→int conversions) is fixed at qlsc's own boundaries (Python validation plus wrapper checks); upstream's v2.0.6 fix is in the PostgreSQL layer, so qlsc's checks remain the only ones on its paths. The same applies to the vendored functions documented "strictly 0<=ra<360" that do not enforce it (`q3c_get_facenum` and everything reaching it, including both polygon paths via `q3c_get_facenum_poly`): the wrapper unwraps ra before the call. **`qlsc_c_module.c` carries a per-function ra-precondition audit in a block comment at the top of the file — re-run that audit whenever a Q3C function is newly exported or the vendored code is synced**; `source/tests/test_ra_wrapping.py` is its guard. Never edit these files for style or convenience, but a bug that actually impacts qlsc (reachable from the wrapper) DOES get patched locally and immediately — correctness outranks byte-for-byte purity. Document any such divergence here and in the commit, and reconcile it on the next upstream sync. A bug that is dormant in qlsc (unreachable code) is documented only, keeping the files byte-identical. To sync a newer upstream release: copy the four files from https://github.com/segasai/q3c verbatim, check whether any function signature used by `qlsc_c_module.c` changed, and re-apply any local patch upstream has not fixed. **`scripts/sync_q3c.py` mechanizes all of that** — `python3 scripts/sync_q3c.py v2.0.7 --dry-run` first, then without `--dry-run`; it applies each patch by exact match only and stops loudly rather than fuzzing one (including when upstream fixed the bug with different text, which it cannot recognize — then decide by hand, as the v2.0.6 sync did), and it prints the human work (upstream issue re-triage, CLAUDE.md, CHANGELOG) it cannot do. Its `PATCHES` tuple is the second copy of the divergence list above; keep the two in step (`test_registered_patches_are_documented_and_guarded` checks each entry's guard test is named here and exists). Tests: `source/tests/test_sync_q3c.py`. The suite imports the script, so it must stay inside the `python_requires` floor (`test_python_version_floor.py` gates that for the whole tree — it is why `Patch` is a plain class and not a dataclass). All the actual spherical-cube math (ipix encoding/decoding, radial query bin selection) lives here. The central runtime object is the `q3c_prm` struct ("hprm"), which holds precomputed bit-interleaving tables for a given `nside`. Note: as of qlsc v1.1, `QLSC.ipix2ang` returns the pixel *center*, matching upstream's contract change in Q3C v2.0.1 (`ipix2ang_center` is a synonym); pixel corners are available via `ipix2xy` and `ipix2polygon`.

2. **CPython extension** — `source/qlsc/c_code/qlsc_c_module.c`, compiled by `setup.py` into the extension module `qlsc.q3c`. It wraps the Q3C functions and passes the `q3c_prm` struct around as an opaque PyCapsule (`init_q3c()` creates it; every geometry call takes it back). Functions accept scalars or NumPy arrays (arrays must be `dtype=double`; other dtypes get copied).

3. **Python API** — `source/qlsc/qlsc_base.py` defines the two public classes:
   - **`QLSC(depth)`** — pure geometry: conversions among (ra,dec), ipix, and per-face (x,y) coordinates (`ang2ipix`, `ipix2ang_center`, `ipix2polygon`, `ipix_up`/`ipix_down`, …), plus the two scheme-independent containment tests `point_in_polygon` and `point_in_ellipse` that the index's region queries use as their exact filters. `depth` ∈ [0,30]; each depth level subdivides every bin into four; depth 30 is the Q3C/PostgreSQL default. The face-plane (x,y) coordinate system (origin at face center, axes −1..+1) is the recommended working space for nontrivial calculations.
   - **`QLSCIndex(qlsc, filepath)`** — a region-search index backed by SQLite (in-memory by default, or a file; table `qlsc_ipix`). Points go in via `add_point`/`add_points` (with optional string keys). The three region queries — `radial_query()`, `polygon_query()`, `ellipse_query()` — share one shape: coarse ipix ranges from the C layer (`q3c.radial_query`/`poly_query`/`ellipse_query`, always computed at depth 30 and mapped down to the index depth by `_index_depth_ranges`), SQL over those ranges, then an exact per-candidate filter in Python (`sindist`, `QLSC.point_in_polygon`, `QLSC.point_in_ellipse`). The coarse ranges are always a superset — that is the invariant everything rests on, and it is why an over-covering coarse step is harmless. For the cone and the ellipse the search region is grown (floor 1e-6 deg, ×10 on empty output or a range overflow) because q3c produces no or incomplete ranges for degenerate regions; a polygon has no such knob, so an empty range set falls back to a full scan with a logged warning. Two contracts to preserve: (1) **`radial_query` alone squeezes its plain result**, so one match is `(2,)` not `(1,2)`; that wart is queued for 2.0 (see `AI/FEATURE_IDEAS.md`) and every query added since returns `(n,2)` always — the choice is carried explicitly by the `MatchArrayShape` enum, not by accident. (2) **Q3C's ellipse is a conic in the gnomonic projection about its center** (semi-axes `tan(a)` and `tan(a)·axis_ratio`, verified numerically in `test_ellipse_geometry.py`), so it only means anything inside the near hemisphere: `axis_ratio < 1` is capped at `ELLIPSE_MAXIMUM_SEMI_MAJOR_AXIS` = 90 deg with a loud `ValueError`, while `axis_ratio == 1` is a circle that `ellipse_query` answers by *calling* `radial_query` — never by q3c's conic, which rejects everything past 90 deg and is not monotonic in the semi-major axis beyond it. The polygon exact filter runs as one batched C call (`q3c.check_sphere_points_in_poly`, q3c's `invocation=1` path) so the polygon is projected onto the cube faces once rather than once per candidate. The SQLite backing is explicitly documented as an implementation detail. Connections are opened/closed per operation (deliberate, for multithreading) — except the in-memory case, where `memory_db_connection` must stay open for the life of the object or the database vanishes.

   Supporting modules: `qlsc_functions.py` (module-level `distance`/`sindist`/`xy2facenum` wrappers), `utilities.py` (private, package-internal math helpers — docstring warns they are not general purpose), `generate.py` (`sunflower_points_on_sphere`, evenly distributed test points used heavily by the tests).

Ipix numbering follows a z-order curve, so numerically close ipix values are spatially close — this is what makes range-based SQL queries efficient and is the core invariant of the whole scheme.

`astropy` and `pandas` are optional soft dependencies (guarded imports in `qlsc_base.py`): angles/radii may be passed as astropy `Quantity`, and the availability flags gate those paths. The only hard dependency is NumPy.

## Conventions

- Python source is indented with tabs, not spaces. Match it.
- Docstrings are reStructuredText/Sphinx style (`:param:`, `:returns:`) — the API docs are generated from them.
