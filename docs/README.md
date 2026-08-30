# Building and checking the documentation

The published site is at <https://qlsc.readthedocs.io>. Read the Docs builds it
from `.readthedocs.yaml` at the repository root, which installs
`docs/requirements.txt` and runs Sphinx over `docs/conf.py`. The C extension is
*not* built there: `conf.py` mocks `qlsc.q3c` so that autodoc can import the
package without a compiler.

## Build

```bash
cd docs
make html                       # output in docs/_build/html
python -m sphinx -b html . _build/html -W --keep-going   # warnings as errors
```

The build should be warning-free. One `AstropyDeprecationWarning` about
`astropy.samp` may appear if astropy is installed in the build environment; it is
emitted before Sphinx starts counting, does not fail `-W`, and does not occur on
Read the Docs, where astropy is not installed.

## Checking the examples

Every `code-block:: text` that follows a `code-block:: python` in these pages is
real output from really running that code. `check_example_outputs.py` re-runs the
code and diffs it, so that promise cannot quietly stop being true:

```bash
cd docs
python check_example_outputs.py             # every page, offline
python check_example_outputs.py --network   # also the page that queries VizieR
python check_example_outputs.py --list      # coverage, and every skip with its reason
python check_example_outputs.py -v examples/cone_search.rst
```

It exits non-zero on any mismatch. When it reports one, **re-run the code and
paste what it really prints** — never edit the expected output by hand, which is
how the two drifted apart in the first place.

Blocks that cannot run hermetically are marked in the page source with a reST
comment giving the reason (`.. check-outputs: run-only -- ...`,
`skip`, `setup`, and the page-level `.. check-outputs-page: needs-network -- ...`).
The script's module docstring lists what it does *not* check — unpaired code
blocks, numbers quoted in prose, and the figures — which remain a manual duty.

It needs the packages the examples use: numpy, astropy, matplotlib and a built
`qlsc`. Those are not in `docs/requirements.txt`, because the site itself does not
need them.

## Figures

Every image under `_static/figures` is produced by `generate_figures.py` from
seeded synthetic catalogs; none are drawn by hand. Regenerate them all, or one at
a time:

```bash
cd docs
python generate_figures.py
python generate_figures.py sky_density cube_faces
```

The seeds match the ones in the pages, so the numbers a figure shows agree with
the numbers the surrounding text quotes. That agreement is *not* checked
automatically — if you change a figure's catalog, re-read the pages that cite it.
