
import ast
import pathlib

import numpy as np
import pytest

from qlsc import QLSC

# A gate over docs/generate_figures.py: no figure may take a single pixel's area as the
# area of a pixel.
#
# This exists because the mistake was made three times in one file and is invisible in the
# output. QLSC pixels vary in area by a factor of about 5.2, and pixel 0 sits at a cube
# VERTEX, which makes ipix_area(0, depth) the SMALLEST area in the scheme — about 2.55
# times below the mean. Both density figures divided their counts by it, overstating every
# density by that factor while the page beside them said the mean had been used, and the
# pixel-scale figure plotted sqrt(ipix_area(0, d)) as "the pixel side" against a page that
# tabulates the mean.
#
# Two scales are defensible, and generate_figures.py names both: mean_pixel_area_sq_deg
# for anything characterising a depth, and pixel_areas_sq_deg for anything showing a
# density. ipix_area with a literal pixel number is neither.
#
# Blind spot, stated rather than assumed covered: this checks the SPELLING, not the
# arithmetic. A figure that divides by the mean where it should divide per pixel passes
# here; that the figures agree with the prose beside them remains a human-review duty,
# and docs/check_example_outputs.py covers only the printed code blocks, not the images.

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
FIGURE_SCRIPT = REPO_ROOT / "docs" / "generate_figures.py"

def ipix_area_calls(tree:ast.AST):
	'''
	Yield every call of an ``ipix_area`` method in the tree, with its line number.
	'''
	for node in ast.walk(tree):
		if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) \
		   and node.func.attr == "ipix_area":
			yield node

def test_the_figure_script_is_present():
	'''
	The gate below is vacuous if the file moves; say so rather than passing quietly.
	'''
	assert FIGURE_SCRIPT.is_file(), f"expected the figure script at {FIGURE_SCRIPT}"

def test_no_figure_uses_one_pixels_area_as_the_pixel_area():
	'''
	Test that no ipix_area call in the figure script names a literal pixel.

	A literal pixel number is always the tell: a figure asking for one specific pixel's
	area is either characterising the scheme (and wants the mean) or scaling a count (and
	wants that count's own pixel). Reading the pixel number from a variable is how the
	per-pixel helper does it, and is fine.
	'''
	tree = ast.parse(FIGURE_SCRIPT.read_text(encoding="utf-8"), str(FIGURE_SCRIPT))

	problems = []
	for node in ipix_area_calls(tree):
		first = node.args[0] if node.args else None
		if isinstance(first, ast.Constant):
			problems.append(f"{FIGURE_SCRIPT.name}:{node.lineno}: ipix_area({first.value!r}, ...) "
							f"takes one pixel's area as the area of a pixel")

	assert not problems, (
		"\n  ".join(["a figure is scaling by a single pixel's area:"] + problems)
		+ "\nUse mean_pixel_area_sq_deg() to characterise a depth, or pixel_areas_sq_deg() "
		  "to turn counts into a density.")

def test_the_gate_detects_a_violation():
	'''
	Test that the check above actually fires; a gate that cannot fail guards nothing.
	'''
	sample = ast.parse("area = q.ipix_area(0, depth)\ngood = q.ipix_area(int(p), q.depth)\n")
	literal = [n for n in ipix_area_calls(sample) if isinstance(n.args[0], ast.Constant)]

	assert len(literal) == 1
	assert literal[0].args[0].value == 0

def test_pixel_zero_really_is_the_smallest_pixel():
	'''
	Test the fact this gate rests on: ipix_area(0, depth) is not a representative area.

	If this ever stops being true — if pixel 0 stops being a cube-vertex pixel, or the
	pixels become genuinely equal-area — the gate above is guarding against nothing and
	should be revisited rather than left in place.
	'''
	ratios = {}
	for depth in (2, 3, 5):
		q = QLSC(depth=depth)
		areas = np.array([q.ipix_area(p, depth) for p in range(q.nbins)])
		mean = 4. * np.pi / q.nbins

		# pixel 0 is a cube-vertex pixel, and those are the smallest in the scheme;
		# there are 24 of them, so it is tied for the minimum rather than uniquely it
		assert q.ipix_area(0, depth) == pytest.approx(areas.min(), rel=1e-9)
		assert q.ipix_area(0, depth) < mean
		ratios[depth] = (areas.min() / mean, areas.max() / areas.min())

	# The spread widens with depth towards the factor of about 5.2 the documentation
	# quotes; it is not yet there at these depths, which is why this asserts the trend
	# and the depth-5 value the figures are drawn at rather than one constant.
	assert ratios[2][1] < ratios[3][1] < ratios[5][1]
	assert ratios[5][1] > 4.8

	# the number the prose quotes: at depth 5, dividing by a cube-vertex pixel's area
	# rather than the mean overstates a density by about 2.55
	assert 1. / ratios[5][0] == pytest.approx(2.55, abs=0.02)

def test_the_named_helpers_mean_what_they_say():
	'''
	Test the two helpers the figures use, since the figures themselves are images and
	nothing else checks their arithmetic.
	'''
	import importlib.util

    # the script is not importable as part of the package; load it by path
	spec = importlib.util.spec_from_file_location("qlsc_generate_figures", FIGURE_SCRIPT)
	if spec is None or spec.loader is None:
		pytest.skip("the figure script could not be loaded")
	module = importlib.util.module_from_spec(spec)
	try:
		spec.loader.exec_module(module)
	except ImportError as e:
		pytest.skip(f"the figure script needs matplotlib and astropy, which are not installed here ({e})")

	q = QLSC(depth=4)
	sq_deg = (180. / np.pi) ** 2

	assert module.mean_pixel_area_sq_deg(q) == pytest.approx(4. * np.pi / q.nbins * sq_deg)

	pixels = [0, 1, q.nbins // 2]
	areas = module.pixel_areas_sq_deg(q, pixels)
	assert areas.shape == (3,)
	for area, pixel in zip(areas, pixels):
		assert area == pytest.approx(q.ipix_area(pixel, q.depth) * sq_deg)

	# a scalar is accepted as one pixel, as it is everywhere else in qlsc
	assert module.pixel_areas_sq_deg(q, 3).shape == (1,)
