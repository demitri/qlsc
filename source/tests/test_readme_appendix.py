
import math
import pathlib
import re

import pytest

from qlsc import QLSC

# The README's appendix tabulates, for every depth, the number of bins and the mean
# bin area in three units. Those are all computable, so they are checked here rather
# than trusted.
#
# This exists because the table was wrong: the column headed "Bin Area (square µas)"
# held square MILLIarcseconds, a factor of 10^6 out, and the same mistaken figure was
# quoted twice in the prose as the resolution of a depth-30 pixel. Nothing caught it
# because nothing compared the table to the code.
#
# The README lives at the repository root and is not shipped in the sdist, so this test
# skips (with a reason) when run from an unpacked distribution.

README_PATH = pathlib.Path(__file__).resolve().parent.parent.parent / "README.md"

if not README_PATH.is_file():
	pytest.skip(f"README.md is not present at {README_PATH}; it ships in the repository, "
				"not in the qlsc package.", allow_module_level=True)

README = README_PATH.read_text()

STERADIANS_PER_SQUARE_DEGREE = (math.pi / 180) ** 2
SQUARE_MAS_PER_SQUARE_DEGREE = (3600 * 1000) ** 2

def appendix_rows() -> list:
	'''
	Parse the appendix table into (depth, bins per face, bins on sphere, area in
	steradians, area in square mas, area in square degrees) tuples.
	'''
	rows = []
	for line in README.splitlines():
		match = re.match(r"^\|\s*(\d+)\s*\|(.+)\|\s*$", line)
		if match is None:
			continue
		fields = [field.strip().replace(",", "") for field in match.group(2).split("|")]
		if len(fields) != 5:
			continue
		rows.append((int(match.group(1)), int(fields[0]), int(fields[1]),
					 float(fields[2]), float(fields[3]), float(fields[4])))
	return rows

def test_the_table_was_found():
	'''
	Test that the table parsed at all; every check below is vacuous otherwise.
	'''
	rows = appendix_rows()
	assert len(rows) == 31, f"expected one row per depth 0-30, parsed {len(rows)}"
	assert [row[0] for row in rows] == list(range(31))

@pytest.mark.parametrize("row", appendix_rows(), ids=lambda row: f"depth{row[0]}")
def test_appendix_row(row):
	'''
	Test one row of the appendix against the scheme it describes.
	'''
	depth, bins_per_face, bins_on_sphere, steradians, square_mas, square_degrees = row
	qlsc = QLSC(depth=depth)

	assert bins_per_face == qlsc.nside * qlsc.nside
	assert bins_on_sphere == qlsc.nbins

	mean_area = 4 * math.pi / qlsc.nbins
	assert steradians == pytest.approx(mean_area, rel=1e-6)
	assert square_degrees == pytest.approx(mean_area / STERADIANS_PER_SQUARE_DEGREE, rel=1e-6)
	assert square_mas == pytest.approx(
		mean_area / STERADIANS_PER_SQUARE_DEGREE * SQUARE_MAS_PER_SQUARE_DEGREE, rel=1e-6), \
		"the area column is square milliarcseconds; it was once labelled square microarcseconds"

def test_the_area_column_is_not_square_microarcseconds():
	'''
	Test the specific error directly: read as square microarcseconds the column is out
	by a factor of a million.
	'''
	depth, _, _, steradians, square_mas, _ = appendix_rows()[30]
	assert depth == 30

	square_degrees = steradians / STERADIANS_PER_SQUARE_DEGREE
	true_square_uas = square_degrees * (3600 * 1e6) ** 2

	assert square_mas == pytest.approx(true_square_uas / 1e6, rel=1e-6)
	assert not square_mas == pytest.approx(true_square_uas, rel=1e-3)

def test_the_prose_resolution_figure_matches_the_table():
	'''
	Test that the depth-30 figure quoted in the prose is the one in the table, and in
	the units the table is in. It was quoted twice as "~0.08 µ"", which is the
	milliarcsecond value wearing the microarcsecond label.
	'''
	_, _, _, steradians, square_mas, _ = appendix_rows()[30]

	assert square_mas == pytest.approx(0.077, abs=0.005), "the depth-30 mean bin area is ~0.08 square mas"
	assert math.sqrt(square_mas) == pytest.approx(0.278, abs=0.005), "~0.28 mas on a side"

	# both prose mentions now name milliarcseconds
	assert README.count("square milliarcseconds") == 2
	assert 'µ"' not in README and 'μ"' not in README, \
		'the prose still quotes the depth-30 pixel size in microarcseconds'

def test_the_readme_does_not_claim_equal_area_pixels():
	'''
	Test that the README no longer says the pixels "are equal to within a few percent".
	They differ by a factor of 3*sqrt(3); the sentence is replaced by that number.
	'''
	assert "equal to within a few percent" not in README
	assert "5.196" in README

def test_the_claimed_area_ratio_is_the_measured_one():
	'''
	Test the figure the README now quotes for the spread, against the pixels themselves.
	'''
	qlsc = QLSC(depth=12)
	largest = qlsc.ipix_area(qlsc.ang2ipix(0., 0.))                 # centre of a cube face
	smallest = qlsc.ipix_area(qlsc.ang2ipix(44.9999999, 35.2643896)) # a cube vertex

	assert largest / smallest == pytest.approx(3 * math.sqrt(3), abs=3e-3)

	coarse = QLSC(depth=2)
	areas = [coarse.ipix_area(ipix) for ipix in range(coarse.nbins)]
	assert max(areas) / min(areas) == pytest.approx(2.5, abs=0.05), \
		'the README says the spread is "already a factor of 2.5 at depth=2"'
