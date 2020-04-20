
import pytest

from pyq3c import Q3C

# bin_level, ra, dec, facenum
expected_results = [
	(30, 0, 0, 1),
	(30, 45, 0, 2),
	(30, 95, 0, 2),
	(30, 134.9, 0, 2),
	(30, 136, 0, 3),
	(30, 220, 0, 3),
	(30, 275, 0, 4),
	(30, 10, 89, 0),
	(30, 225, 89, 0),
	(30, 0, -45, 1),
	(30, 0, -46, 5),
	(30, 99, -89, 5)
]

@pytest.mark.parametrize("bin_level, ra, dec, facenum", expected_results)
def test_facenum(bin_level, ra, dec, facenum):
	'''
	Test Q3C ang2ipix.
	'''
	qlsc = Q3C(bin_level=bin_level)
	assert facenum == qlsc.face_number(ra, dec)
