
import pytest

from qlsc import QLSC

# depth, ra, dec, facenum
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

@pytest.mark.parametrize("depth, ra, dec, facenum", expected_results)
def test_facenum(depth, ra, dec, facenum):
	'''
	Test QLSC ang2ipix.
	'''
	q = QLSC(depth=depth)
	assert facenum == q.face_number(ra, dec)
