
import pytest

from pyq3c import Q3C

# ra, dec, facenum
expected_results = [
	(0, 0, 1),
	(45, 0, 2),
	(95, 0, 2),
	(134.9, 0, 2),
	(136, 0, 3),
	(220, 0, 3),
	(275, 0, 4),
	(10, 89, 0),
	(225, 89, 0),
	(0, -45, 1),
	(0, -46, 5),
	(99, -89, 5)
]

@pytest.mark.parametrize("ra, dec, facenum", expected_results)
def test_ang2ipix(ra, dec, facenum):
	'''
	Test Q3C ang2ipix.
	'''
	q3c = Q3C()
	assert facenum == q3c.face_number(ra, dec)
