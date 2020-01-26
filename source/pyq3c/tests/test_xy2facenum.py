
import pytest

# See for C to assrt_approx_equal:
# https://stackoverflow.com/questions/16839658/printf-width-specifier-to-maintain-precision-of-floating-point-value

import pyq3c

# x, y, facenum_in, facenum
expected_results = [
	(0.0000, 0.0000, 5, 5),
	(0.0000, 0.0000, 4, 4),
	(-0.5, 0.5, 5, 5),
	(-0.9, 0.9, 1, 1),
	(0.0, -1.9, 1, 5),
	(1.01, 1.0, 0, 2)
]

@pytest.mark.parametrize("x, y, facenum_in, facenum", expected_results)
def test_xy2facenum(x, y, facenum_in, facenum):
	'''
	Test Q3C xy2facenum.
	'''
	assert facenum == pyq3c.xy2facenum(x, y, facenum_in)
