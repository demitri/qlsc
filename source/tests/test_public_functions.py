
import math

import pytest
from numpy.testing import assert_approx_equal

# These test the module-level convenience wrappers exported from the package
# (the existing distance/sindist tests call the C-level qlsc.q3c functions
# directly, so they would not catch a broken wrapper).
from qlsc import distance, sindist, xy2facenum

# ra1, dec1, ra2, dec2, distance - all degrees
expected_distances = [
	(0.0000, 0.0000, 30.0000, 0.0000, 30.0),
	(15.0000, 88.0000, 15.0000, -88.0000, 176.0),
	(45.4500, 48.8000, 127.8900, -27.7000, 105.844400629850568407164)
]

@pytest.mark.parametrize("ra1, dec1, ra2, dec2, d", expected_distances)
def test_distance(ra1, dec1, ra2, dec2, d):
	'''
	Test the qlsc.distance wrapper.
	'''
	assert_approx_equal(d, distance(ra1, dec1, ra2, dec2))

@pytest.mark.parametrize("ra1, dec1, ra2, dec2, d", expected_distances)
def test_sindist(ra1, dec1, ra2, dec2, d):
	'''
	Test the qlsc.sindist wrapper against the sin^2(d/2) identity.
	'''
	assert_approx_equal(math.sin(math.radians(d)/2)**2, sindist(ra1, dec1, ra2, dec2))

# x, y, reference facenum, expected facenum
# (x,y) within [-1,1] stays on the reference face; beyond it crosses onto
# the neighboring face (see the face table in the README).
expected_faces = [
	(0.0, 0.0, 1, 1),
	(0.0, 1.5, 1, 0),   # above face 1 -> top face
	(0.0, -1.5, 1, 5),  # below face 1 -> bottom face
	(1.5, 0.0, 1, 2),   # +x from face 1 -> face 2
	(-1.5, 0.0, 1, 4)   # -x from face 1 -> face 4
]

@pytest.mark.parametrize("x, y, facenum, expected", expected_faces)
def test_xy2facenum(x, y, facenum, expected):
	'''
	Test the qlsc.xy2facenum wrapper.
	'''
	assert expected == xy2facenum(x, y, facenum)

@pytest.mark.parametrize("facenum", [-1, 6, 99])
def test_xy2facenum_invalid_face(facenum):
	'''
	Test that out-of-range face numbers raise a ValueError.
	'''
	with pytest.raises(ValueError):
		xy2facenum(0., 0., facenum)

def test_embedded_q3c_version():
	'''
	Test that the embedded Q3C version is observable and pinned; a mismatch here
	means a Q3C sync forgot to update Q3C_VERSION in setup.py (or vice versa).
	'''
	from qlsc import q3c
	assert q3c.version() == "2.0.5"

def test_range_overflow_error_type():
	'''
	Test that the dedicated overflow exception type exists and subclasses RuntimeError.
	'''
	from qlsc import q3c
	assert issubclass(q3c.RangeOverflowError, RuntimeError)
