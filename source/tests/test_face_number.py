
import pytest
import numpy as np

from qlsc import QLSC, q3c

# Tests for the two ways of asking which cube face something is on:
# QLSC.ipix2face(ipix) and QLSC.face_number(ra, dec) / face_number(ipix=...).
#
# Both were broken. ipix2face divided through a float, which at depths above about 26
# put the last pixel of every face on the next face - at depth 30 the last pixel of the
# scheme came back as face 6, which is not a face. face_number's ipix branch could not
# be reached at all: ra and dec were required positional parameters, so face_number(ipix=n)
# raised TypeError before any of its logic ran, and the branch was guarded by "if ipix:",
# which is false for face 0's first pixel.

depths = [0, 1, 2, 5, 17, 26, 27, 30]

@pytest.mark.parametrize("depth", depths)
def test_ipix2face_returns_an_int_in_range(depth):
	'''
	Test that the face number is an integer in [0,5], not a float.
	'''
	qlsc = QLSC(depth=depth)
	for ipix in (0, qlsc.nbins // 2, qlsc.nbins - 1):
		face = qlsc.ipix2face(ipix)
		assert isinstance(face, int), f"ipix2face returned {type(face).__name__}, not int"
		assert 0 <= face <= 5

@pytest.mark.parametrize("depth", depths)
def test_ipix2face_at_every_face_boundary(depth):
	'''
	Test the first, second and LAST pixel of every face.

	The last pixel of a face is the one that went wrong: ipix // math.pow(nside,2)
	rounds the numerator up to the next multiple of the divisor once the values pass
	2^53, so ipix2face reported the pixel one past the end of each face.
	'''
	qlsc = QLSC(depth=depth)
	pixels_per_face = qlsc.nside * qlsc.nside

	for face in range(6):
		first = face * pixels_per_face
		last = first + pixels_per_face - 1
		assert qlsc.ipix2face(first) == face
		assert qlsc.ipix2face(last) == face, \
			f"depth={depth}: the last pixel of face {face} was placed on face {qlsc.ipix2face(last)}"
		if pixels_per_face > 1:
			assert qlsc.ipix2face(first + 1) == face

def test_ipix2face_never_reports_a_face_that_does_not_exist():
	'''
	Test the specific symptom at depth 30: the very last pixel of the scheme was
	reported as being on face 6.
	'''
	qlsc = QLSC(depth=30)
	assert qlsc.ipix2face(qlsc.nbins - 1) == 5

@pytest.mark.parametrize("depth", depths)
def test_ipix2face_agrees_with_the_coordinate_of_the_pixel(depth):
	'''
	Test the integer arithmetic against the independent geometric answer: the face
	q3c computes from the pixel's own centre coordinate.
	'''
	qlsc = QLSC(depth=depth)
	pixels_per_face = qlsc.nside * qlsc.nside
	rng = np.random.default_rng(depth)

	sample = [0, qlsc.nbins - 1]
	sample += [int(value) for value in rng.integers(0, qlsc.nbins, 60)]
	for face in range(6):
		sample += [face * pixels_per_face, face * pixels_per_face + pixels_per_face - 1]

	for ipix in sample:
		ra, dec = qlsc.ipix2ang(ipix)
		assert qlsc.ipix2face(ipix) == q3c.facenum(qlsc._hprm, ra, dec), \
			f"depth={depth}: ipix {ipix} placed differently by the pixel number and by its coordinate"

# ipix
invalid_ipix_values = [-1, 1.5, "3", None]

@pytest.mark.parametrize("ipix", invalid_ipix_values)
def test_ipix2face_rejects_invalid_values(ipix):
	'''
	Test that a non-integer or out-of-range pixel number raises rather than being
	truncated or clamped into some other face.
	'''
	qlsc = QLSC(depth=5)
	with pytest.raises(ValueError):
		qlsc.ipix2face(ipix)

def test_ipix2face_rejects_a_value_past_the_end_of_the_scheme():
	'''
	Test the upper bound, which is the one a float division was most likely to slip past.
	'''
	qlsc = QLSC(depth=5)
	with pytest.raises(ValueError):
		qlsc.ipix2face(qlsc.nbins)

def test_face_number_accepts_ipix_as_a_keyword():
	'''
	Test that face_number(ipix=...) works at all.

	ra and dec were required positional parameters, so every call of this form raised
	TypeError before reaching the ipix handling, and the only way in was the undocumented
	face_number(None, None, ipix=n).
	'''
	qlsc = QLSC(depth=5)
	pixels_per_face = qlsc.nside * qlsc.nside
	for face in range(6):
		assert qlsc.face_number(ipix=face * pixels_per_face) == face

def test_face_number_accepts_the_first_pixel():
	'''
	Test ipix=0 specifically: the branch was guarded with "if ipix:", so pixel zero -
	the first pixel of face 0 - fell through to the coordinate branch and failed there.
	'''
	qlsc = QLSC(depth=5)
	assert qlsc.face_number(ipix=0) == 0

@pytest.mark.parametrize("depth", depths)
def test_face_number_by_ipix_matches_ipix2face(depth):
	'''
	Test that the two entry points give the same answer, at every face boundary.
	'''
	qlsc = QLSC(depth=depth)
	pixels_per_face = qlsc.nside * qlsc.nside
	for face in range(6):
		for ipix in (face * pixels_per_face, face * pixels_per_face + pixels_per_face - 1):
			assert qlsc.face_number(ipix=ipix) == qlsc.ipix2face(ipix) == face

# ra, dec, expected face
coordinate_cases = [
	(0., 89., 0),
	(0., 0., 1),
	(90., 0., 2),
	(180., 0., 3),
	(270., 0., 4),
	(0., -89., 5)
]

@pytest.mark.parametrize("ra, dec, face", coordinate_cases)
def test_face_number_by_coordinate(ra, dec, face):
	'''
	Test the coordinate form, positionally and by keyword.
	'''
	qlsc = QLSC(depth=5)
	assert qlsc.face_number(ra, dec) == face
	assert qlsc.face_number(ra=ra, dec=dec) == face

def test_face_number_by_coordinate_is_periodic_in_ra():
	'''
	Test that the coordinate form still accepts any right ascension, as it did before.
	'''
	qlsc = QLSC(depth=5)
	for shift in (-720., -360., 0., 360., 720.):
		assert qlsc.face_number(10. + shift, 20.) == qlsc.face_number(10., 20.)

def test_face_number_requires_a_complete_coordinate():
	'''
	Test that half a coordinate raises rather than being silently completed.
	'''
	qlsc = QLSC(depth=5)
	for arguments in ({}, {"ra": 10.}, {"dec": 20.}):
		with pytest.raises(ValueError):
			qlsc.face_number(**arguments)

def test_face_number_refuses_a_coordinate_and_a_pixel_together():
	'''
	Test that giving both forms raises rather than one silently winning.
	'''
	qlsc = QLSC(depth=5)
	with pytest.raises(ValueError):
		qlsc.face_number(10., 20., ipix=0)
	with pytest.raises(ValueError):
		qlsc.face_number(ra=10., dec=20., ipix=100)

def test_face_number_validates_dec():
	'''
	Test that a declination off the sphere raises.
	'''
	qlsc = QLSC(depth=5)
	for dec in (91., -91., np.nan):
		with pytest.raises(ValueError):
			qlsc.face_number(10., dec)

def test_face_number_validates_ipix():
	'''
	Test that an out-of-range or non-integer pixel number raises in this entry point too.
	'''
	qlsc = QLSC(depth=5)
	for ipix in (-1, qlsc.nbins, 1.5):
		with pytest.raises(ValueError):
			qlsc.face_number(ipix=ipix)
