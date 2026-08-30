
import pytest
import numpy as np

from qlsc import QLSC

q30 = QLSC(depth=30)

small_square = np.array([[5.,-5.],[5.,5.],[-5.,5.],[-5.,-5.]])

# ra, dec, expected
small_square_cases = [
	(0., 0., True),
	(4.9, 0., True),
	(0., 4.9, True),
	(20., 0., False),
	(0., 20., False),
	(180., 0., False),
	(0., -90., False)
]

@pytest.mark.parametrize("ra, dec, expected", small_square_cases)
def test_point_in_polygon(ra, dec, expected):
	'''
	Test point_in_polygon with a small square about the origin.
	'''
	assert expected == q30.point_in_polygon(ra=ra, dec=dec, polygon=small_square)

def test_point_in_polygon_face_boundary():
	'''
	Test a polygon that straddles a cube face boundary (ra=45), which exercises
	the multi-face projection path.
	'''
	poly = np.array([[40.,-5.],[50.,-5.],[50.,5.],[40.,5.]])
	assert q30.point_in_polygon(ra=45., dec=0., polygon=poly) is True
	assert q30.point_in_polygon(ra=42., dec=3., polygon=poly) is True
	assert q30.point_in_polygon(ra=45., dec=20., polygon=poly) is False
	assert q30.point_in_polygon(ra=35., dec=0., polygon=poly) is False

def test_point_in_polygon_large():
	'''
	Test a large (60 degrees across) polygon that still fits on one face projection.
	'''
	poly = np.array([[30.,-30.],[30.,30.],[-30.,30.],[-30.,-30.]])
	assert q30.point_in_polygon(ra=0., dec=0., polygon=poly) is True
	assert q30.point_in_polygon(ra=180., dec=0., polygon=poly) is False

def test_point_in_polygon_too_large():
	'''
	Test that oversized polygons raise ValueError rather than returning wrong results.

	These polygons overflow the projected face bounding box on three or four sides,
	which is the geometry that overflowed q3c's undersized internal buffers (the
	points[4] stack buffer patched locally to points[8] - see CLAUDE.md); the
	projection loop runs before the polygon is rejected, so these cases exercise
	the patched code path.
	'''
	# the shape that reproduced the q3c stack overflow under ASan
	huge = np.array([[80.,0.],[0.,80.],[280.,0.],[0.,-80.]])
	with pytest.raises(ValueError):
		q30.point_in_polygon(ra=10., dec=10., polygon=huge)

	# a diamond whose bounding box overflows all four face sides
	diamond = np.array([[46.,0.],[0.,46.],[314.,0.],[0.,-46.]])
	with pytest.raises(ValueError):
		q30.point_in_polygon(ra=0., dec=0., polygon=diamond)

def test_point_in_polygon_pixel_round_trip():
	'''
	Test that a pixel's polygon (from ipix2polygon) contains the pixel's own center.
	'''
	for depth in (2, 10, 30):
		q = QLSC(depth=depth)
		for (ra, dec) in [(12., 34.), (300., -55.), (0., 90.)]:
			ipix = q.ang2ipix(ra, dec)
			poly = np.asarray(q.ipix2polygon(ipix))
			cra, cdec = q.ipix2ang(ipix)
			assert q.point_in_polygon(ra=cra, dec=cdec, polygon=poly) is True, \
				f"pixel {ipix} (depth {depth}) polygon does not contain its own center"

def test_point_in_polygon_validation():
	'''
	Test input validation.
	'''
	with pytest.raises(ValueError):
		q30.point_in_polygon(ra=0., dec=0., polygon=small_square[:2]) # too few vertices
	with pytest.raises(ValueError):
		q30.point_in_polygon(ra=0., dec=0., polygon=np.zeros((101,2))) # too many vertices
	with pytest.raises(ValueError):
		q30.point_in_polygon(ra=0., dec=0., polygon=np.zeros((4,3))) # wrong shape
	with pytest.raises(ValueError):
		q30.point_in_polygon(ra=0., dec=0., polygon=np.array([[0.,95.],[10.,0.],[0.,-10.]])) # vertex dec out of range
	with pytest.raises(ValueError):
		q30.point_in_polygon(ra=0., dec=0.) # missing polygon
