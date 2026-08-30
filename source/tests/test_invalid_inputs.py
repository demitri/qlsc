
import pytest
import numpy as np

from qlsc import QLSC
from qlsc import q3c

q30 = QLSC(depth=30)
poly_ra = np.array([5., 5., -5.])
poly_dec = np.array([-5., 5., 5.])

# every C function that dereferences an hprm capsule, with representative arguments
# (q3c.facenum accepts an hprm argument but never uses it, so it is not listed)
capsule_functions = [
	lambda c: q3c.nside(c),
	lambda c: q3c.ang2ipix(c, 10., 20.),
	lambda c: q3c.ang2ipix_xy(c, 10., 20.),
	lambda c: q3c.ipix2ang(c, 42),
	lambda c: q3c.ipix2xy(c, 42),
	lambda c: q3c.pixarea(c, 42, 1),
	lambda c: q3c.radial_query(c, 10., 20., 1.),
	lambda c: q3c.check_sphere_point_in_poly(c, 0., 0., poly_ra, poly_dec)
]

@pytest.mark.parametrize("call", capsule_functions)
@pytest.mark.parametrize("bad_capsule", [None, 42, "not a capsule"])
def test_invalid_hprm_capsule_raises(call, bad_capsule):
	'''
	Test that C-level functions raise an exception (rather than crash the
	interpreter) when handed something that is not an hprm capsule.
	'''
	# PyCapsule_GetPointer raises ValueError for a non-capsule object; TypeError
	# can come from argument parsing. Anything broader could mask a wrong failure.
	with pytest.raises((ValueError, TypeError)):
		call(bad_capsule)

def test_point_in_polygon_non_finite():
	'''
	Test that non-finite inputs to point_in_polygon raise rather than
	silently returning False.
	'''
	poly = np.array([[5.,-5.],[5.,5.],[-5.,5.]])
	with pytest.raises(ValueError):
		q30.point_in_polygon(ra=np.nan, dec=0., polygon=poly)
	with pytest.raises(ValueError):
		q30.point_in_polygon(ra=0., dec=np.nan, polygon=poly)
	with pytest.raises(ValueError):
		q30.point_in_polygon(ra=0., dec=0., polygon=np.array([[np.nan,-5.],[5.,5.],[-5.,5.]]))
	with pytest.raises(ValueError):
		q30.point_in_polygon(ra=0., dec=0., polygon=np.array([[np.inf,-5.],[5.,5.],[-5.,5.]]))
