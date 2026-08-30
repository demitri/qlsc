
import pytest
import numpy as np
from numpy.testing import assert_allclose #, assert_approx_equal

# See for C to assert_approx_equal:
# https://stackoverflow.com/questions/16839658/printf-width-specifier-to-maintain-precision-of-floating-point-value

from qlsc import QLSC

# depth, ipix, [[ra1,dec1], ... , [ran,decn]]

# bin level 0
qlsc0_expected_results = [
	(0, 0, np.array([[315., 35.26438968], [45., 35.26438968], [135., 35.26438968], [225., 35.26438968]])),
	(0, 1, np.array([[315.0, -35.264389682754654], [45.0, -35.264389682754654], [45.0, 35.264389682754654], [315.0, 35.264389682754654]])),
	(0, 2, np.array([[45.0, -35.264389682754654], [135.0, -35.264389682754654], [135.0, 35.264389682754654], [45.0, 35.264389682754654]])),
	(0, 3, np.array([[135.0, -35.264389682754654], [225.0, -35.264389682754654], [225.0, 35.264389682754654], [135.0, 35.264389682754654]])),
	(0, 4, np.array([[225.0, -35.264389682754654], [315.0, -35.264389682754654], [315.0, 35.264389682754654], [225.0, 35.264389682754654]])),
	(0, 5, np.array([[225.0, -35.264389682754654], [135.0, -35.264389682754654], [45.0, -35.264389682754654], [315.0, -35.264389682754654]]))
]

# bin level 2
qlsc2_expected_results = [
	(2, 21, np.array([[26.56505118, -41.8103149], [45., -35.26438968], [ 45., -19.47122063], [ 26.56505118, -24.09484255]])),
	(2, 0, np.array([[315., 35.26438968], [333.43494882, 41.8103149], [315., 54.73561032], [296.56505118, 41.8103149]])),
	(2, 10, np.array([[243.43494882, 41.8103149], [225., 54.73561032], [206.56505118, 41.8103149], [225., 35.26438968]])),
	(2, 60, np.array([[180., 0.], [206.56505118, 0.], [206.56505118, 24.09484255], [180., 26.56505118]]))
]

# bin level 10
qlsc10_expected_results = [
	(10, 1123456, np.array([[336.3706222693432, -38.21416175691502], [336.4646171068093, -38.23410604406816], [336.4646171068093, -38.17074948425823], [336.3706222693432, -38.150815876492366]])),
	(10, 4072375, np.array([[182.57210462788328, 34.25305856945815], [182.68377515946898, 34.25067408080925], [182.68377515946898, 34.32697941232311], [182.57210462788328, 34.32936639388736]]))
]

# bin level 30 (max)
qlsc30_expected_results = [
	(30, 6437157808037343522, np.array([[303.99999460974766, -89.00000003520915], [303.9999980286991, -89.00000012365851], [304.0000030975036, -89.00000006399867], [303.99999967855194, -88.99999997554931]])),
	(30, 2017612633061982208, np.array([[0.0, 0.0], [1.0672170578144923e-07, 0.0], [1.0672170578144923e-07, 1.0672170578144923e-07], [0.0, 1.0672170578144923e-07]])),
	(30, 6629364818475238354, np.array([[43.99999803104112, -89.00000005741951], [44.00000242914403, -88.99999998330696], [43.9999981819457, -88.99999990656117], [43.99999378384282, -88.99999998067372]]))
	
]

expected_results = qlsc0_expected_results + qlsc2_expected_results + qlsc10_expected_results + qlsc30_expected_results

@pytest.mark.parametrize("depth, ipix, polygon", expected_results)
def test_ipix2polygon(depth, ipix, polygon):
	'''
	Test QLSC ipix2polygon.
	'''
	q = QLSC(depth=depth)
	assert_allclose(polygon, q.ipix2polygon(ipix))

@pytest.mark.parametrize("depth, ipix, polygon", expected_results)
def test_ipix2polygon_duplicate_endpoint(depth, ipix, polygon):
	'''
	Test that duplicate_endpoint closes the polygon.

	It used to do nothing at all: the branch built the five row form and the line after
	it stripped the extra row back off unconditionally, so both forms returned four
	rows and the parameter was silently ignored.
	'''
	q = QLSC(depth=depth)

	open_form = q.ipix2polygon(ipix)
	closed_form = q.ipix2polygon(ipix, duplicate_endpoint=True)

	assert open_form.shape == (4,2)
	assert closed_form.shape == (5,2)
	assert_allclose(closed_form[:-1], open_form)
	assert_allclose(closed_form[-1], closed_form[0])

def test_ipix2polygon_shapes_at_every_depth():
	'''
	Test the shape contract of both forms across depths, independently of the
	coordinate values above.
	'''
	for depth in (0, 1, 5, 17, 30):
		q = QLSC(depth=depth)
		ipix = q.ang2ipix(123.4, -45.6)
		assert q.ipix2polygon(ipix).shape == (4,2)
		assert q.ipix2polygon(ipix, duplicate_endpoint=False).shape == (4,2)
		assert q.ipix2polygon(ipix, duplicate_endpoint=True).shape == (5,2)

def test_a_closed_polygon_is_accepted_by_the_polygon_functions():
	'''
	Test that the closed form can still be handed to the containment test and the index
	query, as the docstring promises: the repeated vertex makes a zero-length edge,
	which the crossing test never counts.
	'''
	from qlsc import QLSCIndex

	q = QLSC(depth=10)
	ipix = q.ang2ipix(41., 17.)
	open_form = q.ipix2polygon(ipix)
	closed_form = q.ipix2polygon(ipix, duplicate_endpoint=True)
	centre_ra, centre_dec = q.ipix2ang_center(ipix)

	assert q.point_in_polygon(ra=centre_ra, dec=centre_dec, polygon=closed_form) is True
	assert q.point_in_polygon(ra=centre_ra, dec=centre_dec, polygon=open_form) is True

	index = QLSCIndex(qlsc=q)
	index.add_point(ra=centre_ra, dec=centre_dec, key="centre")
	index.add_point(ra=180., dec=0., key="elsewhere")
	assert set(index.polygon_query(polygon=closed_form, return_key=True).key) == {"centre"}
	assert set(index.polygon_query(polygon=open_form, return_key=True).key) == {"centre"}
