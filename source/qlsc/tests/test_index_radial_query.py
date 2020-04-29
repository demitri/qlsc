
import pytest
import numpy as np
from numpy.testing import assert_approx_equal # use foe scalars
from numpy.testing import assert_allclose     # use for arrays

from qlsc import QLSC, QLSCIndex
from qlsc.generate import sunflower_points_on_sphere

def test_radial_query():
	'''
	Test a radial query.
	'''

	q30 = QLSC(depth=30)
	idx = QLSCIndex(qlsc=q30)
	
	sample_points = sunflower_points_on_sphere(n=5e5)
	
	idx.add_points(points=sample_points)
	
	matches = idx.radial_query(ra=12, dec=-33, radius=2000/3600)
	
	expected_matches = np.array([[ 11.9470522717, -32.8642186861],
                                 [ 12.0478643726, -33.301040754 ],
                                 [ 12.2109817407, -33.0308141307],
                                 [ 11.5200054079, -32.9671430623],
                                 [ 12.4749112364, -33.1977251167],
                                 [ 12.1101696666, -32.5953220055],
                                 [ 11.7839348769, -33.1339331435],
                                 [ 11.6831228027, -32.6979355995],
                                 [ 12.6380286046, -32.9278156427],
                                 [ 11.6208174821, -33.4044789167],
                                 [ 12.3117938416, -33.4684691438],
                                 [ 12.3740991356, -32.7614136187]])

	
	assert_allclose(matches, expected_matches, err_msg="expected matches not found")

