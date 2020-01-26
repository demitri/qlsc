
import pytest

import pyq3c

# ra, dec, depth, ipix
expected_results = [
	(0, 0, 1, 2017612633061982208),
	(0, 0, 2, 2017612633061982211),
	(12, 34, 1, 2178706836447141476)
]

# ra, dec, depth, ipix
depth_out_of_range_expected_results = [
	(12, 34, 0, 6790392854874529380) # the ipix value acutally return in the Q3C PostgreSQL plugin (?)
]

@pytest.mark.parametrize("ra, dec, depth, ipix", expected_results)
def test_ipixcenter(ra, dec, depth, ipix):
	'''
	Test Q3C ipixcenter.
	'''
	q3c = pyq3c.Q3C()
	assert ipix == q3c.ipixcenter(ra, dec, depth)

@pytest.mark.parametrize("ra, dec, depth, ipix", depth_out_of_range_expected_results)
def test_ipixcenter_depth_range(ra, dec, depth, ipix):
	'''
	Test Q3C ipixcenter depth in correct range.
	'''
	with pytest.raises(ValueError):
		q3c = pyq3c.Q3C()
		q3c.ipixcenter(ra, dec, depth) # depth should be int > 0
