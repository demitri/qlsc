
import pytest

from pyq3c import Q3C

# bin_level, ra, dec, ipix
expected_results = [
	(0, 0, 0, 1),
	(0, -22, 0, 1),
	(0, -45, 0, 1),
	(0, 27, 0, 1),
	(0, 45, 0, 2),
	(0, 46, 0, 2),
	(0, 256, 0, 4),
	(0, 12, 34, 1),
	(0, 12, 87, 0),
	(0, 12, -73, 5),
	(0, 12, -3, 1),
	(30, 0, 0, 2017612633061982208),
	(30, 180, 0, 4323455642275676160),
	(30, -1, 0, 1825388318008017157),
	(30, 359, 0, 1825388318008017157),
	(30, 0, -90, 6629298651489370112),
	(30, 57.3, 88, 480333897118235422)
]

@pytest.mark.parametrize("bin_level, ra, dec, ipix", expected_results)
def test_ang2ipix(bin_level, ra, dec, ipix):
	'''
	Test Q3C ang2ipix.
	'''
	qlsc = Q3C(bin_level=bin_level)
	assert ipix == qlsc.ang2ipix(ra, dec)
