
import pytest
from numpy.testing import assert_approx_equal

from qlsc import QLSC

# ipix, depth, facenum
expected_results = [
	(2882303761517117440, 30, 2.0943951024e+00), # ra=45, dec=0
	(2017612633061982208, 30, 2.0943951024e+00), # ra=0, dec=0
	(2882303761517117440,  1, 4.9065384358e-18), # ra=45, dec=0
	(2882303761517117440, 10, 1.2862215819e-12), # ra=45, dec=0
	(480243908601203291,   1, 1.3871447804e-17), # ra=10, dec=89
	(6227361032241724404,  1, 1.0580610731e-17), # ra=134, dec=-66
	(6227361032241724404,  4, 6.7715908182e-16), # ra=134, dec=-66
	(6227361032241724404, 11, 1.1094554873e-11), # ra=134, dec=-66
]

@pytest.mark.parametrize("ipix, depth, area", expected_results)
def test_pixarea(ipix, depth, area):
	'''
	Test QLSC pixarea.
	'''
	q = QLSC()
	assert_approx_equal(area, q.pixarea(ipix, depth))
