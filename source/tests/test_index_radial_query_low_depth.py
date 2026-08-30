
import pytest
import numpy as np

from qlsc import QLSC, QLSCIndex
from qlsc.generate import sunflower_points_on_sphere

# Regression test for radial queries on indices built at low depths.
# The underlying q3c_radial_query produces degenerate (empty) ipix ranges for
# some queries at small nside values (upstream Q3C only ever runs at depth 30),
# silently returning zero matches; QLSCIndex now performs the pixel query at
# depth 30 and maps the ranges down to the index's resolution.

sample_points = sunflower_points_on_sphere(n=20_000)

def brute_force_query(ra0:float, dec0:float, radius:float) -> np.ndarray:
	'''
	Return the sample points within `radius` degrees of (ra0,dec0),
	computed directly with the haversine formula.
	'''
	ra, dec = np.deg2rad(sample_points[:,0]), np.deg2rad(sample_points[:,1])
	r0, d0 = np.deg2rad(ra0), np.deg2rad(dec0)
	h = np.sin((dec-d0)/2)**2 + np.cos(dec)*np.cos(d0)*np.sin((ra-r0)/2)**2
	return sample_points[2*np.arcsin(np.sqrt(h)) <= np.deg2rad(radius)]

# depth, ra, dec, radius; the depth=2 and depth=4 cases produced zero matches
# prior to the depth-30 query mapping
query_cases = [
	(0, 0., 0., 34.9),
	(1, 315., 43., 20.),
	(2, 0., 0., 34.9),
	(4, 0., 90., 10.),
	(5, 12., -33., 5.)
]

@pytest.mark.parametrize("depth, ra, dec, radius", query_cases)
def test_radial_query_low_depth(depth, ra, dec, radius):
	'''
	Test radial queries on low-depth indices against a brute-force calculation.
	'''
	idx = QLSCIndex(qlsc=QLSC(depth=depth))
	idx.add_points(points=sample_points)

	matches = idx.radial_query(ra=ra, dec=dec, radius=radius)
	expected_matches = brute_force_query(ra, dec, radius)

	matches = matches[np.lexsort((matches[:,1], matches[:,0]))]
	expected_matches = expected_matches[np.lexsort((expected_matches[:,1], expected_matches[:,0]))]

	assert matches.shape == expected_matches.shape, "radial query returned a different number of matches than a brute-force calculation"
	np.testing.assert_allclose(matches, expected_matches, err_msg="radial query matches do not agree with a brute-force calculation")
