
import pytest
import numpy as np
from numpy.testing import assert_allclose

from qlsc import QLSC, QLSCIndex
from qlsc.generate import sunflower_points_on_sphere

# Regression test for queries with large radii just below 35 degrees,
# which missed objects in some areas of the sky (and could crash) prior
# to the Q3C v2.0.5 fix: https://github.com/segasai/q3c/pull/49

sample_points = sunflower_points_on_sphere(n=50_000)

def brute_force_query(ra0:float, dec0:float, radius:float) -> np.ndarray:
	'''
	Return the sample points within `radius` degrees of (ra0,dec0),
	computed directly with the haversine formula.
	'''
	ra, dec = np.deg2rad(sample_points[:,0]), np.deg2rad(sample_points[:,1])
	r0, d0 = np.deg2rad(ra0), np.deg2rad(dec0)
	h = np.sin((dec-d0)/2)**2 + np.cos(dec)*np.cos(d0)*np.sin((ra-r0)/2)**2
	return sample_points[2*np.arcsin(np.sqrt(h)) <= np.deg2rad(radius)]

# ra, dec, radius (degrees)
query_centers = [
	(0., 0., 34.9),
	(45., 35., 34.9),
	(315., -45., 34.9),
	(180., 89., 34.9),
	(45., 35., 34.0)
]

@pytest.mark.parametrize("ra, dec, radius", query_centers)
def test_radial_query_large_radius(ra, dec, radius):
	'''
	Test radial queries with radii just below 35 degrees against a brute-force calculation.
	'''
	q30 = QLSC(depth=30)
	idx = QLSCIndex(qlsc=q30)
	idx.add_points(points=sample_points)

	matches = idx.radial_query(ra=ra, dec=dec, radius=radius)
	expected_matches = brute_force_query(ra, dec, radius)

	# sort both lexicographically by (ra,dec) so they can be compared row by row
	matches = matches[np.lexsort((matches[:,1], matches[:,0]))]
	expected_matches = expected_matches[np.lexsort((expected_matches[:,1], expected_matches[:,0]))]

	assert matches.shape == expected_matches.shape, "radial query returned a different number of matches than a brute-force calculation"
	assert_allclose(matches, expected_matches, err_msg="radial query matches do not agree with a brute-force calculation")
