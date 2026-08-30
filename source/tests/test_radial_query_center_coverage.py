
import pytest
import numpy as np

from qlsc import QLSC
from qlsc import q3c

# Regression test for the res_depth patch to the vendored q3cube.c (upstream q3c
# issue #51). res_depth counts classification passes, of which only the first
# res_depth - 1 subdivide; a plain log2(nside/n0) skips the final classification
# pass, and at nside/n0 == 1 skips classification entirely, reading the squares'
# status fields uninitialized - tiny queries then silently miss their own center
# (100% of the time at depth 0). The source-invariant test cannot catch this;
# only behavior can.

@pytest.mark.parametrize("depth", [0, 1, 2, 5, 30])
def test_tiny_radius_queries_cover_their_center(depth):
	'''
	Test that direct C-level tiny-radius queries always cover the center's pixel.
	'''
	rng = np.random.default_rng(depth)
	q = QLSC(depth=depth)
	for _ in range(400):
		ra = float(rng.uniform(0, 360))
		dec = float(np.degrees(np.arcsin(rng.uniform(-1, 1)))) # uniform on the sphere
		center = q.ang2ipix(ra, dec)
		fulls, partials = q3c.radial_query(q._hprm, ra, dec, 1e-7)
		ranges = np.vstack((fulls, partials))
		assert any(a <= center < b for a, b in ranges), \
			f"depth={depth} ({ra},{dec}): center pixel {center} not covered by any candidate range"
