
import math

import pytest
import numpy as np

from qlsc import QLSC
from qlsc import q3c

# QLSC.ipix_neighbors returns the pixels that touch a given pixel: eight in the
# interior of a cube face, seven for the 24 pixels sitting on a corner of the cube
# itself, and four at depth 0 where a face is a single pixel.


def _ipix_at(q, facenum, x, y):
	'''
	Return the ipix of the pixel containing the given point on the given cube face.
	'''
	ra, dec = q3c.xy2ang(x, y, facenum)
	return int(q3c.ang2ipix(q._hprm, ra, dec))


def _cube_corner_pixels(q):
	'''
	Return the ipix values of the pixels that sit on a corner of the cube: the four
	corner pixels of each of the six faces, which meet three at a time at the eight
	cube vertices.
	'''
	d = 2. / q.nside
	corners = set()
	for facenum in range(6):
		for sx in (-1, 1):
			for sy in (-1, 1):
				corners.add(_ipix_at(q, facenum, sx * (1 - d/2), sy * (1 - d/2)))
	return sorted(corners)


def _corner_vectors(q, ipix):
	'''
	Return the four corners of a pixel as 3D unit vectors, so that corners computed
	on two different cube faces can be compared.
	'''
	polygon = q.ipix2polygon(int(ipix))
	ra = np.deg2rad(polygon[:,0])
	dec = np.deg2rad(polygon[:,1])
	return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], axis=1)


def _pixels_touch(q, a, b):
	'''
	Return True if two pixels share at least one corner, which for a grid of
	equal-depth pixels is exactly the condition for them to touch.

	The tolerance is a thousandth of a pixel: two corners are either the same point
	(agreeing to ~1e-16, even when computed on different faces) or a whole pixel apart.
	'''
	tolerance = (math.pi / 2. / q.nside) * 1e-3
	separations = np.linalg.norm(_corner_vectors(q, a)[:,None,:] - _corner_vectors(q, b)[None,:,:], axis=2)
	return separations.min() < tolerance


def _neighbors_by_dense_probe(q, ipix):
	'''
	Find the neighbors of a pixel independently of ipix_neighbors: sample a dense grid
	of points covering the pixel and everything around it, collect the pixels those
	points land in, and keep the ones that share a corner with the target pixel.
	'''
	d = 2. / q.nside
	facenum, x, y = q.ipix2xy(int(ipix))
	x_center, y_center = x + d/2., y + d/2.

	candidates = set()
	for a in np.linspace(-1.7, 1.7, 41):
		for b in np.linspace(-1.7, 1.7, 41):
			candidates.add(_ipix_at(q, facenum, x_center + a*d, y_center + b*d))
	candidates.discard(int(ipix))

	return {c for c in candidates if _pixels_touch(q, ipix, c)}


# ---------------------------------------------------------------- basic properties

@pytest.mark.parametrize("depth", [1, 2, 5, 12, 26, 30])
def test_result_is_sorted_unique_and_in_range(depth):
	'''
	Test the shape of the result: a sorted array of unique, valid ipix values that
	never includes the pixel itself.
	'''
	q = QLSC(depth=depth)
	rng = np.random.default_rng(seed=depth)

	for ipix in rng.integers(0, q.nbins, size=40):
		ipix = int(ipix)
		neighbors = q.ipix_neighbors(ipix)

		assert isinstance(neighbors, np.ndarray)
		assert np.issubdtype(neighbors.dtype, np.integer)
		assert np.array_equal(neighbors, np.sort(neighbors))
		assert len(set(neighbors.tolist())) == len(neighbors)
		assert ipix not in neighbors
		assert ((0 <= neighbors) & (neighbors < q.nbins)).all()


@pytest.mark.parametrize("depth", [5, 12, 26, 30])
def test_interior_pixels_have_eight_neighbors(depth):
	'''
	Test that a pixel well inside a cube face has exactly eight neighbors, and that
	they all live on the same face.
	'''
	q = QLSC(depth=depth)
	bins_per_face = q.nside * q.nside

	for facenum in range(6):
		for x, y in [(0.1, 0.1), (-0.4, 0.6), (0.75, -0.25)]:
			ipix = _ipix_at(q, facenum, x, y)
			neighbors = q.ipix_neighbors(ipix)

			assert len(neighbors) == 8, f"face {facenum} at ({x},{y}) has {len(neighbors)} neighbors"
			assert set((neighbors // bins_per_face).tolist()) == {facenum}


def test_scalar_like_input():
	'''
	Test that NumPy integers and 0-d arrays are accepted, as elsewhere in the API.
	'''
	q = QLSC(depth=5)
	expected = q.ipix_neighbors(100)

	for form in (np.int64(100), np.int32(100), np.array(100)):
		assert np.array_equal(q.ipix_neighbors(form), expected)


def test_invalid_input():
	'''
	Test that bad input raises instead of returning a plausible-looking answer.
	'''
	q = QLSC(depth=5)

	for bad in (-1, q.nbins, q.nbins + 1000):
		with pytest.raises(ValueError):
			q.ipix_neighbors(bad)

	for bad in (1.5, 10.0, "10", None):
		with pytest.raises((ValueError, TypeError)):
			q.ipix_neighbors(bad)

	# this method is not vectorized; an array must be refused, not silently
	# reduced to its first element
	with pytest.raises(ValueError):
		q.ipix_neighbors(np.array([1, 2]))
	with pytest.raises(ValueError):
		q.ipix_neighbors([1, 2])


# ---------------------------------------------------------------- face boundaries

@pytest.mark.parametrize("depth", [2, 5, 26, 30])
def test_cube_corner_pixels_have_seven_neighbors(depth):
	'''
	Test the documented special case: the pixels on a corner of the cube have seven
	neighbors, not eight, because three faces (and therefore three pixels) meet at a
	cube vertex rather than four.
	'''
	q = QLSC(depth=depth)
	bins_per_face = q.nside * q.nside
	corners = _cube_corner_pixels(q)

	assert len(corners) == 24, "there are 8 cube vertices, each shared by 3 face corner pixels"

	for ipix in corners:
		neighbors = q.ipix_neighbors(ipix)
		assert len(neighbors) == 7, f"cube corner pixel {ipix} at depth {depth} has {len(neighbors)} neighbors"
		# a corner pixel sees its own face and the two others meeting at the vertex
		faces = set((neighbors // bins_per_face).tolist())
		assert len(faces) == 3


@pytest.mark.parametrize("depth", [5, 26, 30])
def test_face_edge_pixels_have_eight_neighbors(depth):
	'''
	Test pixels on each of the four edges of every face (away from the cube corners):
	they keep all eight neighbors, spread over two faces.
	'''
	q = QLSC(depth=depth)
	d = 2. / q.nside
	bins_per_face = q.nside * q.nside

	# positions along an edge: next to a corner, two pixels in, and mid-edge
	along = [-1 + 1.5*d, -1 + 2.5*d, 0.5*d, 1 - 2.5*d, 1 - 1.5*d]

	for facenum in range(6):
		for t in along:
			for x, y in [(-1 + d/2, t), (1 - d/2, t), (t, -1 + d/2), (t, 1 - d/2)]:
				ipix = _ipix_at(q, facenum, x, y)
				neighbors = q.ipix_neighbors(ipix)

				assert len(neighbors) == 8, f"edge pixel {ipix} (face {facenum}, x={x}, y={y}) has {len(neighbors)} neighbors"
				faces = set((neighbors // bins_per_face).tolist())
				assert len(faces) == 2, f"an edge pixel should see exactly two faces; saw {faces}"


@pytest.mark.parametrize("facenum", [0, 5])
def test_pole_pixels(facenum):
	'''
	Test the four pixels around the pole at the center of the top and bottom faces,
	where the face's (x,y) origin sits.
	'''
	q = QLSC(depth=20)
	d = 2. / q.nside

	for sx in (-1, 1):
		for sy in (-1, 1):
			ipix = _ipix_at(q, facenum, sx * d/2, sy * d/2)
			neighbors = q.ipix_neighbors(ipix)
			assert len(neighbors) == 8
			for neighbor in neighbors:
				assert ipix in q.ipix_neighbors(int(neighbor))


def test_depth_zero():
	'''
	Test the degenerate scheme where each face is a single pixel: a face's neighbors
	are the four faces adjacent to it, and never the opposite face.
	'''
	q = QLSC(depth=0)
	opposite = {0: 5, 5: 0, 1: 3, 3: 1, 2: 4, 4: 2}

	for facenum in range(6):
		neighbors = q.ipix_neighbors(facenum)
		assert len(neighbors) == 4, f"face {facenum} has {len(neighbors)} neighbors"
		assert facenum not in neighbors
		assert opposite[facenum] not in neighbors, "the opposite face does not touch this one"


# ---------------------------------------------------------------- symmetry

@pytest.mark.parametrize("depth", [0, 1, 2, 3])
def test_symmetry_exhaustive(depth):
	'''
	Test that neighborhood is symmetric across every pixel of the scheme: if B is a
	neighbor of A then A is a neighbor of B. At these depths every pixel touches a
	face boundary, so this exercises the boundary handling exhaustively.
	'''
	q = QLSC(depth=depth)
	table = {ipix: set(q.ipix_neighbors(ipix).tolist()) for ipix in range(q.nbins)}

	for ipix, neighbors in table.items():
		for neighbor in neighbors:
			assert ipix in table[neighbor], f"depth {depth}: {neighbor} is a neighbor of {ipix} but not the reverse"


@pytest.mark.parametrize("depth", [5, 12, 26, 30])
def test_symmetry_random(depth):
	'''
	Test the same symmetry on random pixels at higher depths.
	'''
	q = QLSC(depth=depth)
	rng = np.random.default_rng(seed=100 + depth)

	for ipix in rng.integers(0, q.nbins, size=50):
		ipix = int(ipix)
		for neighbor in q.ipix_neighbors(ipix):
			assert ipix in q.ipix_neighbors(int(neighbor)), f"depth {depth}: {ipix} -> {neighbor} is not symmetric"


@pytest.mark.parametrize("depth", [2, 5, 26, 30])
def test_symmetry_at_cube_corners(depth):
	'''
	Test symmetry specifically at the cube corners, where a pixel has seven neighbors
	while most of those neighbors have eight. The relation must still hold in both
	directions - this is the case worth checking rather than assuming, since the
	neighbor across the vertex is reached by a different step at each end.
	'''
	q = QLSC(depth=depth)

	for ipix in _cube_corner_pixels(q):
		for neighbor in q.ipix_neighbors(ipix):
			back = q.ipix_neighbors(int(neighbor))
			assert ipix in back, f"depth {depth}: cube corner {ipix} -> {neighbor} is not symmetric"


# ---------------------------------------------------------------- geometry

@pytest.mark.parametrize("depth", [2, 5, 12, 26, 30])
def test_neighbors_share_a_corner(depth):
	'''
	Test the defining property geometrically and independently of how the neighbors
	were found: every reported neighbor shares at least one corner point with the
	pixel, comparing corners as 3D unit vectors so that pixels on different faces
	can be compared.
	'''
	q = QLSC(depth=depth)
	rng = np.random.default_rng(seed=200 + depth)

	sample = [int(v) for v in rng.integers(0, q.nbins, size=12)] + _cube_corner_pixels(q)[:6]

	for ipix in sample:
		for neighbor in q.ipix_neighbors(ipix):
			assert _pixels_touch(q, ipix, int(neighbor)), \
				f"depth {depth}: pixel {neighbor} does not touch {ipix}"


@pytest.mark.parametrize("depth", [2, 5, 12, 26, 30])
def test_no_neighbor_is_missed(depth):
	'''
	Test completeness against an independent method: densely sample the neighborhood
	of a pixel in the face plane, keep every pixel that shares a corner with it, and
	require that set to be exactly what ipix_neighbors returns.
	'''
	q = QLSC(depth=depth)
	rng = np.random.default_rng(seed=300 + depth)

	sample = [int(v) for v in rng.integers(0, q.nbins, size=4)] + _cube_corner_pixels(q)[:3]

	for ipix in sample:
		expected = _neighbors_by_dense_probe(q, ipix)
		assert set(q.ipix_neighbors(ipix).tolist()) == expected, \
			f"depth {depth}: neighbors of {ipix} disagree with the dense probe"


@pytest.mark.parametrize("depth", [12, 26, 30])
def test_neighbor_centers_are_close(depth):
	'''
	Test a coarse sanity check on the returned pixels: the center of every neighbor
	is within two pixel diagonals of the center of the pixel.
	'''
	from qlsc import distance

	q = QLSC(depth=depth)
	rng = np.random.default_rng(seed=400 + depth)
	# the angular size of a pixel: a face spans 90 degrees across nside pixels
	pixel_degrees = 90. / q.nside

	for ipix in rng.integers(0, q.nbins, size=20):
		ipix = int(ipix)
		ra, dec = q.ipix2ang(ipix)
		neighbor_ra, neighbor_dec = q.ipix2ang(q.ipix_neighbors(ipix))

		for n_ra, n_dec in zip(neighbor_ra, neighbor_dec):
			assert distance(ra, dec, n_ra, n_dec) < 2 * math.sqrt(2) * pixel_degrees


def test_region_growing():
	'''
	Test the use case the method exists for: growing a region outward from a pixel.
	The first ring is the neighbors themselves; the second ring adds the neighbors of
	those, and away from any cube corner a 5x5 block of 25 pixels results.
	'''
	q = QLSC(depth=15)
	ipix = _ipix_at(q, 3, 0.13, -0.42) # comfortably inside a face

	region = {ipix} | set(q.ipix_neighbors(ipix).tolist())
	assert len(region) == 9

	grown = set(region)
	for member in region:
		grown.update(q.ipix_neighbors(member).tolist())
	assert len(grown) == 25
