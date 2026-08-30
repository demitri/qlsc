
import json
import math

import pytest
import numpy as np

from qlsc import QLSC, QLSCIndex
from qlsc.generate import sunflower_points_on_sphere

# QLSC.coverage and QLSCIndex.coverage emit the "sky-coverage" version 1 interchange
# format specified by astplot-js (AI/COVERAGE_LAYER_SPEC.md in that repository).
#
# The validator below is written FROM THAT SPECIFICATION, not from qlsc's implementation:
# every check names the rule it enforces, and the checks are the ones a consumer performs
# at ingestion (spec sections 1.1, 1.2 and 4.4). If qlsc's idea of the format and the
# spec's ever diverge, these tests are meant to fail.

# spec 1.1: the complete set of version-1 top-level keys. "Unknown top-level keys are an
# error, not ignored."
SPEC_TOP_LEVEL_KEYS = {"type", "version", "polygons", "units", "frame", "equinox", "edges", "value", "meta"}

# spec 1.2: the complete set of keys a polygon entry may carry.
SPEC_POLYGON_KEYS = {"vertices", "id", "value"}

# spec 1.1: frames the format defines, case-insensitive.
SPEC_FRAMES = {"icrs", "fk5", "fk4", "galactic", "ecliptic"}

# spec 3.4: frames that have no equinox, for which declaring one is an error.
SPEC_FRAMES_WITHOUT_EQUINOX = {"icrs", "galactic"}


def _unit_vector(lon, lat, units):
	'''
	Convert one vertex to a unit 3-vector, as the consumer does before interpolating.
	'''
	if units == "deg":
		lon, lat = math.radians(lon), math.radians(lat)
	return (math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat))


def validate_sky_coverage(document):
	'''
	Assert that a document satisfies every MUST of sky-coverage v1 (spec 1.1, 1.2, 4.4).

	This is the consumer's ingestion contract expressed as assertions; a document that
	passes it is one astplot-js's CoverageLayer accepts.
	'''
	# --- spec 1.1: top-level fields -------------------------------------------------
	assert isinstance(document, dict), "spec 1: the document is a single JSON object"

	unknown = set(document) - SPEC_TOP_LEVEL_KEYS
	assert not unknown, f"spec 1.1: unknown top-level keys are an error; found {sorted(unknown)}"

	assert "type" in document, "spec 1.1: 'type' is required"
	assert document["type"] == "sky-coverage", "spec 1.1: 'type' must be exactly \"sky-coverage\""

	assert "version" in document, "spec 1.1: 'version' is required"
	assert isinstance(document["version"], int) and not isinstance(document["version"], bool), \
		"spec 1.1: 'version' is an integer"
	assert document["version"] == 1, "spec 1.1/4.4: a consumer rejects a version it does not know"

	assert "polygons" in document, "spec 1.1: 'polygons' is required"
	assert isinstance(document["polygons"], list), "spec 4.4: 'polygons' must be an array"

	units = document.get("units", "deg")
	assert units in ("deg", "rad"), f"spec 1.1/4.4: unknown 'units' value '{units}'"

	frame = document.get("frame", "icrs")
	assert isinstance(frame, str), "spec 1.1: 'frame' is a string"
	assert frame.lower() in SPEC_FRAMES, f"spec 1.1/4.4: unknown 'frame' value '{frame}'"

	if "equinox" in document:
		assert frame.lower() not in SPEC_FRAMES_WITHOUT_EQUINOX, \
			f"spec 1.1/3.4/4.4: declaring an 'equinox' for '{frame}' is an error"
		assert isinstance(document["equinox"], (int, float)) and math.isfinite(document["equinox"]), \
			"spec 4.4: 'equinox' must be a finite number"

	edges = document.get("edges", "geodesic")
	assert edges == "geodesic", f"spec 1.1: only \"geodesic\" is defined in v1; any other value is an error (got '{edges}')"

	if "value" in document:
		descriptor = document["value"]
		assert isinstance(descriptor, dict), "spec 1.1: 'value' is an object"
		assert not (set(descriptor) - {"name", "units", "domain"}), \
			f"spec 1.1: the value descriptor is {{name, units, domain}}; found {sorted(descriptor)}"
		if "name" in descriptor:
			assert isinstance(descriptor["name"], str), "spec 1.1: value.name is a string"
		if "units" in descriptor:
			assert isinstance(descriptor["units"], str), "spec 1.1: value.units is a string"
		if "domain" in descriptor:
			domain = descriptor["domain"]
			assert isinstance(domain, list) and len(domain) == 2, "spec 1.1: value.domain is [min, max]"
			assert all(isinstance(v, (int, float)) and math.isfinite(v) for v in domain), \
				"spec 1.1: value.domain bounds are numbers"
			assert domain[0] <= domain[1], "spec 1.1: value.domain is [min, max]"

	if "meta" in document:
		assert isinstance(document["meta"], dict), "spec 1.1: 'meta' is an object"

	# --- spec 1.2 / 4.4: polygon entries --------------------------------------------
	for n,entry in enumerate(document["polygons"]):
		assert isinstance(entry, dict), f"spec 4.4: polygon {n} is not an object"

		unknown = set(entry) - SPEC_POLYGON_KEYS
		assert not unknown, f"spec 1.2: unknown keys inside polygon {n} are an error; found {sorted(unknown)}"

		assert "vertices" in entry, f"spec 1.2: polygon {n} has no 'vertices'"
		vertices = entry["vertices"]
		assert isinstance(vertices, list), f"spec 4.4: polygon {n} 'vertices' must be an array"
		assert len(vertices) >= 3, f"spec 1.2/4.4: polygon {n} has fewer than 3 vertices"

		for v,vertex in enumerate(vertices):
			assert isinstance(vertex, (list, tuple)) and len(vertex) == 2, \
				f"spec 4.4: polygon {n} vertex {v} is not a 2-element pair"
			lon, lat = vertex
			for value in (lon, lat):
				assert isinstance(value, (int, float)) and not isinstance(value, bool), \
					f"spec 4.4: polygon {n} vertex {v} is not a pair of numbers"
				assert math.isfinite(value), f"spec 4.4: polygon {n} vertex {v} is not finite"
			latitude = lat if units == "rad" else math.radians(lat)
			assert -math.pi/2 - 1e-12 <= latitude <= math.pi/2 + 1e-12, \
				f"spec 4.4: polygon {n} vertex {v} has a latitude outside [-90,+90]"

		# spec 4.4: two consecutive vertices antipodal to within ~2e-6 rad are an error,
		# because the great-circle arc between them is not unique. The ring is implicitly
		# closed (spec 1.2), so the last-to-first edge is tested too.
		for v in range(len(vertices)):
			a = _unit_vector(vertices[v][0], vertices[v][1], units)
			b = _unit_vector(vertices[(v+1) % len(vertices)][0], vertices[(v+1) % len(vertices)][1], units)
			cosine = sum(x*y for x,y in zip(a,b))
			assert not (cosine < -1 + 1e-12), \
				f"spec 4.4: polygon {n} vertices {v} and {(v+1) % len(vertices)} are antipodal"

		if "id" in entry:
			assert isinstance(entry["id"], (int, float, str)) and not isinstance(entry["id"], bool), \
				f"spec 4.4: polygon {n} 'id' is neither a number nor a string"

		if "value" in entry and entry["value"] is not None:
			value = entry["value"]
			assert isinstance(value, (int, float)) and not isinstance(value, bool), \
				f"spec 4.4: polygon {n} 'value' is present but not a number"
			assert math.isfinite(value), f"spec 4.4: polygon {n} 'value' is not finite"


def assert_json_native(value, path="document"):
	'''
	Assert that every leaf of the document is a Python type the standard library JSON
	encoder emits, checked by exact type rather than isinstance.

	np.float64 is a subclass of float and so would pass an isinstance test and serialize
	silently; np.int64 is not a subclass of int and would raise from json.dumps. Checking
	the exact type catches both, which is the point: the spec's Python example (1.3) is a
	few lines of json.dumps, and it has to work with no custom encoder.
	'''
	if isinstance(value, dict):
		for key,item in value.items():
			assert type(key) is str, f"{path}: key {key!r} is a {type(key).__name__}, not a str"
			assert_json_native(item, f"{path}[{key!r}]")
	elif isinstance(value, list):
		for n,item in enumerate(value):
			assert_json_native(item, f"{path}[{n}]")
	else:
		assert type(value) in (str, int, float, bool, type(None)), \
			f"{path} is a {type(value).__name__} ({value!r}), not a JSON-native Python type"


def _index(points, depth=30, keys=None):
	index = QLSCIndex(qlsc=QLSC(depth=depth))
	points = np.asarray(points, dtype=np.double)
	index.add_points(ra=points[:,0], dec=points[:,1], keys=keys)
	return index


# ---------------------------------------------------------------------------
# QLSC.coverage
# ---------------------------------------------------------------------------

def test_scalar_ipix_document():
	'''
	Test that a single ipix value produces a valid one-polygon document.
	'''
	q = QLSC(depth=4)
	document = q.coverage(ipix=17)

	validate_sky_coverage(document)
	assert_json_native(document)

	assert len(document["polygons"]) == 1
	assert document["polygons"][0]["id"] == 17
	assert "value" not in document["polygons"][0], "no values were given, so no value may be invented"
	assert "value" not in document, "no values were given, so there is nothing to describe"


def test_array_ipix_document_preserves_order():
	'''
	Test that an array of ipix values produces one polygon each, in the order given.
	'''
	q = QLSC(depth=3)
	ipix = [100, 3, 57, 3]
	document = q.coverage(ipix=ipix)

	validate_sky_coverage(document)
	assert [entry["id"] for entry in document["polygons"]] == ipix, \
		"pixels must appear in the order given, and duplicates must not be merged away"


def test_meta_records_producer_scheme_and_depth():
	'''
	Test the meta block the spec's Python example (1.3) writes.
	'''
	import qlsc

	q = QLSC(depth=7)
	document = q.coverage(ipix=[0, 1])

	assert document["meta"]["producer"] == f"qlsc {qlsc.__version__}"
	assert document["meta"]["scheme"] == "qlsc"
	assert document["meta"]["depth"] == 7


def test_no_equinox_is_emitted():
	'''
	Test that no equinox is written: qlsc emits the ICRS frame, and spec 3.4 makes an
	equinox for ICRS an error rather than something the consumer ignores.
	'''
	document = QLSC(depth=2).coverage(ipix=[0])

	assert document["frame"] == "icrs"
	assert "equinox" not in document


def test_values_become_polygon_values():
	'''
	Test that a parallel values array is attached to the polygons as floats.
	'''
	q = QLSC(depth=3)
	document = q.coverage(ipix=[10, 11, 12], values=[1, 2, 3],
						  value_name="count", value_units="sources", value_domain=[0, 10])

	validate_sky_coverage(document)
	assert_json_native(document)

	assert [entry["value"] for entry in document["polygons"]] == [1.0, 2.0, 3.0]
	assert all(type(entry["value"]) is float for entry in document["polygons"])
	assert document["value"] == {"name": "count", "units": "sources", "domain": [0.0, 10.0]}


def test_numpy_inputs_produce_json_native_output():
	'''
	Test that NumPy ipix and value arrays do not leak NumPy scalars into the document.

	This is the failure the spec's one-line json.dumps hand-off (1.3) hits first:
	np.int64 is not JSON serializable at all, and np.float64 serializes silently only
	because it subclasses float.
	'''
	q = QLSC(depth=5)
	ipix = np.array([5, 6, 7], dtype=np.int64)
	values = np.array([1.5, 2.5, 3.5], dtype=np.float64)

	document = q.coverage(ipix=ipix, values=values, value_name="density")

	assert_json_native(document)
	assert type(document["polygons"][0]["id"]) is int
	assert type(document["meta"]["depth"]) is int

	# the actual hand-off: the standard library encoder, no default= hook, no NaN
	text = json.dumps(document, allow_nan=False)
	validate_sky_coverage(json.loads(text))


def test_empty_ipix_array_is_an_empty_coverage():
	'''
	Test that no pixels produce an empty (but valid) document; spec 1.1 calls an empty
	coverage "a legitimate, explicit state".
	'''
	document = QLSC(depth=4).coverage(ipix=[])

	validate_sky_coverage(document)
	assert document["polygons"] == []


def test_every_pixel_of_a_low_depth_scheme():
	'''
	Test the whole sky at depth 1: every pixel of the scheme, validated.

	Depth 0 and 1 pixels are the largest the scheme can produce, so this is where the
	spec's "smaller than a hemisphere" producer contract (1.2) is under the most strain;
	a depth-0 pixel is one cube face, one sixth of the sphere.
	'''
	for depth in (0, 1):
		q = QLSC(depth=depth)
		document = q.coverage(ipix=list(range(q.nbins)))

		validate_sky_coverage(document)
		assert_json_native(document)
		assert len(document["polygons"]) == q.nbins


def test_vertices_are_the_pixel_corners():
	'''
	Test that the vertices of each entry are exactly what ipix2polygon returns.
	'''
	q = QLSC(depth=6)
	ipix = [0, 1234, 9999, q.nbins-1]
	document = q.coverage(ipix=ipix)

	for pixel,entry in zip(ipix, document["polygons"]):
		expected = q.ipix2polygon(pixel)
		assert np.array_equal(np.asarray(entry["vertices"]), expected)


def test_ra_and_dec_ranges():
	'''
	Test that vertices come out as RA in [0,360) and Dec in [-90,90] (spec 1.2).
	'''
	q = QLSC(depth=2)
	document = q.coverage(ipix=list(range(q.nbins)))

	vertices = np.array([v for entry in document["polygons"] for v in entry["vertices"]])
	assert ((0 <= vertices[:,0]) & (vertices[:,0] < 360)).all()
	assert ((-90 <= vertices[:,1]) & (vertices[:,1] <= 90)).all()


def test_value_description_without_values_is_refused():
	'''
	Test that describing values that were not provided raises rather than producing a
	descriptor for nothing.
	'''
	q = QLSC(depth=3)
	with pytest.raises(ValueError, match="not provided"):
		q.coverage(ipix=[1, 2], value_name="count")


def test_mismatched_values_length_is_refused():
	'''
	Test that a values array of the wrong length raises rather than silently truncating.
	'''
	q = QLSC(depth=3)
	with pytest.raises(ValueError, match="does not match"):
		q.coverage(ipix=[1, 2, 3], values=[1, 2])


def test_non_finite_value_is_refused():
	'''
	Test that a NaN value raises: spec 4.4 requires a present value to be a finite
	number, and json.dumps would otherwise emit the non-standard literal NaN.
	'''
	q = QLSC(depth=3)
	with pytest.raises(ValueError, match="finite"):
		q.coverage(ipix=[1, 2], values=[1.0, float("nan")])


def test_multidimensional_ipix_is_refused():
	'''
	Test that a 2-D ipix array raises rather than being flattened.
	'''
	q = QLSC(depth=3)
	with pytest.raises(ValueError, match="one-dimensional"):
		q.coverage(ipix=np.array([[1, 2], [3, 4]]))


def test_reserved_meta_keys_cannot_be_replaced():
	'''
	Test that a caller cannot overwrite the provenance qlsc records about itself.
	'''
	q = QLSC(depth=3)
	with pytest.raises(ValueError, match="cannot be replaced"):
		q.coverage(ipix=[1], meta={"producer": "somebody else"})


def test_extra_meta_is_carried_through():
	'''
	Test that extra meta keys reach the document; spec 1.1 calls meta the documented
	escape hatch for provenance.
	'''
	q = QLSC(depth=3)
	document = q.coverage(ipix=[1], meta={"query": "select everything"})

	validate_sky_coverage(document)
	assert document["meta"]["query"] == "select everything"
	assert document["meta"]["scheme"] == "qlsc"


# ---------------------------------------------------------------------------
# QLSCIndex.coverage
# ---------------------------------------------------------------------------

def test_index_coverage_round_trip():
	'''
	Round trip: generate points, index them, export the coverage, and check the document
	against the spec, the point count, and the pixel geometry.
	'''
	n_points = 5000
	points = sunflower_points_on_sphere(n_points)
	index = _index(points)

	document = index.coverage(depth=4)

	# --- the document is a valid sky-coverage v1 document ---
	validate_sky_coverage(document)
	assert_json_native(document)
	validate_sky_coverage(json.loads(json.dumps(document, allow_nan=False)))

	assert document["meta"]["depth"] == 4
	assert document["value"]["name"] == "count"

	# --- the counts account for every point, exactly once ---
	counts = [entry["value"] for entry in document["polygons"]]
	assert sum(counts) == n_points
	assert all(count >= 1 for count in counts), "a pixel holding no points must not be in the document"

	# --- the ids are the pixels the points actually fall in, and the counts are right ---
	q4 = QLSC(depth=4)
	expected = dict()
	for ra,dec in points:
		pixel = int(q4.ang2ipix(ra, dec))
		expected[pixel] = expected.get(pixel, 0) + 1

	assert {entry["id"]: entry["value"] for entry in document["polygons"]} == \
		   {pixel: float(count) for pixel,count in expected.items()}

	# --- every polygon's vertices are that pixel's corners ---
	for entry in document["polygons"]:
		assert np.allclose(np.asarray(entry["vertices"]), q4.ipix2polygon(entry["id"]))

	index.close()


def test_index_coverage_default_depth():
	'''
	Test the documented default depth, and that it follows a shallower index down.
	'''
	from qlsc.qlsc_base import DEFAULT_COVERAGE_DEPTH

	points = sunflower_points_on_sphere(500)

	deep = _index(points, depth=30)
	assert deep.coverage()["meta"]["depth"] == DEFAULT_COVERAGE_DEPTH
	deep.close()

	shallow = _index(points, depth=2)
	assert shallow.coverage()["meta"]["depth"] == 2
	shallow.close()


def test_index_coverage_at_the_index_depth():
	'''
	Test the shift-by-zero case: counting at the index's own depth.
	'''
	points = sunflower_points_on_sphere(300)
	index = _index(points, depth=6)

	document = index.coverage(depth=6)
	validate_sky_coverage(document)
	assert sum(entry["value"] for entry in document["polygons"]) == 300

	q6 = QLSC(depth=6)
	pixels = set(int(q6.ang2ipix(ra, dec)) for ra,dec in points)
	assert set(entry["id"] for entry in document["polygons"]) == pixels

	index.close()


def test_index_coverage_finer_depth_is_refused():
	'''
	Test that asking for a depth finer than the index's raises rather than quietly
	returning the index's own depth.
	'''
	index = _index(sunflower_points_on_sphere(50), depth=5)
	with pytest.raises(ValueError, match="finer than"):
		index.coverage(depth=6)
	index.close()


def test_index_coverage_of_an_empty_index():
	'''
	Test that an empty index exports an empty coverage rather than failing.
	'''
	index = QLSCIndex(qlsc=QLSC(depth=10))
	document = index.coverage(depth=3)

	validate_sky_coverage(document)
	assert document["polygons"] == []
	index.close()


def test_index_coverage_file_backed(tmp_path):
	'''
	Test that the export works against a file-backed index too, where the connection is
	opened and closed around the query.
	'''
	points = sunflower_points_on_sphere(400)
	filepath = tmp_path / "coverage.qlscidx"
	index = QLSCIndex(qlsc=QLSC(depth=8), filepath=filepath)
	index.add_points(ra=points[:,0], dec=points[:,1])

	document = index.coverage(depth=3)
	validate_sky_coverage(document)
	assert sum(entry["value"] for entry in document["polygons"]) == 400
	index.close()


def test_index_coverage_value_domain():
	'''
	Test that a pinned colour-scale domain reaches the document (spec 1.1).
	'''
	index = _index(sunflower_points_on_sphere(200), depth=8)
	document = index.coverage(depth=2, value_domain=[0, 250])

	validate_sky_coverage(document)
	assert document["value"]["domain"] == [0.0, 250.0]
	index.close()


def test_index_coverage_depth_must_be_an_integer():
	'''
	Test that a fractional depth raises rather than being truncated to another depth.
	'''
	index = _index(sunflower_points_on_sphere(50), depth=10)
	with pytest.raises(TypeError, match="integer"):
		index.coverage(depth=3.5)
	index.close()


def test_index_coverage_at_the_poles_and_seam():
	'''
	Test coverage of points at the poles and either side of ra=0, where the pixel
	polygons meet the dateline and the projection's degenerate points.
	'''
	points = [[0., 90.], [180., 90.], [0., -90.], [359.999, 0.1], [0.001, -0.1], [0., 0.]]
	index = _index(points, depth=30)

	document = index.coverage(depth=6)
	validate_sky_coverage(document)
	assert sum(entry["value"] for entry in document["polygons"]) == len(points)
	index.close()


def test_scalar_ipix_accepts_a_size_one_values_array():
	'''
	Test the scalar-like convention for values: a size-1 array is one value, as a size-1
	array is a scalar value everywhere else in qlsc.
	'''
	q = QLSC(depth=3)

	for values in (3.0, np.float64(3.0), [3.0], np.array([3.0]), np.array(3.0)):
		document = q.coverage(ipix=17, values=values, value_name="count")
		validate_sky_coverage(document)
		assert [entry["value"] for entry in document["polygons"]] == [3.0]

	# and a length-1 ipix array with a single value is the same one-pixel document
	document = q.coverage(ipix=[17], values=3.0, value_name="count")
	assert [entry["value"] for entry in document["polygons"]] == [3.0]


def test_scalar_ipix_rejects_several_values():
	'''
	Test that more than one value for one pixel still raises.
	'''
	q = QLSC(depth=3)
	with pytest.raises(ValueError, match="does not match"):
		q.coverage(ipix=17, values=[3.0, 4.0])


def test_two_dimensional_values_are_refused():
	'''
	Test that a 2-D values array raises rather than being flattened.
	'''
	q = QLSC(depth=3)
	with pytest.raises(ValueError, match="one-dimensional"):
		q.coverage(ipix=[1, 2, 3, 4], values=np.array([[1., 2.], [3., 4.]]))


@pytest.mark.parametrize("value", ["3", b"3", np.str_("3"), True, np.bool_(False)])
def test_string_and_boolean_values_are_refused(value):
	'''
	Test that a value float() would happily convert but that is not a number is refused.

	Reading the count "3" out of a text column and getting a coverage map back is a data
	problem worth hearing about, and a boolean is a flag mistaken for a measurement.
	'''
	q = QLSC(depth=3)
	with pytest.raises(TypeError, match="must be a number"):
		q.coverage(ipix=[1], values=[value])


def test_string_and_boolean_domain_bounds_are_refused():
	'''
	Test that the same strictness applies to the colour-scale domain.
	'''
	q = QLSC(depth=3)
	with pytest.raises(TypeError, match="must be a number"):
		q.coverage(ipix=[1], values=[1.0], value_domain=["0", "10"])
	with pytest.raises(TypeError, match="must be a number"):
		q.coverage(ipix=[1], values=[1.0], value_domain=[False, True])
