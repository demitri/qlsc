
'''
Assembly of "sky-coverage" version 1 documents.

The format is a producer-neutral JSON interchange for sets of polygons on the sphere,
each optionally carrying a scalar value; it was specified by the `astplot-js
<https://github.com/demitri/astplot-js>`_ project, whose ``CoverageLayer`` renders such
a document as a filled-polygon coverage map or choropleth in any projection.

This module holds the format's fixed vocabulary and the document assembly shared by
:py:func:`qlsc.QLSC.coverage` and :py:func:`qlsc.QLSCIndex.coverage`. Nothing here is
part of the public API: the public entry points are those two methods, which return a
plain Python ``dict`` that :py:func:`json.dumps` serializes with the standard library
encoder alone (no NumPy scalars survive into the document).

Two properties of the format matter to callers and are guaranteed here:

* **Unknown keys are an error for a consumer**, both at the top level and inside a
  polygon entry, so this module emits only the keys version 1 defines. Producer-specific
  information goes in ``meta``, which the renderer never interprets.
* **A ring must be simple and smaller than a hemisphere**; the filled region is the
  smaller of the two a closed curve bounds. This is a producer contract the consumer
  does not check. Every QLSC pixel satisfies it: the largest possible pixel is one whole
  cube face at depth 0, which is one sixth of the sphere.

qlsc emits coordinates in degrees in the ICRS frame with geodesic edges. ICRS has no
equinox, and the format makes declaring one for ICRS an error, so no ``equinox`` key is
ever written.
'''

import math
from typing import Iterable, Union

import numpy as np

from .version import __version__

# ---------------------------------------------------------------------------
# The fixed vocabulary of sky-coverage v1. These are the only values qlsc emits;
# the format defines other frames and units, but qlsc works in ICRS degrees and
# its pixel edges are declared as great-circle arcs between the corners.
# ---------------------------------------------------------------------------
SKY_COVERAGE_TYPE = "sky-coverage"
SKY_COVERAGE_VERSION = 1
SKY_COVERAGE_UNITS = "deg"
SKY_COVERAGE_FRAME = "icrs"
SKY_COVERAGE_EDGES = "geodesic"

# the name this producer records in "meta"
SKY_COVERAGE_SCHEME = "qlsc"

# keys of "meta" that are always written by qlsc and cannot be overridden by a caller
_RESERVED_META_KEYS = ("producer", "scheme", "depth")

def _coverage_value(value, description:str) -> float:
	'''
	Coerce one per-polygon scalar to a plain Python ``float``.

	The format requires the value to be a finite number; a NaN or infinity would
	serialize as ``NaN``/``Infinity``, which is not valid JSON and which the consumer
	rejects, so it is refused here where the offending entry can still be named.

	:param value: the value to convert
	:param description: how to refer to the value in an error message
	:returns: the value as a float
	'''
	try:
		result = float(value)
	except (TypeError, ValueError) as e:
		raise TypeError(f"The {description} must be a number; was given '{value!r}'.") from e
	if not math.isfinite(result):
		raise ValueError(f"The {description} must be finite; was given '{value!r}'.")
	return result

def _value_descriptor(name:str=None, units:str=None, domain:Iterable=None) -> dict:
	'''
	Build the document's optional "value" descriptor.

	Every field of the descriptor is itself optional; a descriptor with no fields
	carries no information and is omitted from the document entirely rather than
	written as an empty object.

	:param name: what the per-polygon values are (e.g. "count"), used as legend text
	:param units: the units of the per-polygon values, used as legend text
	:param domain: the [min,max] value range to pin the consumer's colour scale to
	:returns: the descriptor, or None if no field was given
	'''
	descriptor = dict()

	if name is not None:
		if not isinstance(name, str):
			raise TypeError(f"The 'value_name' must be a string; was given '{name!r}'.")
		descriptor["name"] = name

	if units is not None:
		if not isinstance(units, str):
			raise TypeError(f"The 'value_units' must be a string; was given '{units!r}'.")
		descriptor["units"] = units

	if domain is not None:
		bounds = np.asarray(domain)
		if bounds.shape != (2,):
			raise ValueError(f"The 'value_domain' must be a [min,max] pair; was given shape {bounds.shape}.")
		low = _coverage_value(bounds[0], "'value_domain' minimum")
		high = _coverage_value(bounds[1], "'value_domain' maximum")
		if low > high:
			raise ValueError(f"The 'value_domain' minimum ({low}) must not be greater than its maximum ({high}).")
		descriptor["domain"] = [low, high]

	return descriptor if descriptor else None

def _coverage_meta(depth:int, meta:dict=None) -> dict:
	'''
	Build the document's "meta" object.

	``producer``, ``scheme`` and ``depth`` are always written. A caller-supplied
	``meta`` adds to them, but may not redefine them: a document claiming to have been
	produced by something other than the qlsc version that actually produced it is
	worse than no provenance at all, so the collision raises rather than one value
	quietly winning.

	:param depth: the QLSC depth the polygons were generated at
	:param meta: additional producer metadata to record
	:returns: the "meta" object
	'''
	document_meta = {
		"producer" : f"qlsc {__version__}",
		"scheme"   : SKY_COVERAGE_SCHEME,
		"depth"    : int(depth),
	}

	if meta is not None:
		if not isinstance(meta, dict):
			raise TypeError(f"The 'meta' parameter must be a dictionary; was given '{type(meta).__name__}'.")
		collisions = sorted(set(meta) & set(_RESERVED_META_KEYS))
		if collisions:
			raise ValueError(f"These 'meta' keys are written by qlsc itself and cannot be replaced: {collisions}.")
		document_meta.update(meta)

	return document_meta

def _polygon_entry(ipix:int, vertices:np.ndarray, value:Union[float,None]) -> dict:
	'''
	Build one entry of the document's "polygons" array.

	The ipix number becomes the entry's ``id``, which the format treats as an opaque
	producer-defined identity and surfaces verbatim in hit-testing and readouts.

	:param ipix: the pixel identifier number
	:param vertices: the pixel's corners as an (n,2) array of [ra,dec] pairs in degrees
	:param value: the scalar attached to this pixel, or None for an entry with no value
	:returns: the polygon entry
	'''
	if vertices.ndim != 2 or vertices.shape[1] != 2:
		raise ValueError(f"A polygon's vertices must be an (n,2) array of [ra,dec] pairs; was given shape {vertices.shape}.")
	if vertices.shape[0] < 3:
		raise ValueError(f"A polygon must have at least 3 vertices; was given {vertices.shape[0]}.")

	# .tolist() is what keeps the document serializable by the standard library JSON
	# encoder: it produces Python floats, not the np.float64 values indexing the array
	# would give, and json.JSONEncoder refuses those.
	entry = {"id" : int(ipix), "vertices" : vertices.tolist()}

	if value is not None:
		entry["value"] = value

	return entry

def _sky_coverage_document(polygons:list, depth:int, value_name:str=None, value_units:str=None,
						   value_domain:Iterable=None, meta:dict=None) -> dict:
	'''
	Assemble a complete sky-coverage v1 document around a list of polygon entries.

	:param polygons: the entries produced by :py:func:`_polygon_entry`
	:param depth: the QLSC depth the polygons were generated at
	:param value_name: what the per-polygon values are, for the consumer's legend
	:param value_units: the units of the per-polygon values, for the consumer's legend
	:param value_domain: the [min,max] range to pin the consumer's colour scale to
	:param meta: additional producer metadata to record alongside qlsc's own
	:returns: the document as a plain Python dictionary
	'''
	document = {
		"type"    : SKY_COVERAGE_TYPE,
		"version" : SKY_COVERAGE_VERSION,
		"units"   : SKY_COVERAGE_UNITS,
		"frame"   : SKY_COVERAGE_FRAME,
		"edges"   : SKY_COVERAGE_EDGES,
		# no "equinox" key: ICRS has no equinox and the format makes declaring one an error
		"meta"    : _coverage_meta(depth, meta),
		"polygons": polygons,
	}

	descriptor = _value_descriptor(value_name, value_units, value_domain)
	if descriptor is not None:
		document["value"] = descriptor

	return document
