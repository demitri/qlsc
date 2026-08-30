
'''
Cross-matching two QLSC indices: the spatial join of two catalogs.

This is the module behind :py:func:`qlsc.crossmatch`, which is the public interface.
'''

from typing import Union

import numpy as np

from .qlsc_base import QLSCIndex, _degrees

try:
	# astropy is an optional dependency; the name is needed only to resolve the
	# "Quantity" annotations below. _degrees does the actual conversion, and it checks
	# qlsc_base's own availability flag.
	from astropy.units import Quantity
except ImportError:
	pass

# The fields of the returned table, in order. Keys are object arrays because a point may
# have been stored without one, in which case its key is None; a fixed-width string dtype
# would turn that into the four characters "None".
_PAIR_FIELDS = [("ra_a", np.double), ("dec_a", np.double), ("ra_b", np.double), ("dec_b", np.double)]
_PAIR_FIELDS_WITH_KEYS = [("ra_a", np.double), ("dec_a", np.double), ("key_a", object),
						  ("ra_b", np.double), ("dec_b", np.double), ("key_b", object)]

def _pair_dtype(return_keys:bool) -> np.dtype:
	'''
	Build the dtype of the table :py:func:`crossmatch` returns.

	:param return_keys: whether the key fields are included
	:returns: the dtype
	'''
	fields = _PAIR_FIELDS_WITH_KEYS if return_keys else _PAIR_FIELDS
	return np.dtype(fields + [("distance", np.double)])

def crossmatch(index_a:QLSCIndex=None, index_b:QLSCIndex=None,
			   radius:Union[float, "Quantity"]=None,
			   nearest_only:bool=False, return_keys:bool=True) -> np.recarray:
	'''
	Cross-match two catalogs: find every pair of points, one from each index, that lie
	within ``radius`` of each other.

	This is the spatial join two catalogs need in order to be compared — "which source in
	my catalog is which source in theirs" — performed locally against two
	:py:class:`QLSCIndex` objects, in memory or on disk.

	.. code-block:: python

		pairs = qlsc.crossmatch(my_catalog, reference_catalog, radius=1/3600, nearest_only=True)
		for pair in pairs:
			print(pair.key_a, pair.key_b, pair.distance)

	The result is a :class:`numpy.recarray` with one row per matched pair and the fields
	``ra_a``, ``dec_a``, ``key_a``, ``ra_b``, ``dec_b``, ``key_b``, ``distance`` — the key
	fields only when ``return_keys`` is set, which it is by default. ``distance`` is the
	angular separation of the pair in degrees. The key fields have an object dtype and hold
	``None`` for a point that was stored without a key, which is not the same as a point
	whose key is the string "None".

	**Ordering.** The pairs are grouped by the point of ``index_a``, in ascending
	``(ra_a, dec_a)`` order, and within a group by increasing separation, ties broken by
	ascending ``(ra_b, dec_b)``. The order is therefore fully determined by the data — it
	does not depend on which index the search iterated over, or on the order the databases
	return their rows.

	**nearest_only.** Set it to keep only the closest counterpart of each ``index_a``
	point. A point of ``index_a`` with no counterpart within ``radius`` is absent from the
	result either way; ``nearest_only`` reduces the rows of a point that has several to
	one, it does not add a row for a point that has none. Note that cross-matching an index
	with itself pairs every point with itself at a separation of zero, so ``nearest_only``
	returns those self-pairs, not each point's nearest neighbour.

	**Method.** The search iterates over whichever index holds fewer points and runs a cone
	search of ``radius`` against the other for each of them, which is the same z-order range
	search :py:func:`QLSCIndex.radial_query` performs: the ipix ranges covering the cone are
	computed at depth 30, mapped down to the searched index's depth and merged, the
	candidate rows are selected by those ranges, and every candidate is then filtered
	exactly with ``sindist``. The result is exactly the set of pairs within the radius, with
	no approximation from the pixelization. The two indices do **not** have to share a
	depth: each search uses the depth of the index being searched, and the coarse ranges are
	a superset whichever depth that is. Only the choice of which side to iterate depends on
	the sizes, and it never changes the answer — the relation is symmetric.

	The cost is therefore one cone search per point of the smaller catalog, and it does not
	depend much on the size of the larger one. Measured on in-memory indices at depth 30,
	both cases ran at about 15,000 points of the smaller catalog per second: 50,000 against
	200,000 points at a one arcsecond radius (22,593 pairs) took 3.3 seconds, and 20,000
	against 200,000 points took 1.3 seconds.

	The points of the smaller index are read into memory (its coordinates and, if asked
	for, its keys); the larger index is only ever queried. Cross-matching two catalogs of
	which the smaller has hundreds of millions of points will therefore want more memory
	than this makes available, which is the documented limit of this implementation.

	.. note:: The radius is one value for the whole match. :py:func:`QLSCIndex.radial_query`
	   accepts one radius per target, but there is nothing here for such an array to be
	   parallel *to*: the targets are the points of an index, and the order in which an
	   index yields its points is not part of qlsc's API, so a caller could not line a
	   per-point radius up with them. Giving that a meaning needs an ordering contract on
	   the index first, which is a larger decision than this function; until then, a
	   variable radius is expressed by calling :py:func:`QLSCIndex.radial_query` with an
	   explicit target list.

	:param index_a: the first catalog; the result is grouped by, and ``nearest_only`` reduced to, the points of this index
	:param index_b: the second catalog
	:param radius: the match radius (degrees) in the range [0,180], accepts a `float` value or `astropy.units.Quantity`
	:param nearest_only: if set to True, keep only the closest point of ``index_b`` for each point of ``index_a``
	:param return_keys: if set to True (the default), include the ``key_a`` and ``key_b`` fields
	:returns: a :class:`numpy.recarray` of matched pairs
	'''
	for index,name in ((index_a, "index_a"), (index_b, "index_b")):
		if index is None:
			raise ValueError(f"A value for '{name}' must be specified.")
		if not isinstance(index, QLSCIndex):
			raise TypeError(f"The '{name}' parameter must be a QLSCIndex object; "
							f"was given a '{type(index).__name__}'.")
		index._ensure_open()

	if radius is None:
		raise ValueError("A value for 'radius' must be specified.")
	radius = _degrees(radius, "radius")
	if not (0 <= radius <= 180):
		raise ValueError(f"The value for radius must be in the range [0,180]; was given '{radius}'.")

	dtype = _pair_dtype(return_keys)

	# Iterate over the smaller catalog and query the larger: the work is one cone search
	# per iterated point, so the fewer of those the better, and a cone search of a large
	# index costs no more than one of a small index (it is a z-order range scan, not a
	# table scan). Counting both is two COUNT(*) scans, which is cheap beside the join.
	count_a = index_a.number_of_points
	count_b = index_b.number_of_points
	iterate_a = count_a <= count_b

	if count_a == 0 or count_b == 0:
		# nothing can match; return the empty table rather than running a query per point
		# of an empty catalog
		return np.recarray(0, dtype=dtype)

	source = index_a if iterate_a else index_b
	searched = index_b if iterate_a else index_a

	source_ra, source_dec, source_keys = source._all_points(with_keys=return_keys)

	# return_counts is set not for the counts themselves but because it selects the batch
	# contract outright: without it, a source index holding exactly one point would be
	# read as a scalar query and the result would have no "target_index" field.
	matches, counts = searched.radial_query(ra=source_ra, dec=source_dec, radius=radius,
											return_key=return_keys, return_distance=True,
											return_counts=True)

	if len(matches) == 0:
		return np.recarray(0, dtype=dtype)

	# "target_index" points back into the iterated catalog, once per match
	source_index = matches["target_index"]

	if iterate_a:
		ra_a, dec_a = source_ra[source_index], source_dec[source_index]
		ra_b, dec_b = matches["ra"], matches["dec"]
	else:
		ra_a, dec_a = matches["ra"], matches["dec"]
		ra_b, dec_b = source_ra[source_index], source_dec[source_index]

	if return_keys:
		matched_keys = list(matches["key"])
		source_matched_keys = [source_keys[n] for n in source_index]
		key_a, key_b = (source_matched_keys, matched_keys) if iterate_a else (matched_keys, source_matched_keys)

	distance = matches["distance"]

	# Group the pairs by the index_a point. The points of an index are unique by
	# coordinate (the table has a UNIQUE (ra,dec) constraint), so the coordinate pair
	# identifies the point whether or not this side was the one iterated over - which is
	# what makes the ordering independent of the direction of the search.
	_, group = np.unique(np.column_stack((ra_a, dec_a)), axis=0, return_inverse=True)
	group = np.asarray(group).ravel()

	# last key is primary: group, then separation, then the index_b coordinate as a
	# deterministic tie-break between counterparts at exactly the same separation
	order = np.lexsort((dec_b, ra_b, distance, group))

	if nearest_only:
		# after the sort each group is contiguous and its closest pair is first
		sorted_group = group[order]
		first_of_group = np.ones(sorted_group.shape, dtype=bool)
		first_of_group[1:] = sorted_group[1:] != sorted_group[:-1]
		order = order[first_of_group]

	result = np.recarray(len(order), dtype=dtype)
	result["ra_a"] = ra_a[order]
	result["dec_a"] = dec_a[order]
	result["ra_b"] = ra_b[order]
	result["dec_b"] = dec_b[order]
	result["distance"] = distance[order]
	if return_keys:
		result["key_a"] = [key_a[n] for n in order]
		result["key_b"] = [key_b[n] for n in order]

	return result
