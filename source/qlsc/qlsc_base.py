
import os
import sys
import math
import pathlib
import sqlite3
import logging
import contextlib
from typing import Iterable, Tuple, Union
from urllib.parse import urlparse
from math import pow, sin, asin, sqrt, degrees
from math import radians as deg2rad
import collections

import numpy as np

try:
	import astropy
	astropy_available = True
	from astropy.units import Quantity
	import astropy.units as u
except ImportError:
	astropy_available = False

try:
	import pandas as pd
	pandas_available = True
except ImportError:
	pandas_available = False

# functions defined in "qlsc_c_module.c" as named in the array "qlsc_methods" are in the qlsc.q3c namespace
from . import q3c
from .q3c import sindist

from .utilities import _normalize_ang, _as_float_scalar

# really important for the ipix values!
sqlite3.register_adapter(np.int64, int)

logger = logging.getLogger("qlsc")

def _canonical_position(ra:float, dec:float) -> Tuple[float, float]:
	'''
	Reduce equivalent spherical coordinate representations of the same point to a single form:
	ra is periodic in 360 degrees, and at the poles the ra value is irrelevant.
	'''
	if abs(dec) == 90:
		return (0.0, dec)
	return (ra % 360.0, dec)

def _sindist2deg(value:float) -> float:
	'''
	Convert the q3c "sindist" of two points to their angular separation in degrees.

	q3c's sindist is the haversine quantity sin^2(d/2), so d = 2*asin(sqrt(sindist)).
	Rounding can push the value a hair outside [0,1] for coincident or antipodal
	pairs, which asin would reject; the value is clamped to that range first. Any
	larger excursion is a real error and is left to raise.
	'''
	if not (-1e-9 <= value <= 1 + 1e-9):
		raise ValueError(f"A sindist value outside [0,1] cannot be an angular separation; was given '{value}'.")
	return degrees(2. * asin(sqrt(min(max(value, 0.), 1.))))

def _polygon_ra_dec(polygon:Iterable) -> Tuple[np.ndarray, np.ndarray]:
	'''
	Validate a polygon given as an array of [ra,dec] vertex pairs and return its ra
	and dec values as separate contiguous double arrays, as the C layer requires.

	The vertex count is not checked here: the C layer knows the limit
	(``Q3C_MAX_N_POLY_VERTEX``) and reports it.

	:param polygon: an array of the polygon's vertices as [ra,dec] pairs, shape (n,2)
	:returns: a tuple of (ra, dec) arrays
	'''
	poly = np.asarray(polygon, dtype=np.double)
	if poly.ndim != 2 or poly.shape[1] != 2:
		raise ValueError(f"The polygon must be an array of [ra,dec] pairs with shape (n,2); was given shape {poly.shape}.")
	if not np.isfinite(poly).all():
		raise ValueError("The polygon vertex values must all be finite.")
	if not ((-90 <= poly[:,1]).all() and (poly[:,1] <= 90).all()):
		raise ValueError("The polygon vertex dec values must be in the range [-90,90].")

	return np.ascontiguousarray(poly[:,0]), np.ascontiguousarray(poly[:,1])

def _degrees(value:Union[float, "Quantity"], name:str) -> float:
	'''
	Interpret an angle given either as a number of degrees or as an astropy
	:class:`~astropy.units.Quantity`.

	:param value: the angle
	:param name: the parameter name, for error messages
	:returns: the angle in degrees as a float
	'''
	if astropy_available and isinstance(value, Quantity):
		value = value.to(u.deg).value
	return _as_float_scalar(value, name)

def _match_result(match_ra:list, match_dec:list, match_key:list, match_distance:list,
				  return_key:bool, return_distance:bool):
	'''
	Assemble the return value shared by the region query methods.

	:param match_ra: the ra values of the matches
	:param match_dec: the dec values of the matches
	:param match_key: the key values of the matches, if they were collected
	:param match_distance: the angular separations of the matches, if they were collected
	:param return_key: whether the "key" field was requested
	:param return_distance: whether the "distance" field was requested
	:returns: an array of shape (n,2) if neither flag is set, otherwise a :class:`numpy.recarray`
	'''
	if return_key or return_distance:
		arrays = [match_ra, match_dec]
		names = ["ra", "dec"]
		if return_key:
			arrays.append(match_key)
			names.append("key")
		if return_distance:
			arrays.append(match_distance)
			names.append("distance")
		# np.rec is the public spelling of this module; np.core.records is
		# deprecated in NumPy 1.25+ and gone from the NumPy 2 namespace.
		return np.rec.fromarrays(arrays, names=",".join(names))
	else:
		return np.squeeze(np.dstack((match_ra,match_dec)))

class QLSC:
	'''
	A class that describes a quadrilaterized spherical cube (QLSC).

	The quadrilaterized sphere is divided into six faces. The bin level defines the number times a
	face is subdivided. A bin level of 1 divides each face into four "pixels"; bin level two divides
	each of those pixels again into four, etc.

	At the lowest resolution where bin level=0 there is 1 bin per face for a total of 6 bins.
	The next bin level 1 divides each face into four bins, i.e. a total of 24 bins.
	Higher bin levels further divide each bin four more times.

	A scheme can be defined by the bin level where:

		.. code-block::

			no. bins per face        = (2**depth)**2
			number of bins on sphere = no. bins per face * 6

	Each bin has a pixel number called an "ipix", which is an integer encoded with the location of the bin.

	The default value of `depth` = 30 that divides the sphere into 6,442,450,944 bins (i.e. nside=1,073,741,824),
	is the one used by the `Q3C PostgreSQL plugin <https://github.com/segasai/q3c>`_.
	This value is used if neither `nside` or `bin_size` are specified.

	:param depth: number of times each bin in each face is subdivided by four
	'''
	def __init__(self, depth:int=30):

		# k = depth in Chan paper

		# The Q3C code only supports up to depth=30
		if not (0 <= depth <= 30):
			raise ValueError(f"The value for 'depth' must be an integer on [0,30]; was given '{depth}'.")

		# "nside" as defined by Q3C is the number of bins along one side of the cube face.
		# This means the number of bins on a cube face is nside*nside, and the total
		# number of bins is 6*nside*nside.
		#
		nside = int(pow(2, depth))

		self.depth = depth

		# store hprm struct values - hprm is a pointer (struct) to main Q3C structure used by most Q3C calls
		self._hprm = q3c.init_q3c(nside=nside)
		if self._hprm is None:
			raise Exception("Internal Q3C data structure could not be created.")

	def __repr__(self):
		return f"<{self.__class__.__module__.split('.')[0]}.{self.__class__.__name__} at {hex(id(self))}, depth={self.depth}>"

	def ang2ipix(self, ra:Union[float,Iterable]=None, dec:Union[float,Iterable]=None, points:Iterable=None) -> Union[float,np.ndarray]:
		'''
		Convert ra,dec positions to an ipix number.

		This method accepts scalar values or NumPy arrays. Slices of any other NumPy array are accepted.
		Note that arrays are required to be ``dtype.double``; copies will be made of arrays of any other type.

		Right ascension values outside of the range [0,360] are handled efficiently.
		Declination values outside of [-90,90] are handled, but may result in arrays being copied. (This is only important when working with extremely large arrays.)

		Example usage:

		.. code-block:: python

			ipix = ang2ipix(45,60)
			ipix = ang2ipix(ra=45,dec=60)
			ipix = ang2ipix(points=np.array([[45,60],[0,55]])
			ipix = ang2ipix(ra=catalog[:,4], dec=catalog[:,5])

		:param ra: right ascension (degrees)
		:param dec: declination (degrees)
		:param points: an array of [ra,dec] points, shape=(n,2)
		:returns: the ipix integer if a single point is provided, an array of ipix integers if an array of values is provided
		'''

		# .. todo:: support pd.Series for ra, dec.

		# check "ra","dec" are both set
		if points is None:
			if True in [x is None for x in [ra,dec]]: # can't use any/all with possible "0" values
				raise ValueError(f"Values for both 'ra' and 'dec' must be specified.")
		else:
			if True in [x is not None for x in [ra,dec]]:
				raise ValueError("The parameter 'points' was given; 'ra' or 'dec' should not also be specified.")

		# q3c.ang2ipix unwraps ra values to [0,360], but it does not do the same for dec.
		# It will silently truncate to [-90,90].

		# handle simple scalar values
		#
		if np.isscalar(ra) and np.isscalar(dec):
			if not (math.isfinite(ra) and math.isfinite(dec)):
				# NaN passes the range check below unnoticed and reaches undefined
				# float->int conversions in the C code (upstream q3c issue #52)
				raise ValueError(f"The 'ra' and 'dec' values must be finite; was given ra={ra}, dec={dec}.")
			if not (-90 <= dec <= 90):
				ra, dec = _normalize_ang(ra,dec)
			return q3c.ang2ipix(self._hprm, ra, dec)

		if points is not None:
			# break it down:
			ra = points[:,0]
			dec = points[:,1]

		# the C extension requires arrays of dtype double; copy/convert any other type
		ra = np.asarray(ra, dtype=np.double)
		dec = np.asarray(dec, dtype=np.double)

		if not (np.isfinite(ra).all() and np.isfinite(dec).all()):
			raise ValueError("The 'ra' and 'dec' arrays must contain only finite values.")

		if not ((-90 <= dec).all() and (dec <= 90).all()):
			ra, dec = _normalize_ang(ra=ra, dec=dec, copy=True)

		if ra.ndim > 0 and ra.shape[0] == 1: # the ndim check leaves 0-d arrays to the C scalar path
			# A length-1 array would match the C extension's scalar signature (via a
			# __float__ conversion NumPy has deprecated), returning a bare int for
			# array input; call the scalar path explicitly and wrap the result.
			return np.array([q3c.ang2ipix(self._hprm, float(ra[0]), float(dec[0]))], dtype=np.int64)

		return q3c.ang2ipix(self._hprm, ra, dec)

	def ang2ipix_xy(self, ra:float=None, dec:float=None, points=None) -> dict:
		'''
		Convert an ra,dec position to an ipix number; also returns the face number, and (x,y) location on the square face.

		:param ra: right ascension angle (degrees)
		:param dec: declination angle (degrees), must be in [-90,90]
		:param points: an array of ra,dec coordinates (shape: (n,2)) NOT YET SUPPORTED
		:returns: a dictionary with the keys: [`facenum`, `ipix`, `x`, `y`]
		'''

		ra = _as_float_scalar(ra, "ra")
		dec = _as_float_scalar(dec, "dec")

		# The underlying C code forces dec > 90 to 90 and dec < -90 to 90.
		# Either normalize the angle here or raise an exception.
		# The underlying C code *does* unwrap RA.
		if not (-90 <= dec <= 90):
			raise ValueError(f"Values of 'dec' must be in the range [-90,90] (was given '{dec}').")

		return q3c.ang2ipix_xy(self._hprm, ra, dec)

	def _validated_ipix(self, ipix:Union[int,Iterable]) -> Union[int, np.ndarray]:
		'''
		Validate scalar or array-like ipix input for this scheme.

		Returns a Python ``int`` for scalar input (a Python or NumPy integer, or a
		0-d array) and a C-contiguous ``np.int64`` array for array-like input,
		mirroring the scalar/array convention of :py:func:`ang2ipix`.

		Non-integer values (including floats that happen to be whole numbers) and
		values outside [0, nbins) raise ValueError rather than being truncated or
		clamped: both would silently return the coordinates of a different pixel.

		:param ipix: an ipix number or an array-like of ipix numbers
		:returns: an ``int`` for scalar input, an ``np.int64`` array otherwise
		'''
		values = np.asarray(ipix)

		if values.size == 0:
			# nothing to type check or bounds check; preserve the input shape
			return np.empty(values.shape, dtype=np.int64)

		if not np.issubdtype(values.dtype, np.integer):
			raise ValueError(f"The ipix value(s) must be integers; was given a value of type '{values.dtype}'.")

		# perform a bounds check
		out_of_range = values[(values < 0) | (values >= self.nbins)]
		if out_of_range.size > 0:
			raise ValueError(f"An ipix number is out of bounds for this scheme (depth={self.depth}); "
							 f"range: [0,{self.nbins-1}], first offending value: {out_of_range.flat[0]}.")

		if values.ndim == 0:
			return int(values)

		return np.ascontiguousarray(values, dtype=np.int64)

	def ipix2ang(self, ipix:Union[int,Iterable]) -> Tuple[Union[float,np.ndarray], Union[float,np.ndarray]]:
		'''
		Convert an ipix number to the ra,dec coordinate of the center of the pixel, in degrees.

		This method accepts a single ipix value or an array-like of them (a NumPy array,
		list, or tuple of integers). A scalar returns a tuple of two floats; array input
		returns a tuple of two arrays of the same shape, following the same scalar/array
		convention as :py:func:`ang2ipix`.

		This method matches the behavior of the `Q3C PostgreSQL plugin <https://github.com/segasai/q3c>`_
		(v2.0.1 and later). Note that prior to qlsc v1.1, this method returned the "lower left"
		corner of the pixel, matching older versions of the plugin. To work with the pixel
		corners, see :py:func:`ipix2xy` and :py:func:`ipix2polygon`.

		Example usage:

		.. code-block:: python

			ra, dec = qlsc.ipix2ang(294541929365396512)
			ra, dec = qlsc.ipix2ang(catalog_ipix_array)

		:param ipix: an ipix number, or an array-like of ipix numbers
		:returns: the ra,dec of the pixel center in degrees, as a tuple of floats (scalar input) or a tuple of arrays (array input)
		'''
		return q3c.ipix2ang(self._hprm, self._validated_ipix(ipix))

	def ipix2ang_center(self, ipix:Union[int,Iterable]) -> Tuple[Union[float,np.ndarray], Union[float,np.ndarray]]:
		'''
		Convert an ipix number to the ra,dec coordinate of the center of the pixel, in degrees.

		This is a synonym of :py:func:`ipix2ang`, kept for callers who want the "center"
		semantics to be explicit at the call site. Like :py:func:`ipix2ang`, it accepts a
		single ipix value or an array-like of them.

		The pixel "center" here is defined in the x,y coordinates;
		the point projected on the sphere is the one returned.

		:param ipix: an ipix number, or an array-like of ipix numbers
		:returns: the ra,dec of the pixel center in degrees, as a tuple of floats (scalar input) or a tuple of arrays (array input)
		'''
		return self.ipix2ang(ipix)


	def ipix2xy(self, ipix:Union[int,Iterable]) -> Tuple[Union[int,np.ndarray], Union[float,np.ndarray], Union[float,np.ndarray]]:
		'''
		Convert an ipix value to the (x,y) position on the corresponding square face; returns (facenum,x,y).

		The (x,y) returned is the "lower left" corner of the pixel on the face,
		where the face coordinates run from -1 to +1 in both axes.

		This method accepts a single ipix value or an array-like of them (a NumPy array,
		list, or tuple of integers). A scalar returns ``(facenum, x, y)`` as an int and two
		floats; array input returns three arrays of the same shape (the face numbers as
		``np.int64``), following the same scalar/array convention as :py:func:`ang2ipix`.

		:param ipix: an ipix number, or an array-like of ipix numbers
		:returns: (facenum,x,y) as a tuple of scalars (scalar input) or a tuple of arrays (array input)
		'''
		return q3c.ipix2xy(self._hprm, self._validated_ipix(ipix))

	def xy2ang(self, facenum:int=None, x:float=None, y:float=None, points:collections.abc.Iterable=None): # -> Tuple[float,float]:
		'''
		Convert an x,y coordinate pair on the given face number to (ra,dec).

		:param facenum: the number of the cube face (0=top, 1-4=sides, 5=bottom)
		:param x: x coordinate on cube face
		:param y: y coordinate on cube face
		:param points: an array containing (x,y) coordinate pairs; shape = (n,2)
		:returns: tuple of (ra,dec) corresponding to the location on the cube face
		'''
		# .. todo:: validate ra,dec

		if facenum is None:
			raise ValueError("The cube face number ('facenum') must be provided.")
		else:
			if not isinstance(facenum, int):
				# notmally would be ok to cast the value, but this is a safety net to catch 'x' values passed to facenum
				raise ValueError("The value of facenum is expected to be an integer.")
			if not (0 <= facenum <= 5):
				raise ValueError(f"Values of 'facenum' must be an integer between 0 and 5; was given '{facenum}'.")
		if all(v is not None for v in [x, y, points]): # zero is a possible value
			raise ValueError(f"Specify only 'ra' and 'dec', or 'points'")
		if points is None and (x is None or y is None):
			raise ValueError(f"Values for 'x' and 'y' must be specified together; missing one (x={x}, y={y}).")


		if points is None:
			x = _as_float_scalar(x, "x")
			y = _as_float_scalar(y, "y")
			# "or", not "and": either coordinate out of range is an error
			if not (-1 <= x <= 1) or not (-1 <= y <= 1):
				raise ValueError(f"Values of 'x' and 'y' must be in the range [-1,1]; was given ({x},{y}).")

		if points is not None:
			points = np.asarray(points, dtype=np.double)
			if not ((-1 <= points).all() and (points <= 1).all()):
				raise ValueError("Values of 'x' and 'y' must be in the range [-1,1].")
			ang = np.zeros([len(points),2], dtype=np.double)
			idx = 0
			for p in points:
				ang[idx,:] = q3c.xy2ang(p[0], p[1], facenum)
				idx +=1 # seriously, fuck Numpy
			return ang
		else:
			return q3c.xy2ang(x, y, facenum)

	def xy2ipix(self, facenum:int=None, x:float=None, y:float=None) -> int: #, points:collections.abc.Iterable):
		'''
		Convert an x,y coordinate pair on the given face number to the ipix value.

		:param facenum: the number of the cube face (0=top, 1-4=sides, 5=bottom)
		:param x: x coordinate on cube face
		:param y: y coordinate on cube face
		'''
		ra,dec = self.xy2ang(x=x, y=y, facenum=facenum)
		return self.ang2ipix(ra, dec)

	def ipix_down(self, ipix:int) -> Iterable[int]:
		'''
		Returns the four ipix values at the next higher depth.

		:param ipix: the ipix value at the resolution of this object's scheme.
		'''
		if self.depth == 30:
			raise ValueError("QLSC doesn't support depth values greater than 30.")
		elif not (0 <= ipix < self.nbins):
			raise ValueError(f"ipix value out of range for this depth; should be in [0,{self.nbins}).")
		else:
			new_ipix = ipix <<  2 # same as -> int(ipix * 2**(self.nside))
			return list(range(new_ipix, new_ipix+4))

	def ipix_up(self, ipix:int):
		'''
		Returns the ipix value value at the next lower depth.

		:param ipix: the ipix value at the resolution of this object's scheme.
		'''
		if self.depth == 0:
			raise ValueError("Already at lowest depth (0).")
		elif not (0 <= ipix < self.nbins):
			raise ValueError(f"ipix value out of range for this depth; should be in [0,{self.nbins}).")
		else:
			return ipix >> 2

	def ipix_neighbors(self, ipix:int) -> np.ndarray:
		'''
		Return the ipix values of the pixels that touch the given pixel.

		Two pixels "touch" when they share an edge or a single corner, so a pixel in the
		interior of a cube face has the expected eight neighbors. Three special cases
		follow from the geometry of the cube and are not errors:

		* The 24 pixels sitting on a **corner of the cube** (four per face, at every depth
		  above 0) have only **seven** neighbors. Three faces meet at a cube vertex, so
		  three pixels meet there rather than four, and the diagonal step across the vertex
		  reaches no eighth pixel.
		* At ``depth=0`` each face is a single pixel and every face has **four** neighbors —
		  the four faces adjacent to it. (The opposite face is correctly not included.)
		* The pixel itself is never included in the result.

		The returned values are unique and sorted in increasing ipix order. Sorted, because
		the eight directions have no natural ordering once a face boundary is crossed (the
		coordinate axes of two faces meeting at an edge are rotated with respect to each
		other), so any "north, east, ..." ordering would be misleading; unique, because more
		than one direction can reach the same pixel at very low depths. If the direction
		matters, use :py:func:`ipix2xy` on the results and compare positions.

		The calculation is done in the face-plane (x,y) coordinate system. Each of the eight
		directions is probed with a point placed just inside the neighboring pixel, and the
		pixel containing that point is looked up with :py:func:`ang2ipix`. A probe that
		leaves the cube face is placed only a fraction of a pixel past the edge: the
		gnomonic projection of a face remains valid beyond the face itself, but the pixel
		grids of two faces line up exactly only *on* their shared edge and drift apart in
		proportion to the distance past it, so a probe placed a half pixel out could land in
		the wrong row near a cube corner.

		:param ipix: the pixel identifier number; a single value (this method is not vectorized)
		:returns: a sorted array of the unique ipix values of the neighboring pixels
		'''
		value = self._validated_ipix(ipix)
		if not isinstance(value, int):
			raise ValueError(f"This method takes a single ipix value; was given an array of shape {np.shape(ipix)}.")

		nside = self.nside
		bin_length = 2. / nside # the width of a pixel in face (x,y) coordinates
		half = bin_length / 2.

		# How far past a cube face edge a probe point is placed. Crossing an edge by
		# this much moves the transverse coordinate by at most the same amount (the two
		# faces' grids agree exactly on the edge and diverge in proportion to the
		# distance past it), so a fraction of a pixel keeps the probe comfortably inside
		# the intended row while staying many orders of magnitude above the floating
		# point noise of the projection (~1e-16 against ~2e-10 of a pixel at depth 30).
		edge_step = bin_length / 8.

		facenum, x, y = self.ipix2xy(value) # the "lower left" corner of the pixel
		x_center = x + half
		y_center = y + half

		# the pixel's integer position in the face's grid
		i = int(round((x + 1.) / bin_length))
		j = int(round((y + 1.) / bin_length))
		if not (0 <= i < nside and 0 <= j < nside):
			raise RuntimeError(f"The pixel grid position computed for ipix={value} ({i},{j}) is outside the face; this is a bug in qlsc.")

		neighbors = set()

		for step_x in (-1, 0, 1):
			for step_y in (-1, 0, 1):
				if step_x == 0 and step_y == 0:
					continue

				# does this step leave the cube face?
				crosses_x = (step_x < 0 and i == 0) or (step_x > 0 and i == nside - 1)
				crosses_y = (step_y < 0 and j == 0) or (step_y > 0 and j == nside - 1)

				if step_x != 0 and step_y != 0 and crosses_x and crosses_y:
					# A diagonal step across a corner of the cube. Only three pixels meet
					# at a cube vertex - this one and the two reached by the neighboring
					# straight steps - so there is no eighth pixel to find, and probing
					# past the vertex would only rediscover one of the other two.
					continue

				# A step that stays on the face is probed at the center of the neighboring
				# band, the point furthest from any boundary; a step that leaves the face
				# is probed just past the edge, for the reason given above.
				probe_x = x_center + step_x * (half + (edge_step if crosses_x else half))
				probe_y = y_center + step_y * (half + (edge_step if crosses_y else half))

				# q3c.xy2ang is called directly rather than through QLSC.xy2ang, which
				# (correctly, for its own purposes) restricts x and y to [-1,1]: the
				# projection is what carries a probe point past the edge of its face.
				ra, dec = q3c.xy2ang(probe_x, probe_y, facenum)
				neighbors.add(int(q3c.ang2ipix(self._hprm, ra, dec)))

		neighbors.discard(value) # a pixel is not its own neighbor

		return np.array(sorted(neighbors), dtype=np.int64)

	def ipix2face(self, ipix:int) -> int:
		'''
		Return the face number the given ipix value falls on.
		'''
		if not (0 <= ipix < self.nbins):
			raise ValueError(f"ipix value out of range for this depth; should be in [0,{self.nbins}).")
		return ipix // math.pow(self.nside, 2)

	@property
	def nside(self) -> int:
		'''
		Number of bins along one edge of the cube face.
		'''
		return q3c.nside(self._hprm)

	@property
	def nbins(self) -> int:
		'''
		Returns the total number of bins in this pixellation scheme (number of divisions per face * 6).
		'''
		return int(pow(self.nside,2)) * 6

	def face_number(self, ra:float, dec:float, ipix:int=None) -> int:
		'''
		Return the cube face number for the provided coordinate.

		Cube faces are numbered as the following:

		|  face 0 = top
		|  face 1 = -45° ≤ α < 45
		|  face 2 = 45° ≤ α < 135
		|  face 3 = 135° ≤ α < 225
		|  face 4 = 225° ≤ α < 315
		|  face 5 = bottom

		:param ra: right ascension in degrees; must be in [0,360]
		:param dec: declination in degrees; must be in [-90,90]
		:param ipix: the pixel identifier number
		:returns: the cube face number
		'''
		# be careful with all() and any() as 0 is a legitimate value

		if ipix and True in [x is not None for x in [ra,dec]]:
			raise ValueError("Only specify 'ra','dec' OR 'ipix', not all three.")
		if (ra is not None or dec is not None) and (not all([x is not None for x in [ra,dec]])):
			raise ValueError(f"Both parameters 'ra' and 'dec' are required when one is provided (ra={ra}, dec={dec}).")
		if ipix:
			if not (0 <= ipix < self.nbins):
				raise ValueError(f"The ipix number provided ('{ipix}') is outside the range [0,{self.nbins-1}].")
		else:
			# ra,dec provided
			ra = _as_float_scalar(ra, "ra")
			dec = _as_float_scalar(dec, "dec")
			if not -90 <= dec <= 90:
				raise ValueError(f"The value of 'dec' must be normalized to [-90,90]; was given '{dec}'.")
			if not -360 <= ra <= 360:
				raise ValueError(f"The value of 'ra' must be normalized to [0,360]; was given '{ra}'.")

		if ipix:
			ra,dec = self.ipix2ang(ipix)

		return q3c.facenum(self._hprm, ra, dec)

		# Older code:

		# # normalize	 angles
		# if not 0. <= ra < 360.:
		# 	ra = np.fmod(ra,360.)
		# 	if ra < 0.:
		# 		ra += 360.0
		#
		# # .. todo: normalize dec; remember that ra would also need to be changed
		# if not (-90 <= dec <= 90):
		# 	raise ValueError(f"Coordinate values must be normalized; ra on [0,360], dec on [-90,90].")
		#
		# if -45 < dec < 45:
		# 	# faces 1-4
		# 	if ra < 45:
		# 		return 1
		# 	elif ra < 90:
		# 		return 2
		# 	elif ra < 180:
		# 		return 3
		# 	return 4
		# else:
		# 	# north and south faces
		# 	if dec > 0:
		# 		return 0
		# 	else:
		# 		return 5

	def ipix_area(self, ipix:int, depth:int) -> float:
		'''
		Return the area of a given QLSC pixel in steradians for a given ipix.

		NOTE: currently this method returns the average bin size and ignores the ipix and depth values!

		Not all pixels are guaranteed to be the same size, but they are very close,
		i.e. the size of any one pixel is a good approximation to any other.

		pixel area ≈ 4π / no. pixels on sphere

		:param ipix: the pixel identifier number
		:returns: the area of the specified pixel in steradians
		'''

		return 4. * math.pi / self.nbins

		#if not (0 <= depth <= 30):
		#	raise ValueError(f"The depth provided {depth} is out of the range 0-30")

		# The paramter call for the underlying C function is q3c_pixarea(hprm, ipix, depth).
		# It's not clear why "depth" is a parameter since it should be read from hprm.
		# I left the parameter exposed in the Python wrapper, but in this method will just pass depth.

		#return q3c.pixarea(self._hprm, ipix, depth)

	def ang2ipix_at_depth(self, ra:float, dec:float, depth:int) -> int:
		'''
		Return the ipix value at the specified coordinate for any given bin level.
		Returns the ipix value of the pixel center at certain pixel depth covering the specified (ra,dec). [??]

		:param ra: right ascension (degrees); must be in [0,360]
		:param dec: declination (degrees); must be in [-90,90]
		:param depth: ; must be in the integer range [1-30]
		:returns:
		'''
		if not (1 <= depth <= 30):
			raise ValueError(f"The depth provided {depth} is outside the integer range [1,30].")
		return ( q3c.ang2ipix(self._hprm, ra, dec) >> (2*depth) << (2*depth) ) + (1 << (2 * (depth-1))) - 1

	def ipix2polygon(self, ipix:int=None, duplicate_endpoint=False):
		'''
		Returns the points that describe the polygon that define this pixel.

		All lines between each coordinate in the polygon should be drawn on great circles.
		The `Cornish library <http://github.com/demitri/cornish>`_ can be useful for this.

		Unlike :py:func:`ipix2ang` and :py:func:`ipix2xy`, this method takes a single ipix
		value. Vectorizing it would return a (n,4,2) array whose corners must still be
		projected face by face (:py:func:`xy2ang` takes one face number), so the shape and
		the grouping gain nothing over calling this method in a loop.

		:param ipix: the pixel identifier number
		:param duplicate_endpoint: if True, repeat the first coordinate in the polygon as the last element
		:returns: an array of ra,dec coordinates in degrees
		'''
		if ipix is None:
			raise ValueError(f"An 'ipix' value must be specified.")
		if not (0 <= ipix < self.nbins):
			raise ValueError(f"The ipix value given ('{ipix}') is out of the valid integer range for this scheme: [0,{self.nbins-1}].")

		# array of bin edges, both x and y
		#bins = np.linspace(-1,1,endpoint=True, dtype=np.double, num=self.nside*4+1)

		d = 2./self.nside # bin_length since -1 ≤ x ≤ 1 (d=delta)

		facenum, x, y = self.ipix2xy(ipix) # returns lower left corner

		if duplicate_endpoint:
			polygon = np.array([
					[  x ,  y  ],
					[ x+d,  y  ],
					[ x+d, y+d ],
					[  x , y+d ],
					[  x ,  y  ]
				], dtype=np.double)
		else:
			polygon = np.array([
					[  x ,  y  ],
					[ x+d,  y  ],
					[ x+d, y+d ],
					[  x , y+d ]
				], dtype=np.double)

		if duplicate_endpoint == True:
			polygon = polygon[:-1]

		return self.xy2ang(facenum=facenum, points=polygon)

	def point_in_polygon(self, ra:float=None, dec:float=None, polygon:Iterable=None) -> bool:
		'''
		Test whether an ra,dec point lies inside a polygon on the sphere.

		The polygon is an array of [ra,dec] vertices with shape (n,2), where 3 <= n <= 100,
		e.g. the output of :py:func:`ipix2polygon`; the edges between vertices follow
		great circles, and the winding order does not matter. The polygon must be small
		enough to fit within the projection of a single cube face; the exact limit depends
		on where the polygon sits relative to the face centers (up to ~76 degrees across
		for a face-centered polygon, and anything up to ~40 degrees across is safe
		anywhere on the sky). A ValueError is raised when the limit is exceeded.
		The result does not depend on this object's depth.

		:param ra: right ascension of the point to test (degrees)
		:param dec: declination of the point to test (degrees)
		:param polygon: an array of the polygon's vertices as [ra,dec] pairs, shape (n,2)
		:returns: True if the point is inside the polygon
		'''
		if any([x is None for x in [ra, dec, polygon]]):
			raise ValueError("Values for 'ra', 'dec', and 'polygon' must all be specified.")

		# _as_float_scalar also rejects non-finite values
		ra = _as_float_scalar(ra, "ra")
		dec = _as_float_scalar(dec, "dec")

		if not (-90 <= dec <= 90):
			ra, dec = _normalize_ang(ra, dec)

		poly_ra, poly_dec = _polygon_ra_dec(polygon)

		return q3c.check_sphere_point_in_poly(self._hprm, ra, dec, poly_ra, poly_dec)

	def point_in_ellipse(self, ra:float=None, dec:float=None, ra_cen:float=None, dec_cen:float=None,
						 semi_major_axis:Union[float, "Quantity"]=None, axis_ratio:float=None,
						 position_angle:Union[float, "Quantity"]=0.0) -> bool:
		'''
		Test whether an ra,dec point lies inside an ellipse on the sphere.

		The ellipse is the one Q3C defines: centered on (``ra_cen``, ``dec_cen``), with the
		given semi-major axis, axis ratio (semi-minor over semi-major, in ``(0,1]``, where 1
		is a circle) and position angle. This is the exact containment test used to filter
		:py:func:`QLSCIndex.ellipse_query` results, and its answer does not depend on this
		object's depth.

		:param ra: right ascension of the point to test (degrees)
		:param dec: declination of the point to test (degrees)
		:param ra_cen: right ascension of the center of the ellipse (degrees)
		:param dec_cen: declination of the center of the ellipse (degrees)
		:param semi_major_axis: semi-major axis of the ellipse (degrees), accepts a `float` value or `astropy.units.Quantity`
		:param axis_ratio: ratio of the semi-minor to the semi-major axis, in the range (0,1]
		:param position_angle: position angle of the semi-major axis (degrees), accepts a `float` value or `astropy.units.Quantity`
		:returns: True if the point is inside the ellipse
		'''
		if any([x is None for x in [ra, dec, ra_cen, dec_cen, semi_major_axis, axis_ratio]]):
			raise ValueError("Values for 'ra', 'dec', 'ra_cen', 'dec_cen', 'semi_major_axis' and 'axis_ratio' must all be specified.")

		# _as_float_scalar (via _degrees) also rejects non-finite values
		ra = _as_float_scalar(ra, "ra")
		dec = _as_float_scalar(dec, "dec")
		ra_cen = _as_float_scalar(ra_cen, "ra_cen")
		dec_cen = _as_float_scalar(dec_cen, "dec_cen")
		semi_major_axis = _degrees(semi_major_axis, "semi_major_axis")
		axis_ratio = _as_float_scalar(axis_ratio, "axis_ratio")
		position_angle = _degrees(position_angle, "position_angle")

		if not (-90 <= dec <= 90):
			ra, dec = _normalize_ang(ra, dec)
		if not (-90 <= dec_cen <= 90):
			ra_cen, dec_cen = _normalize_ang(ra_cen, dec_cen)

		# the remaining validation (axis ratio range, negative semi-major axis) is
		# performed in the C layer, which raises ValueError
		return q3c.in_ellipse(ra=ra, dec=dec, ra_cen=ra_cen, dec_cen=dec_cen,
							  semi_major_axis=semi_major_axis, axis_ratio=axis_ratio,
							  position_angle=position_angle)


class QLSCIndex:
	'''
	An object that stores (ra,dec) coordinates on a sphere and performs fast cone searches.

	The ``QLSCIndex`` object can be provided with a file path to save the index for repeated use, e.g.

	.. code-block:: python

		qlsc = QLSC(depth=30)
		sdss_catalog = QLSCIndex(qlsc=qlsc, filepath='/path/to/save/file.qlscidx')

	The filename and extension can be anything, but `.qlscidx` is suggested to remember what kind of file it is.
	If no ``filepath`` parameter is provided, an index is created in memory, but note that it will
	no longer exist after the program ends. If an existing file is found at the file path provided,
	an attempt will be made to open it and, if a valid ``QLSCIndex`` file, it will be used.
	The file path can be provided in the form of a URI, e.g. ``file:///Users/meerkat/catalog.qlscidx``.

	The underlying storage is an SQLite database, but this should be considered an implementation
	detail and can change. (This may be documented for expanded use in the future.)

	:param qlsc: an instance of :py:class:`QLSC` set to the desired segmentation level
	:param filepath: the file path (string path, ``pathlib.Path``, or ``file://`` URL) where a persistent index will be stored; the default is to create the index in memory
	:param name: a user-provided name for the index
	:param description: a user-provided description of the index
	'''
	def __init__(self, qlsc:QLSC=QLSC(depth=30), filepath:Union[os.PathLike, str]=None, name:str=None, description:str=None):

		# .. todo:: support 'overwrite' -> check if file already exists

		#if qlsc is None:
		#	# select the highest resolution by default
		#	qslc = QLSC(depth=30)

		if not isinstance(qlsc, QLSC):
			raise Exception(f"The 'qlsc' parameter is not a QLSC object as expected.")

		self.qlsc = qlsc
		self.name = name
		self.description = description
		self.memory_db_connection = None # this must persist
		self._closed = False
		self._query_qlsc = None # depth-30 QLSC used for region queries; created on demand
		#self._db = None # SQLite database connection, defined in _initial_database_connection

		self.database_tablename = "qlsc_ipix"
		self.qlsc_schema_version = "1"

		if filepath:
			if isinstance(filepath, pathlib.PurePath):
				self.db_filepath = filepath
			else:
				if filepath.startswith("file:"): # is URI
					scheme, netloc, path, params, query, fragment = urlparse(filepath)
					self.db_filepath = path
				else:
					self.db_filepath = pathlib.Path(filepath)

			if self.db_filepath.exists():
				logger.debug(f"Found existing db at '{self.db_filepath}'.")
			else:
				if os.access(self.db_filepath.parent, os.W_OK) == False:
					raise Exception(f"The filepath provided was not found, but the location cannot be written to: '{filepath}'")
				logger.debug(f"Provided db filepath ('{self.db_filepath}') not found; creating new db.")
		else:
			# interpret "None" to use an in memory database

			# allows multiple connections to the in-memory database
			# see: https://sqlite.org/inmemorydb.html
			self.db_filepath = ":memory:" # this just keeps creating a file on disk -> "file::memory:?cache=shared"

		self._initial_database_connection()

	def _initial_database_connection(self):
		'''
		Make the initial connection to the database.
		'''

		is_new_database = True

		#if self._is_in_memory_db():
		if (self.db_filepath == ":memory:") or str(self.db_filepath).startswith("file::memory:"):
			# Always keep this connection open to make sure the database is not released
			# i.e., set it here and never touch it again.
			# Use "self._db" to open/close connections as one would with a file-based db.
			#
			# ^--- shared cache plan, but didn't get that to work
			#
			#self._in_memory_db_connection = sqlite3.connect(self.index_path)
			connection = sqlite3.connect(self.db_filepath) #opens a second connection that can be closed without losing the database
			self.memory_db_connection = connection
		else:

			#is_uri = self.index_path.startswith("file:")
			#if is_uri:
			#	scheme, netloc, path, params, query, fragment = urlparse(self.index_path)
			#	self.index_path = path

			is_new_database = not self.db_filepath.exists()
			#logger.debug(f"index path = {self.index_path}")

			try:
				connection = sqlite3.connect(self.db_filepath, timeout=20) #, uri=True)
			except sqlite3.OperationalError as e:
				if is_new_database:
					raise Exception(f"Unable to create database at specified path ('{self.db_filepath}').")
				else:
					raise Exception(f"Found file at path '{self.db_filepath}', but am unable to open as an SQLite database.")

		self._configure_db_connection(connection)

		if is_new_database:
			self._init_sqlite_db(connection)
		else:
			# Opening existing database
			db_ok = self._validate_db_as_qlsc_index(connection)
			if not db_ok:
				connection.close()
				raise Exception(f"The file provided is an SQLite database (good) but not in a format I recognize (bad).")

		if not self._is_in_memory_db():
			connection.close()

	def _validate_db_as_qlsc_index(self, connection:sqlite3.Connection):
		'''
		Validates that the database found is one of ours.
		'''
		# keep this simple for the time being
		try:
			with contextlib.closing(connection.cursor()) as cursor:
				metadata = cursor.execute("SELECT depth, index_name, index_description FROM qlsc_metadata").fetchone()
				if metadata is None:
					raise RuntimeError(f"The file at '{self.db_filepath}' has a QLSC index metadata table but no metadata row; the file appears to be corrupted.")
				if self.qlsc.depth != metadata["depth"]:
					self.qlsc = QLSC(depth=metadata["depth"])
				self.name = metadata["index_name"]
				self.description = metadata["index_description"]
				return True
		except sqlite3.OperationalError as e:
			if "no such table" in str(e):
				return False
			# any other database error (locked, malformed, not a database, ...) is a
			# real problem, not evidence that this isn't a QLSC index - surface it
			raise

	def close(self):
		'''
		Close the index.

		For in-memory indices this releases the underlying database, discarding the data;
		file-based indices are unaffected on disk. The index cannot be used after being
		closed. This method is idempotent, and the object can also be used as a context
		manager, e.g. ``with QLSCIndex(...) as index:``.
		'''
		if self.memory_db_connection is not None:
			self.memory_db_connection.close()
			self.memory_db_connection = None
		self._closed = True

	def __enter__(self):
		return self

	def __exit__(self, exc_type, exc_value, traceback):
		self.close()

	def __del__(self):
		''' Destructor method. '''
		# guard with hasattr in case __init__ raised before defining the attribute
		if hasattr(self, "memory_db_connection"):
			self.close()

	def _ensure_open(self):
		''' Raise if this index has been closed; call at the start of every public operation. '''
		if self._closed:
			raise RuntimeError("This QLSCIndex has been closed and can no longer be used.")

	def _is_in_memory_db(self):
		#return (self.db_filepath == ":memory:") or self.db_filepath.startswith("file::memory:")
		return self.memory_db_connection is not None

	def _configure_db_connection(self, connection:sqlite3.Connection):
		'''
		Configure connection-level settings on the SQLite database.
		'''
		# set database-specific settings
		connection.isolation_level = None    # autocommit mode; transactions can be explicitly created with BEGIN/COMMIT statements
		connection.row_factory = sqlite3.Row # return dictionaries instead of tuples from SELECT statements

	def _init_sqlite_db(self, connection:sqlite3.Connection):
		'''
		Initialize a new index database, e.g. create schema, initialize metadata.

		Methods that accept a 'connection' parameter are not responsible for closing it.

		:param connection:
		'''

		#self.sqlite_db = sqlite3.connect(':memory:')
		#con.create_function("ang2ipix", 1, self.ang2ipix)

		#self.data_source.create_function("qlsc_ang2ipix", 3, self.qlsc.ang2ipix)

# 		def plus(a,b):
# 			return a+b

		#self.data_source.create_function("qlsc_ang2ipix", 2, plus)

		with contextlib.closing(connection.cursor()) as cursor:

			cursor.execute(f'''
				CREATE TABLE {self.database_tablename} (
                    ipix INTEGER,
                    ra REAL,
                    dec REAL,
                    key TEXT,
                    UNIQUE (ra, dec)
                );''')

			cursor.execute(f'''CREATE INDEX ipix_idx ON {self.database_tablename}(ipix); ''')

            # Create a metadata table - COUNT(*) on SQLite databases always does a full scan,
            # so keep track of rows count here.
            # Could add additional information here.
			cursor.execute(f'''
            	CREATE TABLE qlsc_metadata (
	            	row_count INTEGER,
	            	qlsc_version TEXT,
	            	depth INTEGER,
	            	date_created TEXT,
	            	index_name TEXT,
	            	index_description TEXT
	            );''')

			# set initial metadata values -> CURRENT_TIMESTAMP might be used instead of strftime
			cursor.execute(f'''
				INSERT INTO qlsc_metadata
					(row_count, qlsc_version, depth, date_created, index_name, index_description) VALUES
					(0, ?, ?, strftime('%Y-%m-%dT%H:%M:%SZ','now'), ?, ?);
				''', (self.qlsc_schema_version, self.qlsc.depth, self.name, self.description))

			# internally keep count of records since a count(*)
			# in SQLite is a full table scan
# 			cursor.execute(f'''
# 				CREATE TRIGGER row_count_trigger
# 					AFTER INSERT ON {self.database_tablename}
# 					BEGIN
# 						UPDATE qlsc_metadata SET row_count = row_count + 1;
# 					END;
# 				''');

		# Close and reopen the database to reset the "connection.total_changes" value to 0,
		# used to count number of changed rows after inserts.
		# If this is an in-memory database, open a second connection first to make sure that
		# the memory isn't freed (happens when last connection is closed).
# 		if self._is_in_memory_db:
# 			tmp_var = sqlite3.connect(self.index_path)
# 			self._db.close()
# 			self._db = tmp_var
# 		else:
		#self._db.close()
		#self._db = sqlite3.connect(self.index_path)
		#self._configure_db_connection(self._db)

	def _new_db_connection(self) -> sqlite3.Connection:
		'''
		Open and return a new db connection.

		This class minimizes the time a connection is open to the database so
		that it will work in a multithreaded/multiprocessing environment.

		The exception is for in-memory databases where the connection is kept
		open for the life of this object since the database is deleted when the
		last connection is closed.

		It is the responsibility of the calling function to close the connection when finished.
		This is not done automatically in case the object is an in-memory database.
		'''
		# without this check, an operation on a closed in-memory index would
		# silently open a brand new, empty database
		self._ensure_open()
		if self.memory_db_connection is not None:
			return self.memory_db_connection
		connection = sqlite3.connect(self.db_filepath, timeout=20)
		# match the in-memory connection settings (row factory, autocommit)
		self._configure_db_connection(connection)
		return connection

	def add_point(self, ra:float=None, dec:float=None, key:str=None):
		'''
		Add a single coordinate point to the index.

		This method can be called multiple times with the same coordinate, but duplicate
		entries will be discarded. Key values are not required to be unique. If a coordinate
		is added to the index without a key and again later with a key, it will not be added
		to the index (and vice versa).

		:param ra: a single right ascension value, in degrees
		:param dec: a single declination value, in degrees
		:param key: a unique identifier for the coordinate (optional)
		'''

		self._ensure_open()

		if any([x is None for x in [ra,dec]]):
			raise ValueError("Both 'ra' and 'dec' values must be specified.")

		# Coerce scalar-like values: without this, a NumPy array passed here would be
		# stored by sqlite3 as a BLOB via the buffer protocol - silently losing the
		# point and corrupting later queries on the index.
		ra = _as_float_scalar(ra, "ra")
		dec = _as_float_scalar(dec, "dec")

		# Q3C.ang2ipix truncates dec outside of [-90,90]
		if not (-90 <= dec <= 90):
			ra, dec = _normalize_ang(ra,dec)

		ipix = self.qlsc.ang2ipix(ra,dec)

		connection = self._new_db_connection()

		with contextlib.closing(connection.cursor()) as cursor:

			cursor.execute("BEGIN")
			if key is None:
				query = f"INSERT OR IGNORE INTO {self.database_tablename} (ipix, ra, dec) VALUES (?, ?, ?)"
				cursor.execute(query, (ipix, ra, dec)) # pass list of tuples (e.g. use zip)
			else:
				query = f"INSERT OR IGNORE INTO {self.database_tablename} (ipix, ra, dec, key) VALUES (?, ?, ?, ?)"
				cursor.execute(query, (ipix, ra, dec, key))
			cursor.execute("COMMIT")

		if not self._is_in_memory_db():
			connection.close()


	def add_points(self, ra:Iterable=None, dec:Iterable=None, points:Iterable[Tuple]=None, keys:Iterable[str]=None):
		'''
		Add a collection of points to the index.

		This method can be called multiple times with the same coordinate, but duplicate
		entries will be discarded. Key values are not required to be unique. If a coordinate
		is added to the index without a key and again later with a key, it will not be added
		to the index (and vice versa).

		This method is preferred for a large number of points over using :py:func:`add_point` in a loop.

		Note that if the ``dec`` values are not in the range [-90,90], a copy of the ra,dec arrays will
		be made to normalize the points first. This may only be an issue when providing an extremely
		large number of points.

		:param ra: an Iterable (e.g. NumPy array) of right ascension points in degrees
		:param dec: an Iterable of declination points in degrees
		:param points: an Iterable of tuples of (ra,dec) pairs
		:param key: an Iterable of unique identifiers for each point
		'''
		self._ensure_open()

		if all([x is None for x in [ra,dec,points]]):
			raise ValueError("'ra','dec' OR 'points' must be specified.")
		if all([x is True for x in [ra, dec, points]]):
			raise ValueError("Only 'ra','dec' OR 'points' can be specified.")
		elif points is None and any([x is None for x in [ra,dec]]):
			raise ValueError("If 'ra' or dec' is given, then other parameter must also be provided.")

		if keys is not None:
			if points is not None:
				if len(points) != len(keys):
					raise ValueError(f"The number of points provided ({len(points)}) does not match the number of keys ({len(keys)}).")
			else:
				if len(ra) != len(dec) != len(keys):
					raise ValueError(f"The number of ra,dec coordinates provided ({len(ra)},{len(dec)}) does not match the number of keys ({len(keys)}).")

#		if pandas_available:
#			if isinstance(ra, pandas.core.series.Series):
#				ra = ra.values # extract ndarray
#			if isinstance(dec, pandas.core.series.Series):
#				dec = dec.values # extract ndarray
#			if isinstance(key, pandas.core.series.Series):
#				key = key.values # extract ndarray

		if points is not None:
			# break it down
			ra = points[:,0]
			dec = points[:,1]

		ra = np.asarray(ra, dtype=np.double)
		dec = np.asarray(dec, dtype=np.double)

		if not (np.isfinite(ra).all() and np.isfinite(dec).all()):
			raise ValueError("The 'ra' and 'dec' values must all be finite.")

		if not ((-90 <= dec).all() and (dec <= 90).all()):
			# normalize BOTH the ipix calculation inputs and the stored coordinates
			ra, dec = _normalize_ang(ra=ra, dec=dec, copy=True)

		connection = None
		try:
			ipix = self.qlsc.ang2ipix(ra, dec)

			connection = self._new_db_connection()

			with contextlib.closing(connection.cursor()) as cursor:

				cursor.execute("BEGIN")
				if keys is None:
					query = f"INSERT OR IGNORE INTO {self.database_tablename} (ipix, ra, dec) VALUES (?, ?, ?)"
					# use individual lists instead of something like np.dstack- we don't want
					# the 64-bit ipix integers converted to floats and risk changing the value due to rounding
					cursor.executemany(query, zip(ipix, ra, dec)) # pass list of tuples (e.g. use zip)
				else:
					query = f"INSERT OR IGNORE INTO {self.database_tablename} (ipix, ra, dec, key) VALUES (?, ?, ?, ?)"
					cursor.executemany(query, zip(ipix, ra, dec, keys))

				cursor.execute("COMMIT")

		except sqlite3.IntegrityError as e:
			raise NotImplementedError(f"Unexpected integrity error while adding points: {e}") from e
		finally:
			if connection is not None and not self._is_in_memory_db():
				connection.close()

	@property
	def number_of_points(self) -> int:
		''' Return the number of coordinate points in the index. '''

		connection = self._new_db_connection()

		try:
			with contextlib.closing(connection.cursor()) as cursor:
				cursor.execute(f"SELECT count(*) FROM {self.database_tablename}")
				return cursor.fetchone()[0]
		finally:
			if not self._is_in_memory_db():
				connection.close()

	def remove_point(self, ra:float=None, dec:float=None, key:str=None) -> int:
		'''
		Remove points from the index, identified either by their key or by their exact coordinate.

		Provide either ``key`` or both ``ra`` and ``dec``, not both forms. Keys are not
		required to be unique, so removing by key removes every point carrying that key;
		coordinates are unique in the index, so removing by coordinate removes at most one
		point. The coordinate given must be the one that was stored: as with
		:py:func:`add_point`, a declination outside [-90,90] is normalized first, but no
		tolerance is applied — use :py:func:`radial_query` to find a nearby point and remove
		the coordinate it reports.

		Removing something that is not in the index raises ``KeyError`` rather than quietly
		doing nothing.

		:param ra: right ascension of the point to remove, in degrees
		:param dec: declination of the point to remove, in degrees
		:param key: the key of the point(s) to remove
		:returns: the number of points removed
		'''
		self._ensure_open()

		by_coordinate = ra is not None or dec is not None

		if key is None and not by_coordinate:
			raise ValueError("A 'key' or an 'ra','dec' coordinate must be given to identify the point(s) to remove.")
		if key is not None and by_coordinate:
			raise ValueError("Specify either 'key' or 'ra','dec' to identify the point(s) to remove, not both.")

		if by_coordinate:
			if ra is None or dec is None:
				raise ValueError(f"Both 'ra' and 'dec' must be specified when removing by coordinate (ra={ra}, dec={dec}).")
			ra = _as_float_scalar(ra, "ra")
			dec = _as_float_scalar(dec, "dec")
			# match add_point: the stored coordinate is the normalized one
			if not (-90 <= dec <= 90):
				ra, dec = _normalize_ang(ra, dec)
			query = f"DELETE FROM {self.database_tablename} WHERE ra=? AND dec=?"
			parameters = (ra, dec)
			description = f"the coordinate (ra={ra}, dec={dec})"
		else:
			key = self._validated_key(key)
			query = f"DELETE FROM {self.database_tablename} WHERE key=?"
			parameters = (key,)
			description = f"the key '{key}'"

		connection = self._new_db_connection()
		try:
			with contextlib.closing(connection.cursor()) as cursor:
				cursor.execute("BEGIN")
				cursor.execute(query, parameters)
				removed = cursor.rowcount # must be read before the COMMIT resets it
				cursor.execute("COMMIT")
		finally:
			if not self._is_in_memory_db():
				connection.close()

		if removed == 0:
			raise KeyError(f"No point matching {description} was found in this index.")

		return removed

	def _validated_key(self, key) -> str:
		'''
		Validate a key used to look up or remove points.

		The ``key`` column has SQLite's TEXT affinity, so a non-string key would be
		stored (and compared) as its text representation — a silent mismatch waiting
		to happen. Require a string and say so instead.
		'''
		if not isinstance(key, str):
			raise TypeError(f"Index keys are strings (the underlying column has TEXT affinity); "
							f"was given '{key!r}' of type '{type(key).__name__}'.")
		return key

	def keys(self) -> Iterable[str]:
		'''
		Return a generator over the keys of the points in this index.

		Keys are not required to be unique and a point may have been added without one.
		A key shared by several points is therefore yielded once per point, and points
		stored without a key are **not** represented at all — so the number of keys
		yielded can be either larger or smaller than ``len(index)``. Compare the two if
		the difference matters. No particular ordering is guaranteed.

		The generator holds a database connection open while it is being consumed
		(file-backed indices only); the connection is released when the generator is
		exhausted or closed.

		:returns: a generator yielding the key of each keyed point in the index
		'''
		self._ensure_open() # fail here, not on the first next() call
		return self._iterate_keys()

	def _iterate_keys(self):
		''' The generator returned by :py:func:`keys`. '''
		connection = self._new_db_connection()
		try:
			with contextlib.closing(connection.cursor()) as cursor:
				for row in cursor.execute(f"SELECT key FROM {self.database_tablename} WHERE key IS NOT NULL"):
					yield row["key"]
		finally:
			if not self._is_in_memory_db():
				connection.close()

	def __len__(self) -> int:
		'''
		Return the number of coordinate points in the index.

		Note that this makes an empty index falsy, as an empty dict or list is.
		'''
		return self.number_of_points

	def __iter__(self) -> Iterable[str]:
		''' Iterate over the keys in the index; see :py:func:`keys`. '''
		return self.keys()

	def __contains__(self, key:str) -> bool:
		'''
		Return whether any point in the index carries the given key.

		:param key: the key to look for
		'''
		self._ensure_open()
		key = self._validated_key(key)

		connection = self._new_db_connection()
		try:
			with contextlib.closing(connection.cursor()) as cursor:
				cursor.execute(f"SELECT 1 FROM {self.database_tablename} WHERE key=? LIMIT 1", (key,))
				return cursor.fetchone() is not None
		finally:
			if not self._is_in_memory_db():
				connection.close()

	def _query_scheme(self) -> QLSC:
		'''
		Return the depth-30 :class:`QLSC` used to compute candidate ipix ranges.

		Every region query (cone, polygon, ellipse) is computed at depth 30 and the
		resulting ranges mapped down to the index's depth: the underlying Q3C code
		produces degenerate (empty) ranges for some queries at small nside values,
		since upstream Q3C is only ever exercised at depth 30.
		'''
		if self._query_qlsc is None:
			self._query_qlsc = self.qlsc if self.qlsc.depth == 30 else QLSC(depth=30)
		return self._query_qlsc

	def _index_depth_ranges(self, fulls:np.ndarray, partials:np.ndarray) -> list:
		'''
		Map the depth-30 candidate ipix ranges returned by the C layer down to this
		index's depth, then merge them.

		The half-open ranges are widened, never narrowed: the upper bound is rounded
		up so a pixel of the index that is only partly covered by a depth-30 range is
		still a candidate. Q3C's output contains duplicate and overlapping ranges,
		especially after mapping down to lower depths, so the result is sorted and
		merged into disjoint ranges.

		:param fulls: the fully covered depth-30 ipix ranges, shape (n,2)
		:param partials: the partially covered depth-30 ipix ranges, shape (n,2)
		:returns: a sorted list of disjoint [min_ipix, max_ipix) pairs at the index depth
		'''
		# number of ipix bits to shift from depth 30 down to the index depth
		k = 2 * (30 - self.qlsc.depth)

		ranges = list()
		for min_ipix,max_ipix in np.vstack((fulls, partials)):
			if k > 0:
				min_ipix = min_ipix >> k
				max_ipix = (max_ipix + (1 << k) - 1) >> k # round up so partially covered pixels are kept
			if min_ipix < max_ipix: # drop empty [a,a) ranges
				ranges.append([int(min_ipix), int(max_ipix)])

		ranges.sort()
		merged_ranges = list()
		for lo,hi in ranges:
			if merged_ranges and lo <= merged_ranges[-1][1]:
				merged_ranges[-1][1] = max(merged_ranges[-1][1], hi)
			else:
				merged_ranges.append([lo, hi])

		return merged_ranges

	def _candidate_row_query(self, merged_ranges:list, columns:str, description:str) -> str:
		'''
		Build the SQL that selects the candidate rows for a region query.

		If the region produced no usable ipix ranges (which should not happen), every
		row is selected rather than silently returning nothing; the caller's exact
		filter still gives the right answer, only more slowly, and the fallback is
		logged as a warning.

		:param merged_ranges: disjoint [min_ipix, max_ipix) pairs at the index depth
		:param columns: the SQL column list to select
		:param description: the query, for the warning message
		:returns: the SQL query string
		'''
		if merged_ranges:
			wheres = " OR ".join(f"(ipix>={lo} AND ipix<{hi})" for lo,hi in merged_ranges)
			return f"SELECT {columns} FROM {self.database_tablename} WHERE ({wheres})"

		logger.warning(f"The {description} produced no candidate ipix ranges; falling back to scanning the entire index.")
		return f"SELECT {columns} FROM {self.database_tablename}"

	def _select_candidate_rows(self, query:str):
		'''
		Run a candidate row query and yield its rows.

		The connection is opened and closed around the iteration, matching the
		per-operation connection handling used everywhere else in this class.

		:param query: the SQL produced by :py:func:`_candidate_row_query`
		:returns: a generator of :class:`sqlite3.Row`
		'''
		connection = self._new_db_connection()
		try:
			with contextlib.closing(connection.cursor()) as cursor:
				try:
					for row in cursor.execute(query):
						yield row
				except sqlite3.OperationalError as e:
					raise RuntimeError(f"Error in accessing database in region query.\n\nQuery: \"{query}\"\n\nError: {e}") from e
		finally:
			if not self._is_in_memory_db():
				connection.close()

	def _find_candidate_ranges(self, ra:float, dec:float, radius:float):
		'''
		Return the (fulls, partials) candidate ipix range arrays covering a cone search,
		always computed at depth 30.

		The underlying q3c_radial_query produces degenerate (empty) ipix ranges for some
		queries at small nside values - upstream Q3C is only ever exercised at depth 30 -
		so the search always runs at depth 30; the caller maps the ranges down to the
		index's resolution.

		For zero and very small (< ~1e-7 deg) radii, q3c produces no ranges at all at
		some sky positions, and at seam/pole positions its ranges can omit the pixels
		holding equivalent coordinate representations (e.g. a point stored as ra=360
		queried as ra=0), so the search uses a radius of at least 1e-6 degrees. At rare
		positions a given search radius can also overflow q3c's internal range capacity;
		a larger radius produces fewer, coarser ranges, so the radius is grown and the
		search retried on either failure mode, raising only if even a whole-sphere
		search fails.

		The ranges therefore cover at least the requested cone, possibly more; the
		caller MUST filter the candidate rows by the radius actually requested.
		'''
		query_hprm = self._query_scheme()._hprm

		return self._widened_candidate_ranges(lambda size: q3c.radial_query(query_hprm, ra, dec, size),
											  radius,
											  f"radial query at (ra={ra}, dec={dec}) with radius={radius} deg")

	def _find_candidate_ellipse_ranges(self, ra:float, dec:float, semi_major_axis:float,
									   axis_ratio:float, position_angle:float):
		'''
		Return the (fulls, partials) candidate ipix range arrays covering an ellipse
		search, always computed at depth 30.

		The same failure modes as :py:func:`_find_candidate_ranges` apply, for the same
		reasons, and the same remedy works: growing the semi-major axis while holding the
		axis ratio and position angle fixed scales both axes up about the same center, so
		the widened ellipse contains the requested one.

		The ranges therefore cover at least the requested ellipse, possibly more; the
		caller MUST filter the candidate rows by the ellipse actually requested.

		:param ra: right ascension of the center of the ellipse (degrees)
		:param dec: declination of the center of the ellipse (degrees)
		:param semi_major_axis: semi-major axis of the ellipse (degrees)
		:param axis_ratio: ratio of the semi-minor to the semi-major axis
		:param position_angle: position angle of the semi-major axis (degrees)
		:returns: a tuple of the (fulls, partials) depth-30 ipix range arrays
		'''
		query_hprm = self._query_scheme()._hprm

		return self._widened_candidate_ranges(
			lambda size: q3c.ellipse_query(query_hprm, ra, dec, size, axis_ratio, position_angle),
			semi_major_axis,
			f"ellipse query at (ra={ra}, dec={dec}) with semi_major_axis={semi_major_axis} deg, "
			f"axis_ratio={axis_ratio}, position_angle={position_angle} deg")

	def _widened_candidate_ranges(self, region_query, size:float, description:str):
		'''
		Run a candidate range search, growing the region until it produces usable ranges.

		For zero and very small (< ~1e-7 deg) regions, q3c produces no ranges at all at
		some sky positions, and at seam/pole positions its ranges can omit the pixels
		holding equivalent coordinate representations (e.g. a point stored as ra=360
		queried as ra=0), so the search uses a size of at least 1e-6 degrees. At rare
		positions a given size can also overflow q3c's internal range capacity; a larger
		region produces fewer, coarser ranges, so the size is grown and the search retried
		on either failure mode, raising only if even a whole-sphere search fails.

		The region searched always contains the region asked for, so the ranges are a
		superset and the caller's exact filter still gives the right answer.

		:param region_query: a callable taking a size in degrees and returning the (fulls, partials) arrays
		:param size: the requested size of the region in degrees (a radius or a semi-major axis)
		:param description: the query, for the warning message
		:returns: a tuple of the (fulls, partials) depth-30 ipix range arrays
		'''
		initial_search_size = max(size, 1e-6)
		search_size = initial_search_size
		while True:
			try:
				fulls, partials = region_query(search_size)
			except q3c.RangeOverflowError:
				if search_size >= 180.0:
					raise # even a whole-sphere search failed
				search_size = min(search_size * 10, 180.0)
				continue
			if len(fulls) + len(partials) > 0 or search_size >= 180.0:
				break
			search_size = min(search_size * 10, 180.0)

		if search_size > initial_search_size:
			logger.warning(f"The {description} required widening the candidate search region to {search_size} deg; results are correct, but the query may be slower than usual.")

		return fulls, partials

	def radial_query(self, ra:float, dec:float, radius:Union[float, "Quantity"], return_key:bool=False, return_distance:bool=False) -> np.recarray:
		'''
		Given an ra,dec coordinate and a radius (all in degrees), return the points that fall in the cone search.

		The default return value (both flags left False) is unchanged: an array of the
		matching coordinates with shape (n,2). If either ``return_key`` or
		``return_distance`` is set, the result is instead a :class:`numpy.recarray` whose
		fields are "ra", "dec", then "key" and "distance" for whichever flags were set.

		The "distance" field is the angular separation between the query center and each
		match, in degrees. It costs nothing to ask for: the cone filter already computes
		the sine-squared half-distance of every candidate. Coordinates that are equivalent
		representations of the query center (e.g. a point stored at ra=360 and queried at
		ra=0, or any ra at a pole) report a separation of exactly 0.

		:param ra: right ascension (degrees)
		:param dec: declination (degrees)
		:param radius: radius (degrees) in the range [0,180], accepts a `float` value or `astropy.units.Quantity`
		:param return_key: if set to True, include the key value as provided when added to the index
		:param return_distance: if set to True, include the angular separation from the query center in degrees
		:returns: an array of matches with shape (n,2) if both flags are False; otherwise a :class:`numpy.recarray` with the fields ``ra,dec`` plus ``key`` and/or ``distance``
		'''

		self._ensure_open()

		if astropy_available:
			if isinstance(radius, Quantity):
				radius = radius.to(u.deg).value

		ra = _as_float_scalar(ra, "ra")
		dec = _as_float_scalar(dec, "dec")
		radius = _as_float_scalar(radius, "radius")

		if abs(dec) > 90:
			raise ValueError(f"The value for dec must be in the range [-90,90]; was given '{dec}'.")

		# above 180 degrees the sin^2(radius/2) filter below becomes periodic and meaningless
		if not (0 <= radius <= 180):
			raise ValueError(f"The value for radius must be in the range [0,180]; was given '{radius}'.")

		center_ra = ra
		center_dec = dec

		# gather matches into these lists
		match_ra = list()
		match_dec = list()
		match_key = list()
		match_distance = list()

		#if isinstance(self._db, sqlite3.Connection):

		fulls, partials = self._find_candidate_ranges(center_ra, center_dec, radius)

		#logger.debug(f"{center_ra=}, {center_dec=}, {radius=}")
		#logger.debug(f"{fulls=}, {partials=}")

		merged_ranges = self._index_depth_ranges(fulls, partials)

		columns = "ra,dec,key" if return_key else "ra,dec"
		# If the candidate search produced no usable ranges (which should not happen),
		# every row is scanned rather than silently missing points; the exact filter
		# below preserves correctness either way.
		query = self._candidate_row_query(merged_ranges, columns,
										  f"radial query at (ra={center_ra}, dec={center_dec}) with radius={radius} deg")

		cone_radius = pow(sin(deg2rad(radius)/2.), 2)
		center_canonical = _canonical_position(center_ra, center_dec)

		for row in self._select_candidate_rows(query):
			row_ra, row_dec = row["ra"], row["dec"]
			# Filter out points outside radius (inclusive, so e.g. radius=180 includes
			# an antipodal point). Canonically equal representations of the center
			# (e.g. ra=360 vs ra=0, or any ra at a pole) always match: sindist between
			# equivalent representations leaves tiny nonzero floating point residuals,
			# which would otherwise exclude them for zero or near-zero radii.
			if _canonical_position(row_ra, row_dec) == center_canonical:
				# the same point on the sphere, whatever the residual says
				is_match, row_sindist = True, 0.0
			elif radius > 0:
				row_sindist = sindist(row_ra, row_dec, center_ra, center_dec)
				is_match = row_sindist <= cone_radius
			else:
				# a zero-radius query matches only canonically equal points;
				# sindist is not needed and deliberately not computed
				is_match, row_sindist = False, None
			if is_match:
				match_ra.append(row_ra)
				match_dec.append(row_dec)
				if return_key:
					match_key.append(row["key"])
				if return_distance:
					match_distance.append(_sindist2deg(row_sindist))

		return _match_result(match_ra, match_dec, match_key, match_distance, return_key, return_distance)

	def polygon_query(self, polygon:Iterable=None, return_key:bool=False) -> np.ndarray:
		'''
		Return the points of the index that fall inside a polygon on the sphere.

		The polygon is an array of [ra,dec] vertices with shape (n,2), in the same form
		:py:func:`QLSC.point_in_polygon` accepts (and, for example, the form
		:py:func:`QLSC.ipix2polygon` returns); the edges between vertices follow great
		circles, the winding order does not matter, and the first vertex must not be
		repeated as the last. Both convex and concave polygons are supported.

		The polygon must be small enough to project onto a single cube face — in practice
		up to roughly 23 degrees across — and a ``ValueError`` is raised when it is not.
		This is a tighter limit than :py:func:`QLSC.point_in_polygon`, which only has to
		project the polygon far enough to place one point in it.

		Matches are found by selecting the pixels covering the polygon and then testing
		every candidate with the same exact containment test :py:func:`QLSC.point_in_polygon`
		uses, so the result contains exactly the indexed points inside the polygon. (This
		is what makes concave polygons safe: the coverage step is allowed to be generous,
		and for a concave polygon it is — Q3C calls a pixel of the coarse search fully
		covered when all four of its corners are inside the polygon, which for a concave
		polygon can include area outside it. That produces extra candidates, never missing
		ones, and the extras are removed by the exact test.)

		The default return value is an array of the matching coordinates with shape (n,2).
		If ``return_key`` is set the result is instead a :class:`numpy.recarray` with the
		fields "ra", "dec" and "key". There is no ``return_distance``: a polygon has no
		center to measure a separation from.

		:param polygon: an array of the polygon's vertices as [ra,dec] pairs, shape (n,2)
		:param return_key: if set to True, include the key value as provided when added to the index
		:returns: an array of matches with shape (n,2), or a :class:`numpy.recarray` with the fields ``ra,dec,key``
		'''
		self._ensure_open()

		if polygon is None:
			raise ValueError("A value for 'polygon' must be specified.")

		poly_ra, poly_dec = _polygon_ra_dec(polygon)

		match_ra = list()
		match_dec = list()
		match_key = list()

		# The candidate ranges are always computed at depth 30; q3c_poly_query raises
		# ValueError (via the C layer) when the polygon is too large to project.
		fulls, partials = q3c.poly_query(self._query_scheme()._hprm, poly_ra, poly_dec)

		merged_ranges = self._index_depth_ranges(fulls, partials)

		columns = "ra,dec,key" if return_key else "ra,dec"
		query = self._candidate_row_query(merged_ranges, columns,
										  f"polygon query over a polygon with {len(poly_ra)} vertices")

		# the C function is called directly rather than through QLSC.point_in_polygon:
		# the vertices have already been validated and converted above, and doing that
		# again for every candidate row is the whole cost of this loop. The test itself
		# is identical.
		query_hprm = self._query_scheme()._hprm
		for row in self._select_candidate_rows(query):
			row_ra, row_dec = row["ra"], row["dec"]
			if q3c.check_sphere_point_in_poly(query_hprm, row_ra, row_dec, poly_ra, poly_dec):
				match_ra.append(row_ra)
				match_dec.append(row_dec)
				if return_key:
					match_key.append(row["key"])

		return _match_result(match_ra, match_dec, match_key, None, return_key, False)

	def ellipse_query(self, ra:float=None, dec:float=None, semi_major_axis:Union[float, "Quantity"]=None,
					  axis_ratio:float=None, position_angle:Union[float, "Quantity"]=0.0,
					  return_key:bool=False, return_distance:bool=False) -> np.ndarray:
		'''
		Return the points of the index that fall inside an ellipse on the sphere.

		The ellipse is centered on (``ra``, ``dec``) with the given semi-major axis, axis
		ratio (semi-minor over semi-major, in ``(0,1]``) and position angle. An
		``axis_ratio`` of 1 is a circle and returns the same points as
		:py:func:`radial_query` with that radius.

		Matches are found by selecting the pixels covering the ellipse and then testing
		every candidate with :py:func:`QLSC.point_in_ellipse`, the exact containment test,
		so the result contains exactly the indexed points inside the ellipse.

		The default return value is an array of the matching coordinates with shape (n,2).
		If either ``return_key`` or ``return_distance`` is set, the result is instead a
		:class:`numpy.recarray` whose fields are "ra", "dec", then "key" and "distance" for
		whichever flags were set. The "distance" field is the angular separation between
		the center of the ellipse and each match, in degrees — not a distance to the
		ellipse's boundary.

		:param ra: right ascension of the center of the ellipse (degrees)
		:param dec: declination of the center of the ellipse (degrees)
		:param semi_major_axis: semi-major axis of the ellipse (degrees), in the range [0,180], accepts a `float` value or `astropy.units.Quantity`
		:param axis_ratio: ratio of the semi-minor to the semi-major axis, in the range (0,1]
		:param position_angle: position angle of the semi-major axis (degrees), accepts a `float` value or `astropy.units.Quantity`
		:param return_key: if set to True, include the key value as provided when added to the index
		:param return_distance: if set to True, include the angular separation from the center in degrees
		:returns: an array of matches with shape (n,2) if both flags are False; otherwise a :class:`numpy.recarray` with the fields ``ra,dec`` plus ``key`` and/or ``distance``
		'''
		self._ensure_open()

		if any([x is None for x in [ra, dec, semi_major_axis, axis_ratio]]):
			raise ValueError("Values for 'ra', 'dec', 'semi_major_axis' and 'axis_ratio' must all be specified.")

		ra = _as_float_scalar(ra, "ra")
		dec = _as_float_scalar(dec, "dec")
		semi_major_axis = _degrees(semi_major_axis, "semi_major_axis")
		axis_ratio = _as_float_scalar(axis_ratio, "axis_ratio")
		position_angle = _degrees(position_angle, "position_angle")

		if abs(dec) > 90:
			raise ValueError(f"The value for dec must be in the range [-90,90]; was given '{dec}'.")

		# the same bound radial_query places on its radius: beyond 180 degrees a
		# semi-major axis has no meaning on a sphere
		if not (0 <= semi_major_axis <= 180):
			raise ValueError(f"The value for semi_major_axis must be in the range [0,180]; was given '{semi_major_axis}'.")
		if not (0 < axis_ratio <= 1):
			raise ValueError(f"The value for axis_ratio must be in the range (0,1]; was given '{axis_ratio}'.")

		center_ra = ra
		center_dec = dec

		match_ra = list()
		match_dec = list()
		match_key = list()
		match_distance = list()

		# The candidate ranges are always computed at depth 30.
		fulls, partials = self._find_candidate_ellipse_ranges(center_ra, center_dec, semi_major_axis,
															  axis_ratio, position_angle)

		merged_ranges = self._index_depth_ranges(fulls, partials)

		columns = "ra,dec,key" if return_key else "ra,dec"
		query = self._candidate_row_query(merged_ranges, columns,
										  f"ellipse query at (ra={center_ra}, dec={center_dec}) with "
										  f"semi_major_axis={semi_major_axis} deg, axis_ratio={axis_ratio}, "
										  f"position_angle={position_angle} deg")

		center_canonical = _canonical_position(center_ra, center_dec)

		for row in self._select_candidate_rows(query):
			row_ra, row_dec = row["ra"], row["dec"]
			# As in radial_query, a point that is an equivalent representation of the
			# center (ra=360 vs ra=0, any ra at a pole) is inside any ellipse of
			# non-negative size, whatever floating point residual the test computes.
			if _canonical_position(row_ra, row_dec) == center_canonical:
				is_match = True
			else:
				# the C function is called directly rather than through
				# QLSC.point_in_ellipse: the parameters have already been validated
				# above, and revalidating them for every candidate row is the whole
				# cost of this loop. The test itself is identical.
				is_match = q3c.in_ellipse(ra=row_ra, dec=row_dec,
										  ra_cen=center_ra, dec_cen=center_dec,
										  semi_major_axis=semi_major_axis,
										  axis_ratio=axis_ratio,
										  position_angle=position_angle)
			if is_match:
				match_ra.append(row_ra)
				match_dec.append(row_dec)
				if return_key:
					match_key.append(row["key"])
				if return_distance:
					if _canonical_position(row_ra, row_dec) == center_canonical:
						match_distance.append(0.0)
					else:
						match_distance.append(_sindist2deg(sindist(row_ra, row_dec, center_ra, center_dec)))

		return _match_result(match_ra, match_dec, match_key, match_distance, return_key, return_distance)

