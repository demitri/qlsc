
'''
Note:

This file contains utilities used for calculations in this package.
They are not intended to be general purpose and often specifically cannot
be used as such.
'''

import math
import collections
from typing import Union

import numpy as np
from numpy import deg2rad, rad2deg, cos, sin, arctan2, sqrt, fmod, trunc, sign

def _ang2cartesian(ra=None,dec=None,points=None):
	'''
	Convert ra,dec points in degrees to Cartesian coordinates.
	Calculation is assumed to be taken on a unit sphere (r=1).
	Note that 0 <= α < 360 and -π <= dec <= +π.

	:param ra: right ascension (degrees); either as a value or Numpy array
	:param dec: right ascension (degrees); either as a value or Numpy array
	:param points: a Numpy array of ra,dec points, shape (n,3)
	'''
	
	# Account for declination measured from equatorial plane (elevation angle)
	# vs degrees down from pole (inclination). The calculations here assume
	# an inclination angle.
	
	if any([ra,dec]) and points is not None:
		raise ValueError("Only 'ra','dec' OR 'points' can be specified; not both.")
	
	if ra is None:
		ra = points[:,0]
		dec = points[:,1]
	ra = deg2rad(ra)
	dec = deg2rad(90-dec)
	cos_dec = np.cos(dec)
	sin_dec = np.sin(dec)
	cos_ra = np.cos(ra)
	sin_ra = np.sin(ra)
	return np.array([sin_dec * cos_ra, # x
					 sin_dec * sin_ra, # y
					 cos_dec],         # z
					 dtype=np.double)

def _cartesian2ang(v):
	'''
	Convert Cartesian points into spherical coordinates.
	
	:param v: vector of Cartesian points
	'''
	a = np.zeros((v.shape[0], 2), dtype=np.double)
	x,y,z = v[:,0], v[:,1], v[:,2]
	a[:,0] = arctan2(y, x)
	a[:,1] = math.pi/2 - arctan2(sqrt(x**2+y**2), z)
	return a

def _elevation2inclination(angles, degrees=False):
	'''
	Convert elevation angles (measured from the reference plane) to inclination angles (measured from zenith).
	
	Converts angles that range from [-π,π] to [0,2π].
	'''
	# Ref: https://en.wikipedia.org/wiki/Spherical_coordinate_system#Cartesian_coordinates
	# Convert to Cartesian using one system then back to spherical using the other; assume a unit sphere (r=1).
	
	# this method is here for reference as much as anything else
	if degrees:
		return np.rad2deg(np.arccos(np.sin(np.deg2rad(angles))))
	else:
		return np.arccos(np.sin(angles))

def _inclination2elevation(angles, degrees=False):
	'''
	Convert elevation angles (measured from the zenith) to elevation angles (measured from the reference plane).
	
	Converts angles that range from [0,2π] to [-π,π].
	'''
	# this method is here for reference as much as anything else
	if degrees:
		return np.rad2deg(np.arcsin(np.cos(np.deg2rad(angles))))
	else:
		return np.arcsin(np.cos(angles))

def _rotation_matrix_x(angle) -> np.matrix:
	'''
	Returns a rotation matrix for rotating around the y axis.

	:param angle: angle of rotation in degrees
	'''
	angle = deg2rad(angle)
	return np.matrix([
		[    1,        0,            0        ],
		[    0,  np.cos(angle), np.sin(angle) ],
		[    0, -np.sin(angle), np.cos(angle) ]
		])

def _rotation_matrix_y(angle) -> np.matrix:
	'''
	Returns a rotation matrix for rotating around the y axis.

	:param angle: angle of rotation in degrees
	'''
	angle = deg2rad(angle)
	return np.matrix([
		[  np.cos(angle),  0,  np.sin(angle) ],
		[       0,         1,        0       ],
		[ -np.sin(angle),  0,  np.cos(angle) ]
		])

def _rotation_matrix_z(angle) -> np.matrix:
	'''
	Returns a rotation matrix for rotating around the y axis.

	:param angle: angle of rotation in degrees
	'''
	angle = deg2rad(angle)
	return np.matrix([
		[  np.cos(angle),  0, np.sin(angle) ],
		[ -np.sin(angle),  0, np.cos(angle) ],
		[       0,         0,      1        ]
		])

def _normalize_ra(ra):
	'''
	Normalize right ascension values in degrees to [0,360].
	'''
	is_list = isinstance(ra, list)
	if is_list:
		ra = np.array(ra, dtype=np.double)
		
	if isinstance(ra, np.ndarray):
		if ra.dtype not in [np.float, np.double]:
			ra = ra.astype(np.double)
	
	if isinstance(ra, collections.abc.Iterable):
# 		idx = ra < 0
# 		ra[idx] = np.fmod(ra,360) + 360
# 		id ra >= 360
# 		ra[idx] = fmod(ra,360)

		idx = (ra < 0.) | (ra >=360.)
		ra[idx] -= np.trunc(ra[idx]/360.)*360.
		idx = ra < 0.
		ra[idx] += 360.

	else:
		if ra < 0.:
			ra = fmod(ra,360.) + 360.
		elif ra >= 360.:
			ra = fmod(ra,360.)
	if is_list:
		return list(ra)
	else:
		return ra

def _normalize_ang(ra:Union[collections.abc.Iterable,float]=None,
				   dec:Union[collections.abc.Iterable,float]=None,
				   points:collections.abc.Iterable=None,
				   copy:bool=False): #, radians=False):
	'''
	Normalize ra,dec coordinates on a sphere to 0 <= α < 360 and -π <= δ <= π.
	
	Normalizing the right ascension is straightforward and is independent of declination.
	However, normalizing the declination (e.g. δ = -140) may result in a change
	in right ascension by 180 degrees.
	
	There are two usage patterns; one is to modify the points in the arrays provided:
	
	.. code-block:: python
	
	     _normalize_ang(points) # no return value, points modified in place
	     
	and the other is to return new arrays with the normalized coordinates:
	
	.. code-block:: python
	
	     new_points = _normalize_ang(points) # 'points' is left unmodified

	Note that this function only works on float data types (specifically, ``np.float32``, ``np.float64``).
	Passing an array of integer values will raise an exception since this method does not
	want to modify given arrays unless specifically requested. Thus,
	
	.. code-block:: python
	
	     _normalize_ang(int_points) # raises exception
	     new_points = _normalize_ang(int_points, copy=True) # works, new array dtype = np.double

	The returned form will be a single array with shape (n,2). Also note that the number of values
	returned will match the number of values passed in:

	.. code-block:: python

		ra, dec = _normalize_ang(ra=ra, dec=dec)
		points  = _normalize_ang(points=points)
	
	:param ra: an array of right ascension points (degrees)
	:param dec: an array of declination points (degrees)
	:param points: an array of α,δ points, array shape (n,2)
	:param copy: work on and return a copy of the array without modifying the original
	'''
#	:param radians: True is the coordinate units are in radians, degrees otherwise NOT YET SUPPORTED
	
	# .. todo: add support for radians

	if (ra is not None or dec is not None) and points is not None:
		raise ValueError("Cannot specify 'points' with 'ra' or 'dec'.")
	if True in [x is None for x in [ra,dec]] and points is None:
		raise ValueError("Only one of 'ra' or 'dec' was set; they should be set together.")
	
	return_scalar = False
	if ra is not None:
		if np.isscalar(ra):
			ra = np.array([ra], dtype=np.double)
			dec = np.array([dec], dtype=np.double)
			return_scalar = True
	
	if copy is False:
		# check that the data types are floating point
		if ra is not None:
			if not all([x.dtype in [np.float32,np.float64] for x in [ra,dec]]):
				ValueError("The 'ra','dec' values are not floating point (i.e. np.float32 or np.float64); set 'copy=True' to get around this.")
		else:
			if points is not None and points.dtype not in [np.float32,np.float64]:
				ValueError("The 'points' values are not floating point (i.e. np.float32 or np.float64); set 'copy=True' to get around this.")
	
	# are we returning (ra,dec) or (points)?
	return_ra_dec = None
	
	if points is not None:
		points = np.atleast_2d(points)
		if copy:
			if points.dtype not in [np.float32,np.float64]:
				points = np.copy(points).astype(np.float64)
			else:
				points = np.copy(points)
		ra  = points[:,0]
		dec = points[:,1]
		return_ra_dec = False
	else:
		# points is None, ra.dec defined
		if len(ra) != len(dec):
			raise ValueError(f"The lengths of 'ra' and 'dec' are not equal ({len(ra)}, {len(dec)}).")
		if copy:
			if ra.dtype not in [np.float32,np.float64]:
				ra  = np.copy(ra).astype(np.float64)
			else:
				ra  = np.copy(ra)
			if dec.dtype not in [np.float32,np.float64]:
				dec = np.copy(dec).astype(np.float64)
			else:
				dec = np.copy(dec)
		return_ra_dec = True
			
#	if radians:
#		pi = math.pi
#	else:
	pi = 180.
	
	# normalize dec
	# -------------
	# normalize angles ≥ 360°
	idx = abs(np.trunc(dec/(2*pi))) > 0
	dec[idx] =  fmod(dec[idx], (2*pi))
	
	# normalize angles -360° < δ ≤ -270°
	idx = np.trunc(dec/(1.5*pi)) == -1
	dec[idx] = pi/2 + fmod(dec[idx], (1.5*pi))

	# normalize angles 270° ≤ δ < 360°
	idx = np.trunc(dec/(1.5*pi)) == 1
	dec[idx] = fmod(dec[idx], (1.5*pi)) - pi/2
	
	# normalize angles >= 180°
	idx = abs(np.trunc(dec/pi)) == 1
	dec[idx] = -fmod(dec[idx], pi)
	ra[idx] =  ra[idx] + pi # we've changed quadrants, adjust ra
	
	# normalize angles between -90° and 0°
	idx = np.trunc(dec/(pi/2)) == -1
	dec[idx] = -fmod(dec[idx], pi/2) - pi/2
	ra[idx] =  ra[idx] + pi
	
	# normalize angles between 0° and 90°
	idx = np.trunc(dec/(pi/2)) == 1
	dec[idx] =  pi/2 - fmod(dec[idx], pi/2)
	ra[idx] =  ra[idx] + pi
	
	# normalize ra
	# ------------
	idx = (ra < 0.) | (ra >=360.)
	ra[idx] -= np.trunc(ra[idx]/360.)*360.
	idx = ra < 0.
	ra[idx] += 360.
	
	if return_scalar:
		return ra[0], dec[0]
		
	#return np.squeeze(points)
	if copy:
		if return_ra_dec:
			return ra, dec
		else:
			return np.squeeze(np.dstack((ra, dec)))
	# if not copied, nothing needs to be returned

	# below is test code
	
# 	np.set_printoptions(suppress=True)
# 
# 	n=121
# 	points = np.zeros((n,3), dtype=np.double)
# 	points[:,1] = np.linspace(-600,600,num=n)
# 	points[:,2] = np.linspace(-600,600,num=n)
# 
# 	ra  = points[:,0]
# 	dec = points[:,1]
# 	d4_90 = abs(np.trunc(dec/360)) > 0
# 	dec[d4_90] = fmod(dec[d4_90],360)
# 
# 	d3_90 = np.trunc(dec/270) == -1
# 	dec[d3_90] = 90 + fmod(dec[d3_90],270)
# 
# 	d3_90 = np.trunc(dec/270) == 1
# 	dec[d3_90] = fmod(dec[d3_90],270) - 90
# 
# 	d2_90 = abs(np.trunc(dec/180)) == 1
# 	dec[d2_90] = -fmod(dec[d2_90],180)
# 	ra[d2_90] =  ra[d2_90] + 180
# 
# 	d1_90 = np.trunc(dec/90) == -1
# 	dec[d1_90] = -fmod(dec[d1_90],90) - 90
# 	ra[d1_90] =  ra[d1_90] + 180
# 
# 	d1_90 = np.trunc(dec/90) == 1
# 	dec[d1_90] = 90 - fmod(dec[d1_90],90)
# 	ra[d1_90] =  ra[d1_90] + 180
# 	
# 	print(points)

def great_circle_intersection(circle1_points=None, circle2_points=None,
							  return_spherical=True,
							  closest_to=None):
	'''
	Given two great circles, find the two points where they intersect.
	
	Example usage:
	
	great_circle_intersection(circle1_points=[[0,45],[0,90]], circle2_points=[[0,0],[0,45]])
	
	Ref: https://stackoverflow.com/questions/2954337/great-circle-rhumb-line-intersection

	:param circle1_points: two points on circle 1, degrees
	:param circle2_points: two points on circle 2, degrees
	:param return_spherical: if True, returns points in spherical coordinates (degrees), otherwise Cartesian
	:param closest_to: if a point is provided, only return the intersection point that is closes to this one
	:returns: 
	'''
	
	# array to hold coordinates
	# [[ra1_1, dec1_1], # point 1 circle 1
	#  [ra1_2, dec1_2], # point 2 circle 1
	#  [ra2_1, dec2_1], # point 1 circle 2
	#  [ra2_2, dec2_2], # point 2 circle 2
	
	coords = np.array([circle1_points[0],
					   circle1_points[1],
					   circle2_points[0],
					   circle2_points[1]])

	#print(coords)
	#sys.exit(0)

	coords = deg2rad(coords)
	
	################################################
	#### Intersection of two great circles.
	# Points on great circle 1.
# 	ra1 = deg2rad(-45)
# 	dec1 = deg2rad(45)
# 
# 	ra2 = deg2rad(45)
# 	dec2 = deg2rad(45)
# 
# 	# Points on great circle 2.
# 	cra1 = deg2rad(22.5)
# 	cdec1 = deg2rad(0)
# 
# 	cra2 = deg2rad(22.5)
# 	cdec2 = deg2rad(22.5)

	# 1. Put in Cartesian coordinates

	# [[x,y,z],
	#  [x,y,z],
	#  [x,y,z],
	#  [x,y,z]]
	
	cart = np.zeros((4,3))

	cart[:,0] = np.sin(coords[:,1]) * np.cos(coords[:,0]) # x
	cart[:,1] = np.sin(coords[:,1]) * np.sin(coords[:,0]) # y
	cart[:,2] = np.cos(coords[:,1])                       # z

# 	x1 = sin(dec1) * cos(ra1)
# 	y1 = sin(dec1) * sin(ra1)
# 	z1 = cos(dec1)
# 
# 	x2 = sin(dec2) * cos(ra2)
# 	y2 = sin(dec2) * sin(ra2)
# 	z2 = cos(dec2)
# 
# 
# 	cx1 = sin(cdec1) * cos(cra1)
# 	cy1 = sin(cdec1) * sin(cra1)
# 	cz1 = cos(cdec1)
# 
# 	cx2 = sin(cdec2) * cos(cra2)
# 	cy2 = sin(cdec2) * sin(cra2)
# 	cz2 = cos(cdec2)

	# 2. Get normal to planes containing great circles.
	#    It's the cross product of vector to each point from the origin.

#	N1 = cross([x1, y1, z1], [x2, y2, z2])
#	N2 = cross([cx1, cy1, cz1], [cx2, cy2, cz2])

	N1 = np.cross(cart[0], cart[1])
	N2 = np.cross(cart[2], cart[3])

	# 3. Find line of intersection between two planes.
	#    It is normal to the poles of each plane.

	L = np.cross(N1, N2)

	# 4. Find intersection points.

	X1 = L / np.sqrt(L[0]**2 + L[1]**2 + L[2]**2) 
	X2 = -X1

	if return_spherical:
		ra1  = atan2(X1[1], X1[0])
		#dec1= atan2(np.sqrt(pow(X1[0],2)+pow(X1[1],2)),X1[2]) # for elevation angle defined from z-axis down
		dec1 = atan2(X1[2], np.sqrt(pow(X1[0],2)+pow(X1[1],2))) # for elevation angle defined from x-y plane up

		ra2  = atan2(X2[1], X2[0])
		#dec2= atan2(np.sqrt(pow(X2[0],2)+pow(X2[1],2)),X2[2])
		dec2 = atan2(X2[2], np.sqrt(pow(X2[0],2)+pow(X2[1],2)))
		
		return np.rad2deg(np.array([[ra1,dec1], [ra2,dec2]]))
	else:
		return [X1,X2]

