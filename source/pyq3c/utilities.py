
'''
Note:

This file contains utilities used for calculations in this package.
They are not intended to be general purpose and often specifically cannot
be used as such.
'''

import math
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
					 cos_dec])        # z

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

def _normalize_ang(points, copy=False): #, radians=False):
	'''
	Normalize ra,dec coordinates on a sphere to 0 <= α < 360 and -π <= δ <= π.
	
	Normalizing the right ascension is straightforward and is independent of declination.
	However, normalizing the declination (e.g. δ = -140) may result in a change
	in right ascension by 180 degrees.
	
	:param points: an array of α,δ points, array shape (n,2)
	:param copy: work on and return a copy of the array without modifying the original
	:param radians: True is the coordinate units are in radians, degrees otherwise NOT YET SUPPORTED
	'''
	
	# .. todo: add support for radians
	
	if copy:
		points = np.copy(points)
	
#	if radians:
#		pi = math.pi
#	else:
	pi = 180.
	
	ra  = points[:,0]
	dec = points[:,1]
	
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
	
	return points

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


