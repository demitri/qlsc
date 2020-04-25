
#from . import q3c
from . import q3c

def distance(ra1, dec1, ra2, dec2) -> float:
	'''
	Calculates the angular distance between two points on a sphere.
	
	:param ra1: right ascension, first coordinate (degrees)
	:param dec1: declination, first coordinate (degrees)
	:param ra2: right ascension, second coordinate (degrees)
	:param dec2: declination, second coordinate (degrees)
	'''
	return qlsc.q3c.distance(ra1, dec1, ra2, dec2)

def sindist(ra1, dec1, ra2, dec2) -> float:
	'''
	Calculates the sine of the angular distance between two points on a sphere.
	
	:param ra1: right ascension, first coordinate (degrees)
	:param dec1: declination, first coordinate (degrees)
	:param ra2: right ascension, second coordinate (degrees)
	:param dec2: declination, second coordinate (degrees)
	'''
	return qlsc.q3c.sindist(ra1, dec1, ra2, dec2)

def xy2facenum(x:float, y:float, facenum:int) -> int:
	'''
	Convert an x,y coordinate pair on the given face number to the corresponding cube face number.
	
	If -1 <= x <= 1 and -1 <= y <= 1 (i.e. within the coordinate system of the face),
	the face number returned will be the same as provided. Otherwise, this function
	is useful for determining the face number outside the reference of the given face.
	
	:param x: x coordinate locate to the provided face number
	:param y: y coordinate locate to the provided face number
	:param facenum: the face number that is the reference for the provided (x,y)
	:returns: the face number for the given (x,y) coordinates
	'''
	return qlsc.q3c.xy2facenum(x, y, facenum)
