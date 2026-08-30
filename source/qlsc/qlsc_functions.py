
#from . import q3c
from . import q3c
from .utilities import _as_float_scalar

def distance(ra1, dec1, ra2, dec2) -> float:
	'''
	Calculates the angular distance between two points on a sphere, in degrees.

	:param ra1: right ascension, first coordinate (degrees)
	:param dec1: declination, first coordinate (degrees)
	:param ra2: right ascension, second coordinate (degrees)
	:param dec2: declination, second coordinate (degrees)
	:returns: the angular separation of the two points in degrees, in [0,180]
	'''
	return q3c.distance(_as_float_scalar(ra1, "ra1"), _as_float_scalar(dec1, "dec1"),
	                    _as_float_scalar(ra2, "ra2"), _as_float_scalar(dec2, "dec2"))

def sindist(ra1, dec1, ra2, dec2) -> float:
	'''
	Calculates the haversine of the angular distance between two points on a sphere.

	Despite the name, the value is **not** ``sin(d)``: it is ``sin²(d/2)``, the
	haversine of the separation. Q3C computes it because it is monotonic in ``d`` over
	the whole range [0,180] degrees and costs three trigonometric evaluations, so a
	cone search can compare it against ``sin²(radius/2)`` without ever taking an
	arcsine. Unlike ``sin(d)`` it increases all the way to 180 degrees rather than
	turning over at 90, which is what makes it usable as a distance ordering.

	To recover the separation in degrees: ``2 * degrees(asin(sqrt(sindist(...))))``,
	or just call :py:func:`distance`.

	:param ra1: right ascension, first coordinate (degrees)
	:param dec1: declination, first coordinate (degrees)
	:param ra2: right ascension, second coordinate (degrees)
	:param dec2: declination, second coordinate (degrees)
	:returns: sin²(d/2) for the angular separation d, a value in [0,1]
	'''
	return q3c.sindist(_as_float_scalar(ra1, "ra1"), _as_float_scalar(dec1, "dec1"),
	                   _as_float_scalar(ra2, "ra2"), _as_float_scalar(dec2, "dec2"))

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
	if facenum not in (0, 1, 2, 3, 4, 5):
		raise ValueError(f"The face number must be an integer in [0,5]; was given '{facenum}'.")
	return q3c.xy2facenum(_as_float_scalar(x, "x"), _as_float_scalar(y, "y"), facenum)
