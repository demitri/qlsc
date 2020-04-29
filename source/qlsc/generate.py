
import math
import numpy as np

def sunflower_points_on_sphere(n:int=1000, radians:bool=False, polar_angle_from_zenith:bool=False):
	'''
	Generate a list of points on a sphere using the golden spiral method.
	
	The default values return points that can be taken as ra,dec in degrees.
	
	Points can be returned where the polar angle θ is measured from the zenith
	where 0 ≤ 0 ≤ π (``polar_angle_from_zenith=True``) or else from the reference
	plane where -π/2 ≤ θ ≤ π/2 (``polar_angle_from_zenith=False``).
	The QLSC index requires the latter.
	
	Points are returned in degrees by default, as required by the QLSC index.
	The same points will always be returned for the same values of `n`.
	The points array is sorted in descending order by the polar angle, from the zenith down.
	
	It can be useful to set ``radians=True`` for plotting.
	
	Ref: https://stackoverflow.com/a/44164075/2712652
	
	:param n: number of points to generate
	:param radians: returns points in radians if True, degrees if False
	:param polar_angle_from_zenith: returns 
	:returns: array of points on sphere as polar angle, azimuthal angle: (φ,θ)
	'''
	
	if not isinstance(n, int):
		# Values like "n=1e6" are valid, but Python creates them as floats.
		# If the value given has no decimal component, convert it to an int.
		if n % 1 == 0:
			n = int(n)
		else:
			raise ValueError(f"The value 'n' must be an integer.")
	
	indices = np.arange(0, n, dtype=np.double) + 0.5

	points = np.zeros((n,2), dtype=np.double) # format: [azimuthal angle, polar angle]

	# polar angle
	#  0 ≤ 0 ≤ π 

	# azimuthal angle
	#  0 ≤ φ ≤ 2π

	# note that theta continually increases (always >0)
	# and is not normalized

	# ---
	# polar angle
	# θ (physics convention), φ (math convention)
	# ---
	points[:,1] = abs(np.arccos(1 - 2*indices/n))

	# -----
	# azimuthal angle
	# φ (physics convention), θ (math convention)
	# -----
	points[:,0] = math.pi * (1 + 5**0.5) * indices

	# normalize azimuthal points on [0,2π]
	idx = np.trunc(points[:,0] / 2*math.pi) > 1
	points[:,0] = np.fmod(points[:,0], 2*math.pi)
	
	if polar_angle_from_zenith is False:
		# convert polar angle to -π/2 ≤ 0 ≤ π/2
		points[:,1] = math.pi/2 - points[:,1]
	
	if radians:
		return points # [theta, phi]
	else:
		return np.rad2deg(points)




