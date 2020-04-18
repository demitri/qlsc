
import numpy as np
from numpy import fmod, abs, sign
from pyq3c.utilities import _normalize_ang
np.set_printoptions(suppress=True)

n=121
points = np.zeros((n,3), dtype=np.double)
points[:,1] = np.linspace(-600,600,num=n)
points[:,2] = np.linspace(-600,600,num=n)
#print(points)

old = False

if old:
	ra  = points[:,0]
	dec = points[:,1]
	d4_90 = abs(np.trunc(dec/360)) > 0
	dec[d4_90] = fmod(dec[d4_90],360)

	d3_90 = np.trunc(dec/270) == -1
	dec[d3_90] = 90 + fmod(dec[d3_90],270)

	d3_90 = np.trunc(dec/270) == 1
	dec[d3_90] = fmod(dec[d3_90],270) - 90

	d2_90 = abs(np.trunc(dec/180)) == 1
	dec[d2_90] = -fmod(dec[d2_90],180)
	ra[d2_90] =  ra[d2_90] + 180

	d1_90 = np.trunc(dec/90) == -1
	dec[d1_90] = -fmod(dec[d1_90],90) - 90
	ra[d1_90] =  ra[d1_90] + 180

	d1_90 = np.trunc(dec/90) == 1
	dec[d1_90] = 90 - fmod(dec[d1_90],90)
	ra[d1_90] =  ra[d1_90] + 180
	
	print(points)
else:

	print(_normalize_ang(points))
