#!/usr/local/env python

'''
This script demonstrates how the QLSCIndex class can be used.

In this example, coordinates are stored in an in-memory SQLite database
which is not persisted to disk; it will be deleted when the script exits.
'''

from qlsc import QLSC, QLSCIndex
import numpy as np

input_file = "1e5_wise.txt"

# create a sphere with 6,144 bins (6 * 2**(2*bin_level))
qlsc = QLSC(bin_level=5)

# create an index:
# by default, an in memory SQLite database is used
qlsc_index = QLSCIndex(qlsc=qlsc)

ra_list = list()
dec_list = list()
#i = 0

with open(input_file) as f:
	f.readline() # skip header
	for line in f:
		cntr, ra, dec = line.rstrip("\n").split("\t")
		ra_list.append(float(ra))
		dec_list.append(float(dec))
#		i += 1
#		if i > 4e5:
#			break

data = np.dstack((ra_list,dec_list))
data = np.squeeze(data)
qlsc_index.add_points(ra=data[:,0], dec=data[:,1])

print("data read")

print(qlsc_index.radial_query(36.1052873, -89.4459347, radius=500/3600)) # radius in arcsec


