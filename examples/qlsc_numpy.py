#!/usr/local/env python

from qlsc import QLSC, QLSCIndex
import numpy as np

input_file = "wise_selection.txt"

# create a sphere with 6,144 bins (6 * 2**(2*bin_level))
qlsc = QLSC(bin_level=5)

# create an index
# by default, an in memory SQLite database is used
qlsc_index = QLSCIndex(qlsc=qlsc)

ra_list = list()
dec_list = list()
coord_list = list()
i = 0

with open(input_file) as f:
	f.readline()
	for line in f:
		cntr, ra, dec = line.rstrip("\n").split("\t")
		ra_list.append(float(ra))
		dec_list.append(float(dec))
		#coord_list.append([float(ra),float(dec)])
		i += 1
		if i > 4e5:
			break

#ra_list = np.array(ra_list)
#dec_list = np.array(dec_list)
#data = np.array(list(zip(ra_list, dec_list)),
#				dtype=[('ra', np.double), ('dec', np.double)])

data = np.dstack((ra_list,dec_list))
data = np.squeeze(data)
#print(coord_list)
#print(data)
#data = np.array(data, dtype=[('ra', np.double), ('dec', np.double)])
#print(data[:,1])
qlsc_index.add_points(ra=data[:,0], dec=data[:,1])

print("data loaded")

print(qlsc_index.radial_query(36.1052873, -89.4459347, radius=500/3600)) # radius in arcsec
#print(idx)
#print(data[idx])


