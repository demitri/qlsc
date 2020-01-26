#!/usr/local/env python

import pyq3c

input_file = "/Users/demitri/1m_wise.txt"

q3c = pyq3c.Q3C()

ra_list = list()
dec_list = list()

i = 0

with open(input_file) as f:
	f.readline()
	for line in f:
		cntr, ra, dec = line.rstrip("\n").split("\t")
		ra_list.append(float(ra))
		dec_list.append(float(dec))
		i += 1
		if i > 4e5:
			break

data = np.array(zip(ra_list, dec_list), dtype=[('ra', np.double), ('dec', np.double)]

print("data loaded")

x = data.where(coordinate_in_circle(36.1052873, -89.4459347, radius=500/3600)
print(x)