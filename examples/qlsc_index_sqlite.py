#!/usr/bin/env python

'''
This script demonstrates how the QLSCIndex class can be used.

In this example, coordinates are stored in an SQLite database
that persists on disk, useful for repeated runs.
'''

import os.path
import sqlite3
import numpy as np
from qlsc import QLSC, QLSCIndex

# File contains 10000 points from the WISE catalog.
#input_file = "1e5_wise.txt"
input_file = "../_ignore/1m_wise.txt"

# define a QLSC scheme at a specific bin level
qlsc = QLSC()

qlsc_index = QLSCIndex(qlsc=qlsc)

# The QLSCIndex object will create an SQLite database
# to hold coordinates given to it. Define the database name here (optional).
db_filename = "qlsc_db.qlscidx"
qlsc_index.data_source = db_filename

# The database is not deleted, so repeated runs of this
# script will not reload the points from the file above.

if qlsc_index.number_of_points == 0:
	i = 0

	points = list()
	with open(input_file) as f:
		f.readline()
		for line in f:
			cntr, ra, dec = [float(x) for x in line.rstrip("\n").split("\t")]
			#qlsc_index.add_point(ra, dec)
			points.append((ra,dec))
			i += 1
			if i > 3e5:
				break
				#qlsc_index.add_points(points=points)
				#points = list()
				#i = 0

		print("data finished reading")
	points = np.array(points)
	print(points[0:10])
	qlsc_index.add_points(points=points)

ra = 320.0752
dec=  -86.1787
print(f"ipix = {qlsc.ang2ipix(ra,dec)}")
coordinates = qlsc_index.radial_query(ra, dec, radius=500/3600)

for c in coordinates:
	print(c)
