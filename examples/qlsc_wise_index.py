#!/usr/bin/env python

import os.path
from qlsc import QLSCIndex
import sqlite3

input_file = "1m_wise.txt"

qlsc_index = QLSCIndex()

db_filename = "qlsc_db.sqlite3"
qlsc_index.data_source = db_filename

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
			if i > 1e5:
				break
				#qlsc_index.add_points(points=points)
				#points = list()
				#i = 0

		print("data finished reading")
	qlsc_index.add_points(points=points)


coordinates = qlsc_index.radial_query(144.11901, -77.43726, radius=500/3600)

for c in coordinates:
	print(c)
