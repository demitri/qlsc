#!/usr/bin/env python

import os.path
import pyq3c
import sqlite3

input_file = "1m_wise.txt"

q3c_index = pyq3c.Q3CIndex()

db_filename = "q3c_db.sqlite3"
q3c_index.data_source = db_filename

if q3c_index.number_of_points == 0:
	i = 0

	points = list()
	with open(input_file) as f:
		f.readline()
		for line in f:
			cntr, ra, dec = [float(x) for x in line.rstrip("\n").split("\t")]
			#q3c_index.add_point(ra, dec)
			points.append((ra,dec))
			i += 1
			if i > 1e5:
				break
				#q3c_index.add_points(points=points)
				#points = list()
				#i = 0

		print("data finished reading")
	q3c_index.add_points(points=points)


coordinates = q3c_index.radial_query(144.11901, -77.43726, radius=500/3600)

for c in coordinates:
	print(c)
