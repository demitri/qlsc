
import os
import pathlib
import sqlite3
import contextlib
from typing import Iterable, Tuple
from urllib.parse import urlparse
from math import pow, sin
from math import radians as deg2rad

import numpy as np

try:
	import astropy
	astropy_available = True
except ImportError:
	astropy_available = False

try:
	import pandas
	pandas_available = True
except ImportError:
	pandas_available = False

from . import _q3c_wrapper # as q3c
from ._q3c_wrapper import radial_query_it, sindist

class Q3C:

	def __init__(self, nside:int=1073741824):
		'''
		
		The default value of nside (1073741824) is the one used by the Q3C PostgreSQL plugin.
		
		@param nside number of quadtree subdivisions
		'''		
		#self.sqlite_db = None
		#self.main_dbcursor = None

		# store hprm struct values - hprm is a pointer (struct) to main Q3C structure
		# used by most Q3C calls
		self._hprm = _q3c_wrapper.init_q3c(nside=nside)
		if self._hprm is None:
			raise Exception("Internal Q3C data structure could not be created.")
			
	def ang2ipix(self, ra:float, dec:float) -> int:
		'''
		Convert an ra,dec position to an ipix number.
		
		:param ra: right ascension (degrees)
		:param dec: declination (degrees)
		'''
		return _q3c_wrapper.ang2ipix(self._hprm, ra, dec)
	
	def ipix2ang(self, ipix:int) -> Tuple[float, float]:
		'''
		Convert an ipix number to ra,dec in degrees.
		
		:param ipix: ipix number
		:returns: ra,dec in degrees as a tuple
		'''
		return _q3c_wrapper.ipix2ang(self._hprm, ipix)
				
	@property
	def nsides(self) -> int:
		'''
		'''
		return _q3c_wrapper.nsides(self._hprm)
	
	def face_number(self, ra:float, dec:float) -> int:
		'''
		Return the cube face number for the provided coordinate; in [0-5].
		
		@param ra right ascension in degrees
		@param dec declination in degrees
		'''
		return _q3c_wrapper.facenum(self._hprm, ra, dec)
	
	def pixarea(self, ipix:int, depth:int) -> float: # rename to "ipix_area"?
		'''
		Return the area of a given Q3C pixel in steradians for a given ipix and depth.
		
		depth = 1  -> smallest pixel (1,729,382,256,910,270,000 pixels on sphere)
		depth = 30 -> pixel is one whole face (6 on sphere)
		
		No. pixels on sphere = 6 * pixels per face
		pixels per face = 2 ^ (2 * (30-depth))

		pixel area = 4pi / no. pixels on sphere
		'''
		if not (1 <= depth <= 30):
			raise ValueError(f"The depth provided {depth} is out of the range 1-30")
		
		return _q3c_wrapper.pixarea(self._hprm, ipix, depth)
	
	def ipixcenter(self, ra:float, dec:float, depth:int) -> int:
		'''
		Returns the ipix value of the pixel center at certain pixel depth covering the specified (ra,dec). [??]
	
		:param ra: right ascension (degrees)
		:param dec: declination (degrees)
		:param depth:
		:returns: 
		'''
		if not (1 <= depth <= 30):
			raise ValueError(f"The depth provided {depth} is out of the range 1-30")
		return ( _q3c_wrapper.ang2ipix(self._hprm, ra, dec) >> (2*depth) << (2*depth) ) + (1 << (2 * (depth-1))) - 1

class Q3CIndex:
	
	def __init__(self, q3c:Q3C=Q3C()):
		'''
		'''
		self.q3c = q3c

		self.ra_key = "ra"
		self.dec_key = "dec"
		self.database_tablename = "q3c_table"

		self._data_source = sqlite3.connect(":memory:")
		self._init_sqlite_db()
		
	@property
	def data_source(self):
		'''
		Data structure containing the ra,dec coordinate points.
		'''
		return self._data_source
		
	@data_source.setter
	def data_source(self, new_value):
		'''
		Setter for the 'data_source' property.
		Accepted values:
			* Default value: create an in-memory SQLite database.
			* A path-like object pointing to an SQLite file:
				* If not found, an SQLite file will be created.
				* If a file is present, it will be opened.
			* An SQLite URI, which allows additional read options (e.g. read-only mode).
			* Numpy structured array. This will require the "ra_key" and "dec_key" properties to be set (if not the defaults of "ra", "dec").
			* TODO: astropy.table.Table (https://astropy.readthedocs.io/en/stable/api/astropy.table.Table.html#astropy.table.Table)
			* TODO: pandas.DataFrame (https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html#pandas.DataFrame)
		'''
		# do we have an existing database open?
		if isinstance(new_value, sqlite3.Connection):
			self._data_source.close()
		
		# if a string or path-like object, interpret as an SQLite database
		if isinstance(new_value, (str, pathlib.Path)):

			new_database = True
			is_uri = new_value.startswith("file:")
			if is_uri:
				scheme, netloc, path, params, query, fragment = urlparse(new_value)

			if new_value == ":memory:":
				self._data_source = sqlite3.connect(new_value)
				
			else:
				if is_uri:
					filepath = path
				else:
					filepath = new_value
				new_database = not os.path.exists(filepath)
					
				try:
					self._data_source = sqlite3.connect(new_value, uri=True)
				except sqlite3.OperationalError as e:
					if new_database:
						raise Exception("Unable to create database at specified path.")
					else:
						raise Exception(f"Found file at path 'new_value', but unable to open.")
						
				if new_database:
					self._init_sqlite_db()
						
		elif isinstance(new_value, numpy.ndarray):
			self._data_source = new_value
		elif astropy_available and isinstance(new_value, astropy.table.Table):
			self._data_source = new_value
			raise NotImplementedError("not tested")
		elif pandas_available and isinstance(new_value, pandas.DataFrame):
			raise NotImplementedError("not tested")
		elif new_value is None:	
			raise ValueError("The data source cannot be set to 'None'.")
	
	def __del__(self):
		''' Destructor method. '''
		# if we had an open SQLite connection, close it
		if isinstance(self._data_source, sqlite3.Connection):
			self._data_source.close()
	
	def _init_sqlite_db(self):
		'''
		'''
		#self.sqlite_db = sqlite3.connect(':memory:')
		#con.create_function("ang2ipix", 1, self.ang2ipix)
		assert isinstance(self.data_source, sqlite3.Connection), "_init_sqlite_db called on an object that isn't an SQLite database connection!"
		
		#self.data_source.create_function("q3c_ang2ipix", 3, self.q3c.ang2ipix)

# 		def plus(a,b):
# 			return a+b
		
		#self.data_source.create_function("q3c_ang2ipix", 2, plus)

		with contextlib.closing(self.data_source.cursor()) as cursor:
			cursor.execute(f"CREATE TABLE {self.database_tablename} ({self.ra_key} REAL, {self.dec_key} REAL)") # ipix
			#cursor.execute(f"CREATE INDEX ipix_idx ON {self.database_tablename} (ipix)")

# 			cursor.execute("SELECT q3c_ang2ipix(12, 34)")
# 			ipix = cursor.fetchone()[0]
# 			print(ipix)

	def add_point(self, ra:float=None, dec:float=None):
		'''
		
		'''
		try:
			with self.data_source:
				with contextlib.closing(self.data_source.cursor()) as cursor:
					cursor.execute(f"INSERT INTO {self.database_tablename} ({self.ra_key}, {self.dec_key}) VALUES (?, ?)", (ra, dec)) #, self.ang2ipix(ra,dec))
		except sqlite3.IntegrityError:
			raise NotImplementedError()

	def add_points(self, ra:Iterable=None, dec:Iterable=None, points:Iterable[Tuple]=None):
		'''
		
		:param ra: an Iterable (e.g. list) of right ascension points in degrees
		:param dec: an Iterable (e.g. list) of declination points in degrees
		:param points: an Iterable of tuples of (ra,dec) pairs
		'''
		if not any([ra,dec,points]):
			raise ValueError("'ra','dec' OR 'points' must be specified.")
		if all([ra, dec, points]):
			raise ValueError("Only 'ra','dec' OR 'points' can be specified.")
		elif points is None and not all([ra,dec]):
			raise ValueError("If 'ra' or dec' is given, then other parameter must also be provided.")
		
		try:
			count = 0
			with self.data_source:
				with contextlib.closing(self.data_source.cursor()) as cursor:
					query = f"INSERT INTO {self.database_tablename} ({self.ra_key}, {self.dec_key}) VALUES (?, ?)"
					if ra:
						for i in range(len(ra)):
							cursor.execute(query, (ra[i], dec[i]))
							count += 1
							if count > 50000:
								self.data_source.commit()
								count = 0
					else:
						for ra,dec in points:
							cursor.execute(query, (ra, dec))
							count += 1
							if count > 50000:
								self.data_source.commit()
								count = 0
		except sqlite3.IntegrityError:
			raise NotImplementedError()
	
	@property
	def number_of_points(self) -> int:
		''' Return the number of points in the index. '''
		if isinstance(self.data_source, sqlite3.Connection):
			with contextlib.closing(self.data_source.cursor()) as cursor:
				cursor.execute(f"SELECT COUNT(*) FROM {self.database_tablename}")
				return cursor.fetchone()[0]
				
		else:
			raise NotImplementedError()
		
	
	def radial_query(self, center_ra:float, center_dec:float, radius:float):
		'''
		
		:param ra: right ascension (degrees)
		:param dec: declination (degrees)
		:param ra_center:  (degrees)
		:param dec_center:  (degrees)
		:param radius: radius (degrees)
		'''
				
		matches = list()

		if isinstance(self.data_source, sqlite3.Connection):
			with contextlib.closing(self._data_source.cursor()) as cursor:
				for ra, dec in cursor.execute(f"SELECT {self.ra_key}, {self.dec_key} FROM q3c_table"):
					if self._radial_match(ra, dec, center_ra, center_dec, radius):
					   matches.append((ra,dec))
		
		elif isinstance(self.data_source, numpy.ndarray):
			data = self.data_source
			for idx in range(len(data)):
				ra = data[ra_key][i]
				dec = data[dec_key][i]
				if self._radial_match(ra, dec, center_ra, center_dec, radius):
					matches.append((ra,dec))
			
		elif pandas_available and isinstance(self.data_source, pandas.DataFrame):
			data = self.data_source
			for idx in range(len(data)):
				ra = data[ra_key][i]
				dec = data[dec_key][i]
				if self._radial_match(ra, dec, center_ra, center_dec, radius):
					matches.append((ra,dec))
		
		elif astropy_available and isinstance(self.data_source, astropy.table.Table):
			data = self.data_source
			for idx in range(len(data)):
				ra = data[ra_key][i]
				dec = data[dec_key][i]
				if self._radial_match(ipix, ra, dec, center_ra, center_dec, radius):
					matches.append((ra,dec))
		
		return matches

	def _radial_match(self, ra, dec, center_ra, center_dec, radius):
		'''
		Returns 'True' if the provided ipix value is within radius of (center_ra,center_dec) (all in degrees).
		'''
		# indexing makes this step fast - this is the slow part
		# index is ra,dec in, ipix out
		ipix = self.q3c.ang2ipix(ra, dec)
		hprm = self.q3c._hprm
		return ((ipix>=radial_query_it(hprm,center_ra,center_dec,radius,0,1)  and ipix<radial_query_it(hprm,center_ra,center_dec,radius,1,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,2,1)  and ipix<radial_query_it(hprm,center_ra,center_dec,radius,3,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,4,1)  and ipix<radial_query_it(hprm,center_ra,center_dec,radius,5,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,6,1)  and ipix<radial_query_it(hprm,center_ra,center_dec,radius,7,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,8,1)  and ipix<radial_query_it(hprm,center_ra,center_dec,radius,9,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,10,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,11,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,12,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,13,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,14,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,15,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,16,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,17,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,18,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,19,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,20,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,21,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,22,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,23,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,24,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,25,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,26,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,27,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,28,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,29,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,30,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,31,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,32,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,33,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,34,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,35,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,36,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,37,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,38,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,39,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,40,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,41,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,42,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,43,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,44,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,45,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,46,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,47,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,48,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,49,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,50,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,51,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,52,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,53,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,54,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,55,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,56,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,57,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,58,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,59,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,60,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,61,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,62,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,63,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,64,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,65,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,66,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,67,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,68,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,69,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,70,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,71,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,72,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,73,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,74,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,75,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,76,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,77,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,78,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,79,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,80,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,81,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,82,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,83,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,84,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,85,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,86,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,87,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,88,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,89,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,90,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,91,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,92,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,93,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,94,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,95,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,96,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,97,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,98,1) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,99,1)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,0,0)  and ipix<radial_query_it(hprm,center_ra,center_dec,radius,1,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,2,0)  and ipix<radial_query_it(hprm,center_ra,center_dec,radius,3,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,4,0)  and ipix<radial_query_it(hprm,center_ra,center_dec,radius,5,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,6,0)  and ipix<radial_query_it(hprm,center_ra,center_dec,radius,7,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,8,0)  and ipix<radial_query_it(hprm,center_ra,center_dec,radius,9,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,10,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,11,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,12,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,13,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,14,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,15,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,16,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,17,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,18,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,19,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,20,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,21,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,22,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,23,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,24,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,25,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,26,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,27,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,28,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,29,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,30,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,31,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,32,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,33,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,34,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,35,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,36,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,37,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,38,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,39,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,40,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,41,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,42,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,43,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,44,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,45,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,46,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,47,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,48,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,49,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,50,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,51,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,52,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,53,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,54,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,55,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,56,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,57,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,58,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,59,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,60,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,61,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,62,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,63,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,64,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,65,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,66,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,67,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,68,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,69,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,70,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,71,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,72,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,73,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,74,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,75,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,76,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,77,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,78,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,79,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,80,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,81,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,82,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,83,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,84,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,85,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,86,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,87,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,88,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,89,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,90,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,91,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,92,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,93,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,94,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,95,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,96,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,97,0)) or \
			    (ipix>=radial_query_it(hprm,center_ra,center_dec,radius,98,0) and ipix<radial_query_it(hprm,center_ra,center_dec,radius,99,0))  \
			   ) and sindist(ra,dec,center_ra,center_dec) < pow(sin(deg2rad(radius)/2.), 2)
			
