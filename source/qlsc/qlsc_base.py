
import os
import math
import pathlib
import sqlite3
import contextlib
from typing import Iterable, Tuple
from urllib.parse import urlparse
from math import pow, sin
from math import radians as deg2rad
import collections

import numpy as np

try:
	import astropy
	astropy_available = True
except ImportError:
	astropy_available = False

try:
	import pandas as pd
	pandas_available = True
except ImportError:
	pandas_available = False

from . import q3c
from .q3c import radial_query_it, sindist

# functions defined in "qlsc_c_module.c" as named in the array "qlsc_methods" are in the qlsc.q3c namespace
#from . import q3c
from .utilities import _normalize_ang

class QLSC:
	'''
	A class that describes a quadrilaterized spherical cube (QLSC).
	ggg
	The quadrilaterized sphere is divided into six faces. The bin level defines the number times a
	face is subdivided. A bin level of 1 divides each face into four "pixels"; bin level two divides
	each of those pixels again into four, etc.

	At the lowest resolution where bin level=0 there is 1 bin per face for a total of 6 bins.
	The next bin level 1 divides each face into four bins, i.e. a total of 24 bins.
	Higher bin levels further divide each bin four more times.
	
	A scheme can be defined by the bin level where:

		.. code-block:: 

			no. bins per face        = (2**bin_level)**2
			number of bins on sphere = no. bins per face * 6
	
	Each bin has a pixel number called an "ipix", which is an integer encoded with the location of the bin.
	
	The default value of `bin_level` = 30 that divides the sphere into 6,442,450,944 bins (i.e. nside=1,073,741,824),
	is the one used by the `Q3C PostgreSQL plugin <https://github.com/segasai/q3c>`_.
	This value is used if neither `nside` or `bin_size` are specified.
	
	:param bin_level: number of times each bin in each face is subdivided by four
	'''		
#	:param nside: number of quadtree subdivisions per face (of six faces)
	def __init__(self, bin_level:int=30): # nside:int=None

#		if all([nside,bin_level]):
#			raise ValueError("Only one of 'nside' or 'bin_level' can be specified in the initialization.")
#		elif nside is None and bin_level is None:
#			nside = 1073741824
		
		# The Q3C code only supports up to bin_level=30
		if not (0 <= bin_level <= 30):
			raise ValueError(f"The value for 'bin_level' must be an integer on [0,30]; was given '{bin_level}'.")
#		if nside and not (0 <= nside <= int(pow(pow(2, 30), 2)):
#			raise ValueError(f"The value for 'nside' must be an integer on [0,30]; was given '{bin_level}'.")
		
		# "nside" as defined by Q3C is the number of bins along one side of the cube face.
		# This means the number of bins on a cube face is nside*nside, and the total
		# number of bins is 6*nside*nside.
		#
		nside = int(pow(2, bin_level))
		
		self.bin_level = bin_level
		
		#self.sqlite_db = None
		#self.main_dbcursor = None

		# store hprm struct values - hprm is a pointer (struct) to main Q3C structure
		# used by most Q3C calls
		self._hprm = q3c.init_q3c(nside=nside)
		if self._hprm is None:
			raise Exception("Internal Q3C data structure could not be created.")
			
	def __repr__(self):
		return f"<{self.__class__.__module__.split('.')[0]}.{self.__class__.__name__} at {hex(id(self))}, bin_level={self.bin_level}>"
	
	def ang2ipix(self, ra:float=None, dec:float=None, points:Iterable=None):
		'''
		Convert an ra,dec position to an ipix number.
		
		Example usage:
		
		.. code-block:: python
		
			ipix = ang2ipix(45,60)
			ipix = ang2ipix(ra=45,dec=60)
			ipix = and2ipix(points=np.array([[45,60],[0,55]])
		
		:param ra: right ascension (degrees); must be in [0,360]
		:param dec: declination (degrees); must be in [-90,90]
		:param points: an array of [ra,dec] points, shape=(n,2)
		:returns: the ipix integer if a single point is provided, an array of ipix integers if an array of points is provided
		'''
		if points is None:
			# check "ra","dec" are both set
			if True in [x is None for x in [ra,dec]]: # can't use any/all with possible "0" values
				raise ValueError(f"Values for 'ra','dec' must be specified. {ra},{dec}")
		else:
			raise NotImplementedError("ang2ipix does yet support a list of points.")
			# check neither "ra","dec" are set
			if ra is not None or dec is not None:
				raise ValueError("Only specify 'ra','dec' OR 'points'.")
		
		print(f"flags: {ra.flags}")
		print(f"(ra,dec) = {type(ra)},{type(dec)}")
		return q3c.ang2ipix(self._hprm, ra, dec)

		if isinstance(ra, collections.abc.Iterable):
			if isinstance(points, np.ndarray):
				points = points.astype(np.double, copy=False) # was getting int64 when using, e.g. (0,45)
			if not (min(dec) >= -90 and max(dec) <= 90):
				raise ValueError("ang2ipix: dec out of range [-90,90]")
			if not (min(ra) >= 0 and max(ra) <= 360):
				raise ValueError("ang2ipix: ra out of range [0,360]")
			# ra is unwrapped in the C code, but should be checked here anyway
		else:
			if not ((-90 <= dec <= 90) and (0 <= ra <= 360)):
				ra, dec = _normalize_ang(np.array([ra,dec]))
#				raise ValueError("ang2ipix: dec out of range [-90,90]")
#			if not (0 <= ra <= 360):
#				raise ValueError("ang2ipix: ra out of range [0,360]")
		
#		if points is None:
		print(f"flags: {ra.flags}")
		return q3c.ang2ipix(self._hprm, ra, dec)
#		else:
# 			# .. todo:: modify C wrapper code to take numpy arrays directly
# 			v = np.vectorize(q3c.ang2ipix)
# 			return v(self._hprm, points[:,0], points[:,1])
# 			#return q3c.ang2ipix(self._hprm, points[:,0], points[:,1])

	def ang2ipix_xy(self, ra:float=None, dec:float=None, points=None) -> dict:
		'''
		Convert an ra,dec position to an ipix number; also returns the face number, and (x,y) location on the square face.
		
		:param ra: right ascension angle (degrees)
		:param dec: declination angle (degrees), must be in [-90,90]
		:param points: an array of ra,dec coordinates (shape: (n,2)) NOT YET SUPPORTED
		:returns: a dictionary with the keys: [`facenum`, `ipix`, `x`, `y`]
		'''
		
		# The underlying C code forces dec > 90 to 90 and dec < -90 to 90.
		# Either normalize the angle here or raise an exception.
		# The underlying C code *does* unwrap RA.
		if not (-90 <= dec <= 90):
			raise ValueError(f"Values of 'dec' must be in the range [-90,90] (was given '{dec}').")
			
		return q3c.ang2ipix_xy(self._hprm, ra, dec)
	

	def ipix2ang(self, ipix:int) -> Tuple[float, float]:
		'''
		Convert an ipix number to the ra,dec coordinate of the "lower left" corner of the pixel, in degrees.
		
		See the method :py:func:`ipix2ang_center` to return the center of the pixel.
		
		This method matches that of the `Q3C PostgreSQL plugin <https://github.com/segasai/q3c>`_.
		
		:param ipix: ipix number
		:returns: ra,dec in degrees as a tuple of the lower left corner of the pixel
		'''
		if isinstance(ipix, int) is False:
			raise ValueError(f"The ipix value must be an integer; was given '{type(ipix)}'.")

		# perform a bounds check
		if not (0 <= ipix < self.nbins):
			raise ValueError(f"The ipix number is out of bounds for this scheme; range: [0,{self.nbins-1}].")

		return q3c.ipix2ang(self._hprm, ipix)

	def ipix2ang_center(self, ipix:int) -> Tuple[float, float]:
		'''
		Convert an ipix number to the ra,dec coordinate of the center of the pixel, in degrees.
		
		See the method :py:func:`ipix2ang` to return the "lower left" corner of the pixel
		which matches what is returned by the `Q3C PostgreSQL plugin <https://github.com/segasai/q3c>`_.
		
		The pixel "center" here is defined in the x,y coordinates;
		the point projected on the sphere is the one returned.
		
		:param ipix: ipix number
		:returns: ra,dec in degrees as a tuple of the center of the pixel
		'''
		if isinstance(ipix, int) is False:
			raise ValueError(f"The ipix value must be an integer; was given '{type(ipix)}'.")

		# perform a bounds check
		if not (0 <= ipix < self.nbins):
			raise ValueError(f"The ipix number is out of bounds for this scheme; range: [0,{self.nbins-1}].")

		# bin_length = 2./self.nside
		hbl = 1. / self.nside # half bin length
		
		facenum, x, y = self.ipix2xy(ipix) # returns lower left corner
		return self.xy2ang(facenum=facenum, x=x+hbl, y=y+hbl)
		
	
	def ipix2xy(self, ipix:int) -> Tuple[int, float, float]:
		'''
		Convert an ipix value to the (x,y) position on the corresponding square face; returns (facenum,x,y).
		
		:param ipix: ipix number
		:returns: (facenum,x,y) as a tuple
		'''
		
		return q3c.ipix2xy(self._hprm, ipix)
	
	def xy2ang(self, facenum:int=None, x:float=None, y:float=None, points:collections.abc.Iterable=None): # -> Tuple[float,float]:
		'''
		Convert an x,y coordinate pair on the given face number to (ra,dec).
		
		:param facenum: the number of the cube face (0=top, 1-4=sides, 5=bottom)
		:param x: x coordinate on cube face
		:param y: y coordinate on cube face
		:param points: an array containing (x,y) coordinate pairs; shape = (n,2)
		:returns: tuple of (ra,dec) corresponding to the location on the cube face
		'''
		# .. todo:: validate ra,dec
		
		if facenum is None:
			raise ValueError("The cube face number ('facenum') must be provided.")
		else:
			if not isinstance(facenum, int):
				# notmally would be ok to cast the value, but this is a safety net to catch 'x' values passed to facenum
				raise ValueError("The value of facenum is expected to be an integer.")
			if not (0 <= facenum <= 5):
				raise ValueError(f"Values of 'facenum' must be an integer between 0 and 5; was given '{facenum}'.")
		if all(v is not None for v in [x, y, points]): # zero is a possible value
			raise ValueError(f"Specify only 'ra' and 'dec', or 'points'")
		if points is None and (x is None or y is None):
			raise ValueError(f"Values for 'x' and 'y' must be specified together; missing one (x={x}, y={y}).")
			
		
		if points is None:
			if not (-1 <= x <= 1) and not (-1 <= y <= 1):
				raise ValueError("Values of 'x' and 'y' must be in the range [-1,1].")

		if points is not None:
			ang = np.zeros([len(points),2], dtype=np.double)
			idx = 0
			for p in points:
				ang[idx,:] = q3c.xy2ang(p[0], p[1], facenum)
				idx +=1 # seriously, fuck Numpy
			return ang
		else:
			return q3c.xy2ang(x, y, facenum)
	
	def xy2ipix(self, facenum:int=None, x:float=None, y:float=None) -> int: #, points:collections.abc.Iterable):
		'''
		Convert an x,y coordinate pair on the given face number to the ipix value.
		
		:param facenum: the number of the cube face (0=top, 1-4=sides, 5=bottom)
		:param x: x coordinate on cube face
		:param y: y coordinate on cube face
		'''
		ra,dec = self.xy2ang(x=x, y=y, facenum=facenum)
		return self.ang2ipix(ra, dec)
	
	@property
	def nside(self) -> int:
		'''
		Number of bins along one edge of the cube face.
		'''
		return q3c.nside(self._hprm)
	
	@property
	def nbins(self) -> int:
		'''
		Returns the total number of bins in this pixellation scheme (number of divisions per face * 6).
		'''
		return int(pow(self.nside,2)) * 6
	
	def face_number(self, ra:float, dec:float, ipix:int=None) -> int:
		'''
		Return the cube face number for the provided coordinate.
		
		Cube faces are numbered as the following:
		
		|  face 0 = top
		|  face 1 = -45° ≤ α < 45
		|  face 2 = 45° ≤ α < 135
		|  face 3 = 135° ≤ α < 225
		|  face 4 = 225° ≤ α < 315
		|  face 5 = bottom
		
		:param ra: right ascension in degrees; must be in [0,360]
		:param dec: declination in degrees; must be in [-90,90]
		:param ipix: the pixel identifier number
		:returns: the cube face number
		'''
		# be careful with all() and any() as 0 is a legitimate value
		
		if ipix and True in [x is not None for x in [ra,dec]]:
			raise ValueError("Only specify 'ra','dec' OR 'ipix', not all three.")
		if (ra is not None or dec is not None) and (not all([x is not None for x in [ra,dec]])):
			raise ValueError(f"Both parameters 'ra' and 'dec' are required when one is provided (ra={ra}, dec={dec}).")
		if ipix:
			if not (0 <= ipix < self.nbins):
				raise ValueError(f"The ipix number provided ('{ipix}') is outside the range [0,{self.nbins-1}].")
		else:
			# ra,dec provided
			if not -90 <= dec <= 90:
				raise ValueError(f"The value of 'dec' must be normalized to [-90,90]; was given '{dec}'.")
			if not -360 <= ra <= 360:
				raise ValueError(f"The value of 'ra' must be normalized to [0,360]; was given '{ra}'.")

		if ipix:
			ra,dec = self.ipix2ang(ipix)

		return q3c.facenum(self._hprm, ra, dec)
	
	def ipix_area(self, ipix:int, depth:int) -> float:
		'''
		Return the area of a given QLSC pixel in steradians for a given ipix.
		
		NOTE: currently this method returns the average bin size and ignores the ipix and depth values!
		
		Not all pixels are guaranteed to be the same size, but they are very close,
		i.e. the size of any one pixel is a good approximation to any other.
		
		pixel area ≈ 4π / no. pixels on sphere
		
		:param ipix: the pixel identifier number
		:returns: the area of the specified pixel in steradians
		'''
		
		return 4. * math.pi / self.nbins
		
		#if not (0 <= depth <= 30):
		#	raise ValueError(f"The depth provided {depth} is out of the range 0-30")
		
		# The paramter call for the underlying C function is q3c_pixarea(hprm, ipix, depth).
		# It's not clear why "depth" is a parameter since it should be read from hprm.
		# I left the parameter exposed in the Python wrapper, but in this method will just pass bin_level.
		
		#return q3c.pixarea(self._hprm, ipix, depth)
	
	def ang2ipix_at_depth(self, ra:float, dec:float, depth:int) -> int:
		'''
		Return the ipix value at the specified coordinate for any given bin level.
		Returns the ipix value of the pixel center at certain pixel depth covering the specified (ra,dec). [??]
	
		:param ra: right ascension (degrees); must be in [0,360]
		:param dec: declination (degrees); must be in [-90,90]
		:param depth: ; must be in the integer range [1-30]
		:returns: 
		'''
		if not (1 <= depth <= 30):
			raise ValueError(f"The depth provided {depth} is outside the integer range [1,30].")
		return ( q3c.ang2ipix(self._hprm, ra, dec) >> (2*depth) << (2*depth) ) + (1 << (2 * (depth-1))) - 1

	def ipix2polygon(self, ipix:int=None, duplicate_endpoint=False):
		'''
		Returns the points that describe the polygon that define this pixel.
		
		All lines between each coordinate in the polygon should be drawn on great circles.
		
		:param ipix: the pixel identifier number
		:param duplicate_endpoint: if True, repeat the first coordinate in the polygon as the last element
		:returns: an array of ra,dec coordinates in degrees
		'''
		if ipix is None:
			raise ValueError(f"An 'ipix' value must be specified.")
		if not (0 <= ipix < self.nbins):
			raise ValueError(f"The ipix value given ('{ipix}') is out of the valid integer range for this scheme: [0,{self.nbins-1}].")
		
		# array of bin edges, both x and y
		#bins = np.linspace(-1,1,endpoint=True, dtype=np.double, num=self.nside*4+1)
		
		d = 2./self.nside # bin_length since -1 ≤ x ≤ 1 (d=delta)
		
		facenum, x, y = self.ipix2xy(ipix) # returns lower left corner
		
		if duplicate_endpoint:
			polygon = np.array([
					[  x ,  y  ],
					[ x+d,  y  ],
					[ x+d, y+d ],
					[  x , y+d ],
					[  x ,  y  ]
				], dtype=np.double)
		else:
			polygon = np.array([
					[  x ,  y  ],
					[ x+d,  y  ],
					[ x+d, y+d ],
					[  x , y+d ]
				], dtype=np.double)

		if duplicate_endpoint == True:
			polygon = polygon[:-1]
		
		return self.xy2ang(facenum=facenum, points=polygon)
		

class QLSCIndex:
	'''
	
	:param qlsc: an instance of :py:class:`QLSC`
	'''
	def __init__(self, qlsc:QLSC):
		
		if qlsc is None:
			raise ValueError("A QLSC object must be provided to define the pixel resolution of the sphere.")
		
		self.qlsc = qlsc

		self.ra_key = "ra"
		self.dec_key = "dec"
		self.database_tablename = "qlsc_table"
		self.qlsc_schema_version = "1"

		self._data_source = None
		self.data_source = ":memory:" #sqlite3.connect(":memory:")
		
	@property
	def data_source(self):
		'''
		Data structure containing the ra,dec coordinate points.
		'''
		return self._data_source
		
	@data_source.setter
	def data_source(self, new_value):
		'''
		Setter for the :py:func:`data_source` property.
		
		Accepted values:
			* Default value: create an in-memory SQLite database.
			* A path-like object pointing to an SQLite file:
				* If not found, an SQLite file will be created.
				* If a file is present, it will be opened.
			* An SQLite URI, which allows additional read options (e.g. read-only mode).
			* TODO: Numpy structured array. This will require the "ra_key" and "dec_key" properties to be set (if not the defaults of "ra", "dec").
			* TODO: astropy.table.Table (https://astropy.readthedocs.io/en/stable/api/astropy.table.Table.html#astropy.table.Table)
			* TODO: pandas.DataFrame (https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html#pandas.DataFrame)
		'''
		
		# .. todo:: when opening a new sqlite database, check the metadata to see if it's ours
		
		if new_value == self._data_source:
			return
			
		if isinstance(new_value, sqlite3.Connection):
			# we're being given a new database - do we have an existing database open?
			if isinstance(self._data_source, sqlite3.Connection):
				self._data_source.close()
		
		# if a string or path-like object, interpret as an SQLite database path
		if new_value == ":memory:" or isinstance(new_value, (str, pathlib.Path)):

			is_new_database = True
			if new_value == ":memory:":
				self._data_source = sqlite3.connect(new_value)
				
			else:
				is_uri = new_value.startswith("file:")
				if is_uri:
					scheme, netloc, path, params, query, fragment = urlparse(new_value)

				if is_uri:
					filepath = path
				else:
					filepath = new_value
				is_new_database = not os.path.exists(filepath)
					
				try:
					self._data_source = sqlite3.connect(new_value, uri=True)
				except sqlite3.OperationalError as e:
					if new_database:
						raise Exception(f"Unable to create database at specified path ('{new_value}').")
					else:
						raise Exception(f"Found file at path '{new_value}', but unable to open as an SQLite database.")
						
			if is_new_database:
				self._init_sqlite_db()
						
# 		elif isinstance(new_value, numpy.ndarray):
# 			self._data_source = new_value
# 		elif astropy_available and isinstance(new_value, astropy.table.Table):
# 			self._data_source = new_value
# 			raise NotImplementedError("not tested")
# 		elif pandas_available and isinstance(new_value, pandas.DataFrame):
# 			raise NotImplementedError("not tested")
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
		assert isinstance(self.data_source, sqlite3.Connection), f"_init_sqlite_db called on an object that isn't an SQLite database connection! {'self.data_source'} (type={type(self.data_source)})"
		
		#self.data_source.create_function("qlsc_ang2ipix", 3, self.qlsc.ang2ipix)

# 		def plus(a,b):
# 			return a+b
		
		#self.data_source.create_function("qlsc_ang2ipix", 2, plus)

		self._data_source.isolation_level = None # require BEGIN/COMMIT for transactions

		with contextlib.closing(self.data_source.cursor()) as cursor:

			cursor.execute(f'''
				CREATE TABLE {self.database_tablename} (
                    ipix INTEGER,
                    {self.ra_key} REAL,
                    {self.dec_key} REAL,
                    key TEXT UNIQUE);''')

			cursor.execute(f'''CREATE INDEX ipix_idx ON {self.database_tablename}(ipix); ''')
                    
            # Create a metadata table - COUNT(*) on SQLite databases always does a full scan,
            # so keep track of rows count here.
            # Could add additional information here.
			cursor.execute(f'''
            	CREATE TABLE qlsc_metadata (
	            	row_count INTEGER,
	            	qlsc_version TEXT,
	            	date_created TEXT
	            );''')
			
			cursor.execute(f'''
				INSERT INTO qlsc_metadata
					(row_count, qlsc_version, date_created) VALUES
					(0, ?, strftime('%Y-%m-%dT%H:%M:%SZ','now'))
				''', (self.qlsc_schema_version))


	def add_point(self, ra:float=None, dec:float=None, key:str=None):
		'''
		Add a single coordinate point to the index.
		
		:param ra: a single right ascension value, in degrees
		:param dec: a single declination value, in degrees
		:param key: a unique identifier for the coordinate (optional)
		'''
		try:
			with self.data_source:
				with contextlib.closing(self.data_source.cursor()) as cursor:
					cursor.execute(f"INSERT INTO {self.database_tablename} (ipix, {self.ra_key}, {self.dec_key}) VALUES (?, ?, ?)",
									(self.qlsc.ang2ipix(ra,dec), ra, dec))
					cursor.execute("UPDATE qlsc_metadata SET row_count = row_count + 1;")
		except sqlite3.IntegrityError:
			raise NotImplementedError()

	def add_points(self, ra:Iterable=None, dec:Iterable=None, points:Iterable[Tuple]=None, key:Iterable[str]=None):
		'''
		Add a collection of points to the index.
		
		:param ra: an Iterable (e.g. NumPy array) of right ascension points in degrees
		:param dec: an Iterable of declination points in degrees
		:param points: an Iterable of tuples of (ra,dec) pairs
		:param key: an Iterable of unique identifiers for each point
		'''
		if all([x is None for x in [ra,dec,points]]):
			raise ValueError("'ra','dec' OR 'points' must be specified.")
		if all([x is True for x in [ra, dec, points]]):
			raise ValueError("Only 'ra','dec' OR 'points' can be specified.")
		elif points is None and any([x is None for x in [ra,dec]]):
			raise ValueError("If 'ra' or dec' is given, then other parameter must also be provided.")
		
		if np.isscalar(ra):
			ra = np.asarray(ra)
#		elif isinstance(ra, list):
#			ra = np.array(list)
			
		if np.isscalar(dec):
			dec = np.asarray(dec)
#		elif isinstance(dec, list):
#			dec = np.array(list)

		if key is not None and key.isscalar(key):
			key = np.asarray(key)
		
#		if pandas_available:
#			if isinstance(ra, pandas.core.series.Series):
#				ra = ra.values # extract ndarray
#			if isinstance(dec, pandas.core.series.Series):
#				dec = dec.values # extract ndarray
#			if isinstance(key, pandas.core.series.Series):
#				key = key.values # extract ndarray
		
		# .. todo: handle if duplicate keys added
		
		try:
#			row_count = None
			#count = 0
			#max_insert_per_transaction = 1e7

			ipix = self.qlsc.ang2ipix(ra, dec)
			print(f"ipix = {ipix}")

			with self.data_source:
				with contextlib.closing(self.data_source.cursor()) as cursor:
					
					cursor.execute("BEGIN")
					if key is None:
						query = f"INSERT INTO {self.database_tablename} (ipix, {self.ra_key}, {self.dec_key}) VALUES (?, ?, ?)"
						cursor.executemany(query, zip((ipix, ra, dec))) # pass list of tuples (e.g. use zip)
					else:
						query = f"INSERT INTO {self.database_tablename} (ipix, {self.ra_key}, {self.dec_key}, key) VALUES (?, ?, ?, ?)"
						cursor.executemany(query, zip((ipix, ra, dec, key)))
					
# 					if ra is not None:
# 						row_count = len(ra)
# 						ipix = self.qlsc.ang2ipix(ra, dec)
# 						cursor.executemany(query, (ipix, ra, dec)) # pass list of tuples (e.g. use zip)
# #						for i in range(len(ra)):
# #							cursor.execute(query, (self.qlsc.ang2ipix(ra[i], dec[i]), ra[i], dec[i]))
# #							count += 1
# #							if count > 1e7:
# #								self.data_source.commit()
# #								count = 0
# 					else:
# 						row_count = len(points)
# 						for ra,dec in points:
# 							cursor.execute(query, (self.qlsc.ang2ipix(ra, dec), ra, dec))
# 							#count += 1
# 							#if count > 1e7:
# 							#	self.data_source.commit()
# 							#	count = 0

					row_count = len(ipix)
							
					assert row_count is not None, "forgot to set row count"
					print(f"row_count={row_count}")
					cursor.execute("UPDATE qlsc_metadata SET row_count = row_count + ?", (row_count,)) # single element tuple
					cursor.execute("COMMIT")

		except sqlite3.IntegrityError:
			raise NotImplementedError()
	
	@property
	def number_of_points(self) -> int:
		''' Return the number of points in the index. '''
		if isinstance(self.data_source, sqlite3.Connection):
			with contextlib.closing(self.data_source.cursor()) as cursor:
				#cursor.execute(f"SELECT COUNT(*) FROM {self.database_tablename}")
				cursor.execute(f"SELECT row_count FROM qlsc_metadata")
				row_count = cursor.fetchone()[0]
				print(f"row_count: {row_count}")
				return row_count
		else:
			raise NotImplementedError()
		
	
	def radial_query(self, ra:float, dec:float, radius:float, return_key:bool=False) -> np.ndarray: #, idx:bool=False):
		'''
		Given an ra,dec coordinate and a radius (all in degrees), return the points that fall in the cone search.
		
		.. todo:: create option to return an identifier
		
		:param ra: right ascension (degrees)
		:param dec: declination (degrees)
		:param radius: radius (degrees)
		:param return_key: if set to True, returns the key value as provided when added to the index
		:returns: an array of matches, shape (n,2)
		'''
#		:param center_ra:  (degrees)
#		:param center_dec:  (degrees)
#		:param idx: if True, returns an array of Boolean values whether each element matches or not
		
		if abs(dec) > 90:
			raise ValueError(f"The value for dec must be in the range [-90,90]; was given '{dec}'.")
		
		matches = list()
		#matches = np.zeros((len(data), 2)) # fill with boolean values
											# - problem is this can get huge, also, we're not necessarily working with ndarrays
		
		center_ra = ra
		center_dec = dec
		
		if isinstance(self.data_source, sqlite3.Connection):

# 			with contextlib.closing(self._data_source.cursor()) as cursor:
# 				for ra, dec in cursor.execute(f"SELECT {self.ra_key}, {self.dec_key} FROM qlsc_table"):
# 					if self._radial_match(ra, dec, center_ra, center_dec, radius):
# 					   matches.append((ra,dec))
# 				print(matches)
# 				matches = list()
# 					   #matches.append(idx)
# 			print("","===================","")

			#np.array(cursor.fetchall())

			fulls, partials = q3c.radial_query(self.qlsc._hprm, center_ra, center_dec, radius)
			ipix_statements = list()
			for min_ipix,max_ipix in fulls:
				ipix_statements.append(f"(ipix>={min_ipix} AND ipix<{max_ipix})")
			for min_ipix,max_ipix in partials:
				ipix_statements.append(f"(ipix>={min_ipix} AND ipix<{max_ipix})")
			wheres = "({0})".format(" OR ".join(ipix_statements))
			query = f"SELECT ra,dec,key FROM {self.database_tablename} WHERE {wheres}"
			#print(query)
			cone_radius = pow(sin(deg2rad(radius)/2.), 2)
			with contextlib.closing(self._data_source.cursor()) as cursor:
				print(cursor.execute(query))
				for ra, dec,key in cursor.execute(query):
					# filter out points outside radius
					if sindist(ra, dec, center_ra, center_dec) < cone_radius:
						matches.append((ra,dec,key))
		
# 		elif isinstance(self.data_source, numpy.ndarray):
# 			data = self.data_source
# 			for i in range(len(data)):
# 				ra = data[ra_key][i]
# 				dec = data[dec_key][i]
# 				if self._radial_match(ra, dec, center_ra, center_dec, radius):
# 					#matches.append((ra,dec))
# 					matches.append(i)
# 
# 		elif pandas_available and isinstance(self.data_source, pandas.DataFrame):
# 			data = self.data_source
# 			for i in range(len(data)):
# 				ra = data[ra_key][i]
# 				dec = data[dec_key][i]
# 				if self._radial_match(ra, dec, center_ra, center_dec, radius):
# 					#matches.append((ra,dec))
# 					matches.append(i)
# 		
# 		elif astropy_available and isinstance(self.data_source, astropy.table.Table):
# 			data = self.data_source
# 			for i in range(len(data)):
# 				ra = data[ra_key][i]
# 				dec = data[dec_key][i]
# 				if self._radial_match(ipix, ra, dec, center_ra, center_dec, radius):
# 					#matches.append((ra,dec))
# 					matches.append(i)
		
#		if idx:
#			# return list of indices that match
#			print(f"matches = {matches}")
#			return np.array(matches)
#		else:
			# return array of matches
#			return np.take(data, matches) # select the elements that match

		return np.array([matches])

	def _radial_match(self, ra, dec, center_ra, center_dec, radius):
		'''
		Returns 'True' if the provided ipix value is within radius of (center_ra,center_dec) (all in degrees).
		'''
		# indexing makes this step fast - this is the slow part
		# index is ra,dec in, ipix out
		ipix = self.qlsc.ang2ipix(ra, dec)
		hprm = self.qlsc._hprm
# 		print(f"iter={0}, full=1 : {radial_query_it(hprm,center_ra,center_dec,radius,0,1)}")
# 		print(f"iter={1}, full=1 : {radial_query_it(hprm,center_ra,center_dec,radius,1,1)}")
# 		print(f"iter={0}, full=0 : {radial_query_it(hprm,center_ra,center_dec,radius,0,0)}")
# 		print(f"iter={1}, full=1 : {radial_query_it(hprm,center_ra,center_dec,radius,1,0)}")
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
			
