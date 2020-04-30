
import os
import sys
import math
import pathlib
import sqlite3
import logging
import contextlib
from typing import Iterable, Tuple, Union
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

# functions defined in "qlsc_c_module.c" as named in the array "qlsc_methods" are in the qlsc.q3c namespace
from . import q3c
from .q3c import radial_query_it, sindist

from .utilities import _normalize_ang

# really important for the ipix values!
sqlite3.register_adapter(np.int64, int)

logger = logging.getLogger("qlsc_logger")

class QLSC:
	'''
	A class that describes a quadrilaterized spherical cube (QLSC).

	The quadrilaterized sphere is divided into six faces. The bin level defines the number times a
	face is subdivided. A bin level of 1 divides each face into four "pixels"; bin level two divides
	each of those pixels again into four, etc.

	At the lowest resolution where bin level=0 there is 1 bin per face for a total of 6 bins.
	The next bin level 1 divides each face into four bins, i.e. a total of 24 bins.
	Higher bin levels further divide each bin four more times.
	
	A scheme can be defined by the bin level where:

		.. code-block:: 

			no. bins per face        = (2**depth)**2
			number of bins on sphere = no. bins per face * 6
	
	Each bin has a pixel number called an "ipix", which is an integer encoded with the location of the bin.
	
	The default value of `depth` = 30 that divides the sphere into 6,442,450,944 bins (i.e. nside=1,073,741,824),
	is the one used by the `Q3C PostgreSQL plugin <https://github.com/segasai/q3c>`_.
	This value is used if neither `nside` or `bin_size` are specified.
	
	:param depth: number of times each bin in each face is subdivided by four
	'''
	def __init__(self, depth:int=30):

		# k = depth in Chan paper

		# The Q3C code only supports up to depth=30
		if not (0 <= depth <= 30):
			raise ValueError(f"The value for 'depth' must be an integer on [0,30]; was given '{depth}'.")
		
		# "nside" as defined by Q3C is the number of bins along one side of the cube face.
		# This means the number of bins on a cube face is nside*nside, and the total
		# number of bins is 6*nside*nside.
		#
		nside = int(pow(2, depth))
		
		self.depth = depth
		
		# store hprm struct values - hprm is a pointer (struct) to main Q3C structure used by most Q3C calls
		self._hprm = q3c.init_q3c(nside=nside)
		if self._hprm is None:
			raise Exception("Internal Q3C data structure could not be created.")
			
	def __repr__(self):
		return f"<{self.__class__.__module__.split('.')[0]}.{self.__class__.__name__} at {hex(id(self))}, depth={self.depth}>"
	
	def ang2ipix(self, ra:Union[float,Iterable]=None, dec:Union[float,Iterable]=None, points:Iterable=None) -> Union[float,np.ndarray]:
		'''
		Convert ra,dec positions to an ipix number.
		
		This method accepts scalar values or NumPy arrays. Slices of any other NumPy array are accepted.
		Note that arrays are required to be ``dtype.double``; copies will be made of arrays of any other type.
		
		Right ascension values outside of the range [0,360] are handled efficiently.
		Declination values outside of [-90,90] are handled, but may result in arrays being copied. (This is only
		important when working with extremely large arrays.)
		
		Example usage:
		
		.. code-block:: python
		
			ipix = ang2ipix(45,60)
			ipix = ang2ipix(ra=45,dec=60)
			ipix = ang2ipix(points=np.array([[45,60],[0,55]])
			ipix = ang2ipix(ra=catalog[:,4], dec=catalog[:,5])
		
		:param ra: right ascension (degrees)
		:param dec: declination (degrees)
		:param points: an array of [ra,dec] points, shape=(n,2)
		:returns: the ipix integer if a single point is provided, an array of ipix integers if an array of values is provided
		'''
		# check "ra","dec" are both set
		if points is None:
			if True in [x is None for x in [ra,dec]]: # can't use any/all with possible "0" values
				raise ValueError(f"Values for both 'ra' and 'dec' must be specified.")
		else:
			if True in [x is not None for x in [ra,dec]]:
				raise ValueError("The parameter 'points' was given; 'ra' or 'dec' should not also be specified.")

		# q3c.ang2ipix unwraps ra values to [0,360], but it does not do the same for dec.
		# It will silently truncate to [-90,90].
		
		# handle simple scalar values
		#		
		if np.isscalar(ra) and np.isscalar(dec):
			if not (-90 <= dec <= 90):
				ra, dec = _normalize_ang(ra,dec)
			return q3c.ang2ipix(self._hprm, ra, dec)

		if points is not None:
			# break it down:
			ra = points[:,0]
			dec = points[:,1]

		if not ((-90 <= dec).all() and (dec <= 90).all()):
			ra, dec = _normalize_ang(ra=ra, dec=dec, copy=True)

		return q3c.ang2ipix(self._hprm, ra, dec)

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
			
	def ipix_down(self, ipix:int) -> Iterable[int]:
		'''
		Returns the four ipix values at the next higher depth.
		
		:param ipix: the ipix value at the resolution of this object's scheme.
		'''
		if self.depth == 30:
			raise ValueError("QLSC doesn't support depth values greater than 30.")
		elif not (0 <= ipix < self.nbins):
			raise ValueError(f"ipix value out of range for this depth; should be in [0,{self.nbins}).")
		else:
			new_ipix = ipix <<  2 # same as -> int(ipix * 2**(self.nside))
			return list(range(new_ipix, new_ipix+4))
		
	def ipix_up(self, ipix:int):
		'''
		Returns the ipix value value at the next lower depth.

		:param ipix: the ipix value at the resolution of this object's scheme.
		'''
		if self.depth == 0:
			raise ValueError("Already at lowest depth (0).")
		elif not (0 <= ipix < self.nbins):
			raise ValueError(f"ipix value out of range for this depth; should be in [0,{self.nbins}).")
		else:
			return ipix >> 2
	
	def ipix2face(self, ipix:int) -> int:
		'''
		Return the face number the given ipix value falls on.
		'''
		if not (0 <= ipix < self.nbins):
			raise ValueError(f"ipix value out of range for this depth; should be in [0,{self.nbins}).")
		return ipix // math.pow(self.nside, 2)
	
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
		# I left the parameter exposed in the Python wrapper, but in this method will just pass depth.
		
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
	An object that stores (ra,dec) coordinates on a sphere and performs fast cone searches.
	
	The ``QLSCIndex`` object can be provided with a file path to save the index for repeated use, e.g.
	
	.. code-block:: python
		
		qlsc = QLSC(depth=30)
		sdss_catalog = QLSCIndex(qlsc=qlsc, filepath='/path/to/save/file.qlscidx')
		
	The filename and extension can be anything, but `.qlscidx` is suggested to remember what kind of file it is.
	If no ``filepath`` parameter is provided, an index is created in memory, but note that it will
	no longer exist after the program ends. If an existing file is found at the file path provided,
	an attempt will be made to open it and, if a valid ``QLSCIndex`` file, it will be used.
	The file path can be provided in the form of a URI, e.g. ``file:///Users/meerkat/catalog.qlscidx``.
	
	The underlying storage is an SQLite database, but this should be considered an implementation
	detail and can change. (This may be documented for expanded use in the future.)
	
	:param qlsc: an instance of :py:class:`QLSC` set to the desired segmentation level
	:param index_file: the file path where a persistent index will be stored; the default is to create the index in memory
	'''
	def __init__(self, qlsc:QLSC=None, filepath:os.PathLike=None):
		
		if qlsc is None:
			# select the highest resolution by default
			qslc = QLSC(depth=30)
		
		self.qlsc = qlsc
		self._db = None # SQLite database connection
		
		self.database_tablename = "qlsc_ipix"
		self.qlsc_schema_version = "1"
				
		if filepath:
			self.index_path = filepath # interpret "None" to use an in memory database
		else:
			# allows multiple connections to the in-memory database
			# see: https://sqlite.org/inmemorydb.html
			self.index_path = ":memory:" # this just keeps creating a file on disk -> "file::memory:?cache=shared"
		
		self._initial_database_connection()

	def _initial_database_connection(self):
		'''
		Make the initial connection to the database.
		'''

		is_new_database = True

		if self._is_in_memory_db():
			# Always keep this connection open to make sure the database is not released
			# i.e., set it here and never touch it again.
			# Use "self._db" to open/close connections as one would with a file-based db.
			# 
			# ^--- shared cache plan, but didn't get that to work
			#
			#self._in_memory_db_connection = sqlite3.connect(self.index_path)
			self._db = sqlite3.connect(self.index_path) #opens a second connection that can be closed without losing the database
		else:
				
			is_uri = self.index_path.startswith("file:")
			if is_uri:
				scheme, netloc, path, params, query, fragment = urlparse(self.index_path)
				self.index_path = path

			is_new_database = not os.path.exists(self.index_path)
			#logger.debug(f"index path = {self.index_path}")
				
			try:
				self._db = sqlite3.connect(self.index_path) #, uri=True)
			except sqlite3.OperationalError as e:
				if new_database:
					raise Exception(f"Unable to create database at specified path ('{self.index_path}').")
				else:
					raise Exception(f"Found file at path '{self.index_path}', but am unable to open as an SQLite database.")
			
		self._configure_db_connection(self._db)

		if is_new_database:
			self._init_sqlite_db()
		else:
			db_ok = self._validate_db_as_qlsc_index()
			if not db_ok:
				self._db.close()
				raise Exception(f"The file provided is an SQLite database (good) but not in a format I recognize (bad).")
	
	def _validate_db_as_qlsc_index(self):
		'''
		Validates that the database found is one of ours.
		'''
		# keep this simple for the time being
		try:
			with contextlib.closing(self._db.cursor()) as cursor:
				cursor.execute("SELECT * FROM qlsc_metadata")
				metadata = cursor.fetchone()
				logger.debug(metadata)
				return True
		except sqlite3.OperationalError as e:
			if "no such table" in str(e):
				return False
	
	def __del__(self):
		''' Destructor method. '''
		# if we had an open SQLite connection, close it
		#if isinstance(self._data_source, sqlite3.Connection):
		if self._db:
			self._db.close()
	
	def _is_in_memory_db(self):
		return (self.index_path == ":memory:") or self.index_path.startswith("file::memory:")

	def _configure_db_connection(self, dbconn):
		'''
		Configure connection-level settings on SQLite database.
		'''
		# set database-specific settings
		dbconn.isolation_level = None # require BEGIN/COMMIT for transactions
		dbconn.row_factory = sqlite3.Row
	
	def _init_sqlite_db(self):
		'''
		Initialize a new index database, e.g. create schema, initialize metadata.
		'''
		
		#self.sqlite_db = sqlite3.connect(':memory:')
		#con.create_function("ang2ipix", 1, self.ang2ipix)
		
		#self.data_source.create_function("qlsc_ang2ipix", 3, self.qlsc.ang2ipix)

# 		def plus(a,b):
# 			return a+b
		
		#self.data_source.create_function("qlsc_ang2ipix", 2, plus)

		with contextlib.closing(self._db.cursor()) as cursor:

			cursor.execute(f'''
				CREATE TABLE {self.database_tablename} (
                    ipix INTEGER,
                    ra REAL,
                    dec REAL,
                    key TEXT,
                    UNIQUE (ra, dec)
                );''')

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
			
			# set initial metadata values
			cursor.execute(f'''
				INSERT INTO qlsc_metadata
					(row_count, qlsc_version, date_created) VALUES
					(0, ?, strftime('%Y-%m-%dT%H:%M:%SZ','now'));
				''', (self.qlsc_schema_version))

			# internally keep count of records since a count(*)
			# in SQLite is a full table scan
# 			cursor.execute(f'''
# 				CREATE TRIGGER row_count_trigger
# 					AFTER INSERT ON {self.database_tablename}
# 					BEGIN
# 						UPDATE qlsc_metadata SET row_count = row_count + 1;
# 					END;
# 				''');

		# Close and reopen the database to reset the "connection.total_changes" value to 0,
		# used to count number of changed rows after inserts.
		# If this is an in-memory database, open a second connection first to make sure that
		# the memory isn't freed (happens when last connection is closed).
# 		if self._is_in_memory_db:
# 			tmp_var = sqlite3.connect(self.index_path)
# 			self._db.close()
# 			self._db = tmp_var
# 		else:
		#self._db.close()
		#self._db = sqlite3.connect(self.index_path)
		#self._configure_db_connection(self._db)

	def add_point(self, ra:float=None, dec:float=None, key:str=None):
		'''
		Add a single coordinate point to the index.
		
		This method can be called multiple times with the same coordinate, but duplicate
		entries will be discarded. Key values are not required to be unique. If a coordinate
		is added to the index without a key and again later with a key, it will not be added
		to the index (and vice versa).
		
		:param ra: a single right ascension value, in degrees
		:param dec: a single declination value, in degrees
		:param key: a unique identifier for the coordinate (optional)
		'''
		
		if any([x is None for x in [ra,dec]]):
			raise ValueError("Both 'ra' and 'dec' values must be specified.")
		
		# Q3C.ang2ipix truncates dec outside of [-90,90]
		if not (-90 <= dec <= 90):
			ra, dec = _normalize_ang(ra,dec)
		
		ipix = self.qlsc.ang2ipix(ra,dec)
		
		with self._db:
			with contextlib.closing(self._db.cursor()) as cursor:
				
				cursor.execute("BEGIN")
				if key is None:
					query = f"INSERT OR IGNORE INTO {self.database_tablename} (ipix, ra, dec) VALUES (?, ?, ?)"
					cursor.execute(query, (ipix, ra, dec)) # pass list of tuples (e.g. use zip)
				else:
					query = f"INSERT OR IGNORE INTO {self.database_tablename} (ipix, ra, dec, key) VALUES (?, ?, ?, ?)"
					cursor.execute(query, (ipix, ra, dec, key))
				cursor.execute("COMMIT")

	def add_points(self, ra:Iterable=None, dec:Iterable=None, points:Iterable[Tuple]=None, keys:Iterable[str]=None):
		'''
		Add a collection of points to the index.
		
		This method can be called multiple times with the same coordinate, but duplicate
		entries will be discarded. Key values are not required to be unique. If a coordinate
		is added to the index without a key and again later with a key, it will not be added
		to the index (and vice versa).

		This method is preferred for a large number of points over using :py:func:`add_point` in a loop.

		Note that if the ``dec`` values are not in the range [-90,90], a copy of the ra,dec arrays will
		be made to normalize the points first. This may only be an issue when providing an extremely
		large number of points.
		
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
		
		if keys is not None:
			if points is not None:
				if len(points) != len(keys):
					raise ValueError(f"The number of points provided ({len(points)}) does not match the number of keys ({len(keys)}).")
			else:
				if len(ra) != len(dec) != len(keys):
					raise ValueError(f"The number of ra,dec coordinates provided ({len(ra)},{len(dec)}) does not match the number of keys ({len(keys)}).")

#		if pandas_available:
#			if isinstance(ra, pandas.core.series.Series):
#				ra = ra.values # extract ndarray
#			if isinstance(dec, pandas.core.series.Series):
#				dec = dec.values # extract ndarray
#			if isinstance(key, pandas.core.series.Series):
#				key = key.values # extract ndarray
		
		if points is not None:
			# break it down
			ra = points[:,0]
			dec= points[:,1]

		if not ((-90 <= dec).all() and (dec <= 90).all()):
			points = _normalize_ang(points=points, copy=True)
		
		try:
			ipix = self.qlsc.ang2ipix(ra, dec)

			with self._db:
				with contextlib.closing(self._db.cursor()) as cursor:
					
					cursor.execute("BEGIN")
					if keys is None:
						query = f"INSERT OR IGNORE INTO {self.database_tablename} (ipix, ra, dec) VALUES (?, ?, ?)"
						# use individual lists instead of something like np.dstack- we don't want
						# the 64-bit ipix integers converted to floats and risk changing the value due to rounding
						cursor.executemany(query, zip(ipix, ra, dec)) # pass list of tuples (e.g. use zip)
					else:
						query = f"INSERT OR IGNORE INTO {self.database_tablename} (ipix, ra, dec, key) VALUES (?, ?, ?, ?)"
						cursor.executemany(query, zip(ipix, ra, dec, keys))
					
					cursor.execute("COMMIT")

		except sqlite3.IntegrityError:
			raise NotImplementedError()
	
	@property
	def number_of_points(self) -> int:
		''' Return the number of coordinate points in the index. '''
		with contextlib.closing(self._db.cursor()) as cursor:
			cursor.execute(f"SELECT max(rowid) FROM {self.database_tablename}")
			return cursor.fetchone()[0]
	
	def radial_query(self, ra:float, dec:float, radius:float, return_key:bool=False):
		'''
		Given an ra,dec coordinate and a radius (all in degrees), return the points that fall in the cone search.
		
		:param ra: right ascension (degrees)
		:param dec: declination (degrees)
		:param radius: radius (degrees)
		:param return_key: if set to True, returns the key value as provided when added to the index
		:returns: if ``return_key=True``, returns a NumPy record array with keys ``ra,dec,key``; otherwise an array of matches, shape (n,2)
		'''
		
		if abs(dec) > 90:
			raise ValueError(f"The value for dec must be in the range [-90,90]; was given '{dec}'.")
		
		center_ra = ra
		center_dec = dec
		
		# gather matches into these lists
		match_ra = list()
		match_dec = list()
		match_key = list()
		
		if isinstance(self._db, sqlite3.Connection):

			fulls, partials = q3c.radial_query(self.qlsc._hprm, center_ra, center_dec, radius)

			ipix_statements = list()
			for min_ipix,max_ipix in fulls:
				ipix_statements.append(f"(ipix>={min_ipix} AND ipix<{max_ipix})")
			for min_ipix,max_ipix in partials:
				ipix_statements.append(f"(ipix>={min_ipix} AND ipix<{max_ipix})")
			wheres = "({0})".format(" OR ".join(ipix_statements))
			query = f"SELECT ra,dec,key FROM {self.database_tablename} WHERE {wheres}"

			cone_radius = pow(sin(deg2rad(radius)/2.), 2)

			with contextlib.closing(self._db.cursor()) as cursor:
				for ra, dec,key in cursor.execute(query):
					# filter out points outside radius
					if sindist(ra, dec, center_ra, center_dec) < cone_radius:
						if return_key:
							match_ra.append(ra)
							match_dec.append(dec)
							match_key.append(key)
						else:
							match_ra.append(ra)
							match_dec.append(dec)

		if return_key:
			return np.core.records.fromarrays([match_ra, match_dec, match_key], names='ra,dec,key')
		else:
			return np.squeeze(np.dstack((match_ra,match_dec)))

