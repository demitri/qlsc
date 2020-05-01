#
#from setuptools import setup, Extension, find_packages
#from setuptools.command.build_ext import build_ext

import re
import setuptools
from distutils.core import setup, Extension
#from setuptools import setup, Extension
#from setuptools import find_packages

import numpy as np

# directions on including NumPy in a C extension:
# Ref: https://numpy.org/doc/stable/reference/c-api/array.html#importing-the-api

#copt = {}
#lopt = {}

#class build_ext_subclass( build_ext ):
#	pass
#     def build_extensions(self):
#         c = self.compiler.compiler_type
#         if c in copt:
#            for e in self.extensions:
#                e.extra_compile_args = copt[ c ]
#         if c in lopt:
#             for e in self.extensions:
#                 e.extra_link_args = lopt[ c ]
#         build_ext.build_extensions(self)

def get_property(prop:str, project:str):
	'''
	Read the requested property (e.g. '__version__', '__author__') from the specified Python module.
	Ref: https://stackoverflow.com/a/41110107/2712652
	'''
	result = re.search(r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(prop), open(project + '/__init__.py').read())
	return result.group(1)

sources = ["qlsc/c_code/qlsc_c_module.c", "qlsc/c_code/q3c/q3cube.c", "qlsc/c_code/q3c/q3c_poly.c"]
data_files = []
include_dirs = ['q3c', np.get_include()]		# -I directories
library_dirs = []			# -L directories
libraries = []		# libraries to include
define_macros = [('Q3C_VERSION', '"2.0.0"'), # equivalent of a bare #define in C source
				 ('NPY_NO_DEPRECATED_API', 'NPY_1_18_API_VERSION')] # warn if using a deprecated NumPy API, defined in numpy/numpyconfig.h
extra_link_args = [] # e.g. ['-framework', 'OpenGL', '-framework', 'GLUT'])

# The name of the module being build is "q3c".
# Place is under the "qlsl" package below in setup().
# Ref: https://docs.python.org/3/distutils/setupscript.html#extension-names-and-packages
c_extension = Extension(name="q3c",
				sources=sources,
				language='C',
				include_dirs=include_dirs,
				library_dirs=library_dirs,
				define_macros=define_macros,
				libraries=libraries)

description = ("A Python implementation of the quadrilateralized spherical cube scheme.")
long_description = '''The quadrilateralized spherical cube (QLSC) is a geospatial indexing scheme for segmenting a sphere into pixels with the aim of optimized spatial indexing and queries. QLSC is an implementation of this scheme in a Python package. Parts of it are based on code from Sergey Koposov’s Q3C, a PostgreSQL extension that implements QLSC indexing. In addition to sphere segmentation, this package provides the catalog indexing functionality of Q3C without the need to install a PostgreSQL database.

This package also enables fast local cone searches on astronomical catalogs without the need to install any other dependencies, supporting essentially unlimited numbers of coordinates (i.e. if you're familiar with Q3C, the same can be performed without PostgreSQL).'''

# list of classifiers: https://pypi.org/classifiers/
classifiers = [
    "Development Status :: 5 - Production/Stable",
    "License :: OSI Approved :: GNU General Public License v2 (GPLv2)",
    "Topic :: Scientific/Engineering :: Astronomy",
    "Intended Audience :: Science/Research"
]

exec(open('qlsc/version.py').read())
setup(
    name="qlsc",
    version=__version__,
    #version=get_property('__version__', 'qlsc'),
    description=description,
    long_description=long_description,
    #long_description_content_type='text/markdown; charset=UTF-8; variant=GFM',
    #license="GPL",
    #classifiers=classifiers,
    url="https://github.com/demitri/qlsc",
    author="Demitri Muna",
    author_email="demitri@trillianverse.org",
    #setup_requires=['wheel'], # needed to package for distribution
    #install_requires=[],
    #packages=find_packages(),
    data_files=data_files,
    ext_package="qlsc", # will compile the methods from the extension to the namespace "qlsc"
    ext_modules=[c_extension], # alternative: cythonize(etc), needs "from Cython.Build import cythonize"
    include_dirs=['q3c'],
    packages=['qlsc']
)
#    cmdclass={"build_ext": build_ext_subclass}
#)
