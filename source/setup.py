#
#from setuptools import setup, Extension, find_packages
#from setuptools.command.build_ext import build_ext

import re
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
				language=['C'],
				include_dirs=include_dirs,
				library_dirs=library_dirs,
				define_macros=define_macros,
				libraries=libraries)

description = ("A Python library around the Q3C code.")
long_description = "-- add long description here --"

# list of classifiers: https://pypi.org/classifiers/
classifiers = [
    "Development Status :: 4 - Beta",
    "License :: OSI Approved :: GNU General Public License (GPL)",
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
    author_email="demitri@scicoder.org",
    #setup_requires=[],
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
