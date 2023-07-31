#
#from setuptools import setup, Extension, find_packages
#from setuptools.command.build_ext import build_ext

import re
import setuptools
from distutils.core import setup, Extension
#from setuptools import setup, Extension
#from setuptools import find_packages
from setuptools.command.build_ext import build_ext as _build_ext
#from setuptools import dist

#dist.Distribution().fetch_build_eggs(['numpy>=1.19'])

class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)

        import numpy as np
        self.include_dirs.append(np.get_include())


sources = ["qlsc/c_code/qlsc_c_module.c", "qlsc/c_code/q3c/q3cube.c", "qlsc/c_code/q3c/q3c_poly.c"]
data_files = []
include_dirs = ['q3c']		# -I directories
library_dirs = []			# -L directories
libraries = []		# libraries to include
define_macros = [('Q3C_VERSION', '"2.0.0"'), # equivalent of a bare #define in C source
				 ('NPY_NO_DEPRECATED_API', 'NPY_1_18_API_VERSION')] # warn if using a deprecated NumPy API, defined in numpy/numpyconfig.h
extra_link_args = [] # e.g. ['-framework', 'OpenGL', '-framework', 'GLUT'])

# The name of the module being built is "q3c".
# Place it under the "qlsc" package below in setup().
# Ref: https://docs.python.org/3/distutils/setupscript.html#extension-names-and-packages
c_extension = Extension(name="q3c",
						sources=sources,
						language='C',
						include_dirs=include_dirs,
						library_dirs=library_dirs,
						define_macros=define_macros,
						libraries=libraries,
						extra_compile_args=["-std=c99"])

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
    description=description,
    long_description=long_description,
    #long_description_content_type='text/markdown; charset=UTF-8; variant=GFM',
    #license="GPL",
    #classifiers=classifiers,
    url="https://github.com/demitri/qlsc",
    project_urls={
    	"Documentation" : "https://qlsc.readthedocs.io/en/latest/",
    	"Source Code" : "https://github.com/demitri/qlsc"
	},
    #project_urls={
    #        "Bug Tracker": "https://bugs.example.com/HelloWorld/",
    #        "Documentation": "https://docs.example.com/HelloWorld/",
    #        "Source Code": "https://code.example.com/HelloWorld/",
    #    },
    author="Demitri Muna",
    author_email="demitri@trillianverse.org",
    setup_requires=['wheel', 'numpy'], # needed to package for distribution
    install_requires=['numpy'],
    #packages=find_packages(),
    data_files=data_files,
    ext_package="qlsc", # will compile the methods from the extension to the namespace "qlsc"
    ext_modules=[c_extension], # alternative: cythonize(etc), needs "from Cython.Build import cythonize"
    include_dirs=['q3c'],
    packages=setuptools.find_packages(), #['qlsc'],
    include_package_data = True, # add files specified in MANIFEST.in, specifically header files
    python_requires='>=3.6',
    cmdclass={"build_ext": build_ext}
)
