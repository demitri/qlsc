#
#from setuptools import setup, Extension, find_packages
#from setuptools.command.build_ext import build_ext

import re
from distutils.core import setup, Extension
#from setuptools import setup, Extension
#from setuptools import find_packages

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

sources = ["pyq3c/pyq3c_module.c", "q3c/q3cube.c", "q3c/q3c_poly.c"]
data_files = []
include_dirs = ['q3c']		# -I directories
library_dirs = []			# -L directories
libraries = []		# libraries to include
define_macros = [('Q3C_VERSION', '"1.8.1"')] # equivalent of a bare #define in C source
extra_link_args = [] # e.g. ['-framework', 'OpenGL', '-framework', 'GLUT'])

ext = Extension(name="pyq3c._q3c_wrapper",
				sources=sources,
				language=['C'],
				include_dirs=include_dirs,
				library_dirs=library_dirs,
				define_macros=define_macros,
				libraries=libraries)

description = ("A Python library around the Q3C code.")
long_description = "-- add long description here --"

classifiers = [
    "Development Status :: 2 - Pre-Alpha",
    "License :: OSI Approved :: GNU General Public License (GPL)",
    "Topic :: Scientific/Engineering :: Astronomy",
    "Intended Audience :: Science/Research"
]

exec(open('pyq3c/version.py').read())
setup(
    name="pyq3c",
    version=__version__,
    #version=get_property('__version__', 'pyq3c'),
    description=description,
    long_description=long_description,
    #long_description_content_type='text/markdown; charset=UTF-8; variant=GFM',
    #license="GPL",
    #classifiers=classifiers,
    url="https://github.com/demitri/pyq3c",
    author="Demitri Muna",
    author_email="demitri@scicoder.org",
    #setup_requires=[],
    #install_requires=[],
    #packages=find_packages(),
    data_files=data_files,
    ext_modules=[ext], # alternative: cythonize(etc), needs "from Cython.Build import cythonize"
    include_dirs=['q3c'],
    packages=['pyq3c']
)
#    cmdclass={"build_ext": build_ext_subclass}
#)
