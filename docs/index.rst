.. QLSC Python Module documentation master file, created by
   sphinx-quickstart on Wed Apr 22 02:59:23 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

QLSC: Quadrilateralized Spherical Cube for Python
=================================================

The quadrilateralized spherical cube (QLSL) is a geospatial indexing
scheme for segmenting a sphere into pixels with the aim of optimized
spatial indexing and queries. *qlsl* is an implemenation of this scheme
in a Python package. Parts of it are based on code from Sergey Koposov’s
`Q3C <https://github.com/segasai/q3c>`_, a PostgreSQL extension that
implements QLSC indexing.

.. image:: _static/cube_subdivisions.png
    :align: center

Note that while this package is designed for astronomical use (it focusses
on right ascension and declination), it could be just as easily be used for
latitude and longitude coordinates, as long as you’re ok with a perfectly
spherical Earth (though QLSL was designed to be used with the real Earth).
Future updates may better facilitate this, but contributions are welcome.

Please refer to `this project's GitHub page <https://github.com/demitri/qlsc>`_ 
for a description of the QLSC scheme and how to use the code.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

API Reference
-------------

.. toctree::
   :maxdepth: 2
   
   api

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
