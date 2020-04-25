# QLSC Python Package

The quadrilateralized spherical cube (QLSL) is a geospatial indexing scheme for segmenting a sphere into pixels with the aim of optimized spatial indexing and queries. *qlsl* is an implemenation of this scheme in a Python package. Parts of it are based on code from Sergey Koposov’s [Q3C](https://github.com/segasai/q3c), a PostgreSQL extension that implements QLSC indexing.

Note that while this package is designed for astronomical use (it focusses on right ascension and declination), it could be just as easily be used for latitude and longitude coordinates, as long as you’re ok with a perfectly spherical Earth (though QLSL was designed to be used with the real Earth). Future updates may better facilitate this, but contributions are welcome.

API documentation can be found here: [https://qlsc.readthedocs.io/en/latest/](https://qlsc.readthedocs.io/en/latest/)

## Installation

#### Install from source

    git clone https://github.com/demitri/qlsc.git
    cd qlsc/source
    python setup.py install
   
## 30 Second Introduction

If you’re not familiar with the QLSC scheme it would be better to start below. There are two classes defined in this package, `QLSC` and `QLSCIndex`. The QLSC scheme projects each of the six faces of a cube onto a sphere. Segmentation is performed on the cube fased and done in levels, where each level divides each bin into four. For example, `bin_level=2` divides each bin into four, then of those into four again, which is 2^(2*`bin_level`) (i.e. 16) bins per cube face, or 96 total bins. ("Pixels" and "bin" are used somewhat interchangably.)

The class can convert between (ra,dec) positions and the pixel number (ipix).

```python
from qlsc import QLSC
    
q = QLSC(bin_level=2)
q.ang2ipix(45,45)     # ra,dec -> ipix number
q.ipix2ang_center(42) # ipix number -> ra,dec at center of pixel
q.ipix2polygon(42)    # returns the points on the sphere describing
                      # the pixel (joined by great circles)
``` 

## The Quadrilateralized Spherical Cube

The quadrilateralized spherical cube was first devised in 1975 in a technical report by F.K. Chan and E.M. O'Neill ([available here](https://ntrl.ntis.gov/NTRL/dashboard/searchResults/titleDetail/ADA010232.xhtml)). It begins with a cube inscribed within a sphere:

![](figures/cube_in_sphere/cube_in_sphere.png)

The six cube faces are then projected onto the sphere via transforms defined in the paper. This is the lowest resolution. Higher resolutions are achieved in steps: each step divides the bin, or pixel, into four. The first step will have four bins per face, the next will have 16, etc. The code refers to these as "bin levels" (number of times the bin has been subdivided), where the diagram above is `bin_level=0`. The Q3C PostgreSQL extension uses a bin level of 30, which is 1,152,921,504,606,846,976 bins per cube face (six times that over the full sphere) corresponding to ~0.08 μ" square per pixel. The QLSC package supports any bin level from 0 to 30.

| Cube Face Number | RA Range | Face Center    |
|:----------------:| :------: | :------------: |
|  0          | top face | α = 0°, δ = 90° |
|  1          | -45° ≤ δ ≤ 45° | α = 0°, δ = 0° |
|  2          | 45° ≤ δ ≤ 135° | α = 90°, δ = 0° |
|  3          | 135° ≤ δ ≤ 225° | α = 180°, δ = 0° |
|  4          | 225° ≤ δ ≤ 315° | α = 270°, δ = 0° |
|  5          | bottom face | α = 0°, δ = -90° |

The diagram shows face 1 divided at `bin_level=2` and the projection of each pixel onto the sphere. Other faces have been hidden for clarity.

![](figures/cube_subdivisions/cube_subdivisions.png)

The distortion may appear to be significant, but at much higher bin levels it is not as severe. Also, it doesn't matter. The scheme is advantageous in that there are no discontinuities at the poles, and the indexing scheme is optimized for fast database queries. While not all pixels are the same size, they are very nearly so, and for the purposes this scheme is used for, this is not a requirement.

For most users, working in ipix values and ra,dec coordinates will accomplish most anything you might need. For those who might be performing more complex calculations, I recommend working in the native coordinates of the face plane. Each face (at any division level) has its 2D coordinate system origin at the cube face center. Both the *x* and *y* axes range from -1 to +1. Methods are provided to translate between *xy* coordinates, ipix value, and ra/dec coordinates. The source code for the [`ipix2polygon()`](https://qlsc.readthedocs.io/en/latest/api.html#qlsc.QLSC.ipix2polygon) method provides an illustrative example.

## References

* [Original Chan, O'Neill 1975 paper](https://ntrl.ntis.gov/NTRL/dashboard/searchResults/titleDetail/ADA010232.xhtml)
* [Abridged F.K. Chan paper - A Quadrilateralized Spherical Cube Earth Data Base](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19810002572.pdf)
* [Koposov, Bartunov Q3C ADASS proceedings paper](https://ui.adsabs.harvard.edu/abs/2006ASPC..351..735K/abstract)
* [Stack Overflow - Is the quadrilateralized spherical cube map projection the same as Snyder's cubic equal area map projection?](https://gis.stackexchange.com/questions/40957/is-the-quadrilateralized-spherical-cube-map-projection-the-same-as-snyders-cubi) (with comments from Ken Chan)


[![Documentation Status](https://readthedocs.org/projects/qlsc/badge/?version=latest)](https://qlsc.readthedocs.io/en/latest/?badge=latest)

