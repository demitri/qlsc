
# The version lives in its own module so that setup.py can read it without importing this
# package, which would pull in NumPy and the compiled extension before they are built.
# Ref: https://stackoverflow.com/a/16084844/2712652
from .version import __version__
#__author__ = "Demitri Muna"

from .qlsc_base import QLSC, QLSCIndex
from .qlsc_functions import distance
from .qlsc_functions import sindist
from .qlsc_functions import xy2facenum
from .qlsc_crossmatch import crossmatch
