"""
tan2fp: transformations between tangent plane and focal plane xy [mm]

given spherical coordinate angles theta (polar), phi (azimuthal)::

    xtan = sin(theta)*cos(phi)
    ytan = sin(theta)*sin(phi)

phi=0 is aligned with +xtan = +RA = -xfp in DESI
"""

#from .echo22 import tan2fp, fp2tan
from .raytracefit import tan2fp, fp2tan
