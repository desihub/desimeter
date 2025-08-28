# -*- coding: utf-8 -*-
"""
Conversions between nominal values for R, S, Z, and N on a nominal Echo22
focal surface. Performed by interpolation of a lookup table. That table was
extracted from DESI-0530.

    R ... (mm) radial distance from optical (Z) axis

    S ... (mm) distance along focal surface from optical axis, in a plane that
          intersects the Z axis

    Z ... (mm) parallel to optical axis, wtih + direction pointing toward primary
          mirror, and Z==0 at the peak of the asphere

    N ... (deg) nutation angle, between optical axis and the central (theta)
          axis of a fiber positioner at that location
"""

import os
import numpy as np
from astropy.table import Table
from importlib.resources import files

data_directory = str(files('desimeter').joinpath('data'))
lookup_file = 'focal_surface_lookup.csv'
lookup_path = os.path.join(data_directory, lookup_file)
table = Table.read(lookup_path, format='csv', comment='#')

R = np.array(table['R'], dtype=np.float64)
Z = np.array(table['Z'], dtype=np.float64)
S = np.array(table['S'], dtype=np.float64)
N = np.array(table['N'], dtype=np.float64)

oob = float('nan') # out of bounds value

def r2s(r):
    '''Convert r (mm) to s (mm) on a nominal Echo22 focal surface.'''
    return np.interp(r, R, S, left=oob, right=oob)

def s2r(s):
    '''Convert s (mm) to r (mm) on a nominal Echo22 focal surface.'''
    return np.interp(s, S, R, left=oob, right=oob)

def r2z(r):
    '''Convert r (mm) to z (mm) on a nominal Echo22 focal surface.'''
    return np.interp(r, R, Z, left=oob, right=oob)

def z2r(z):
    '''Convert z (mm) to r (mm) on a nominal Echo22 focal surface.'''
    return np.interp(z, Z, R, left=oob, right=oob)

def r2n(r):
    '''Convert r (mm) to n (deg) on a nominal Echo22 focal surface.'''
    return np.interp(r, R, N, left=oob, right=oob)

def n2r(n):
    '''Convert n (deg) to r (mm) on a nominal Echo22 focal surface.'''
    return np.interp(n, N, R, left=oob, right=oob)

def s2n(s):
    '''Convert s (mm) to n (deg) on a nominal Echo22 focal surface.
    Composite of s2r --> r2n.'''
    return r2n(s2r(s))  # returns nutation angles in degrees

def s2z(s):
    '''Convert s (mm) to z (mm) on a nominal Echo22 focal surface.
    Composite of s2r --> r2z.'''
    return r2z(s2r(s))

def z2s(z):
    '''Convert z (mm) to s (mm) on a nominal Echo22 focal surface.
    Composite of z2r --> r2s.'''
    return r2s(z2r(z))

def n2s(n):
    '''Convert n (deg) to s (mm) on a nominal Echo22 focal surface.
    Composite of n2r --> r2s.'''
    return r2s(n2r(n))  # takes nutation angles in degrees
