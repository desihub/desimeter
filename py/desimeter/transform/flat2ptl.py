# -*- coding: utf-8 -*-
"""
This module provides coordinate transformation between petal cartesian (x, y)
and positioner local pseudo-cartesian "flat" (x, y) coordinates.

    x_ptl, y_ptl   ... Units mm. Cartesian coordinates, relative to petal.
                       Corresponds to ptlXY in postransforms.py.
                     
    x_flat, y_flat ... Units mm. Pseudo-cartesian, relative to petal.
                       Corresponds to flatXY in postransforms.py.
                       
"Flat" cooridnates are produced by the approximation (Q, S) ~ (Q, R). In other
words:
    
    1. (x_ptl, y_ptl) --> normal polar coordinates, with R = radial position.
    2. Approximate that polar R same as surface coordinate S along asphere.
    3. (Q, S) --> (x_flat, y_flat)
    
The spatial warping in flat xy is < 1 um local to a given positioner. In the
petal or global coordinates, the warping causes up to ~ 500 um variation from
physical reality.

On the DESI instrument, the x and y center position of each robot (OFFSET_X
and OFFSET_Y) are calculated in flat xy space.
"""

import numpy as np
import xy2qs

def flat2ptl(x_flat, y_flat):
    '''Converts (x_flat, y_flat) coordinates to (x_ptl, y_ptl).
    See module notes for description of these coordinate systems.
    Input args may be scalars or 1xN lists, tuples, or numpy arrays.
    Output is a tuple of 2 numpy arrays.
    '''
    Q_rad = np.arctan2(y_flat, x_flat)
    S = np.hypot(flat_x, flat_y)
    
    
    R = pc.S2R_lookup(S)
        ptlX = R * math.cos(Q_rad)
        ptlY = R * math.sin(Q_rad)
        return (ptlX, ptlY)
    
    return x_ptl, y_ptl

def ptl2flat(x_ptl, y_ptl):
    '''Converts (x_ptl, y_ptl) coordinates to (x_ptl, y_ptl).
    See module notes for description of these coordinate systems.
    Input args may be scalars or 1xN lists, tuples, or numpy arrays.
    Output is a tuple of 2 numpy arrays.
    '''
    x_ptl = _to_numpy(x_ptl)
    y_ptl = _to_numpy(y_ptl)

       Q_rad = math.atan2(ptlXY[1], ptlXY[0])
        R = math.hypot(ptlXY[0], ptlXY[1])
        S = pc.R2S_lookup(R)
        flatX = S * math.cos(Q_rad)
        flatY = S * math.sin(Q_rad)
        return (flatX, flatY)
    
    return x_ptl, y_ptl

def _to_numpy(u):
    '''Reduces overhead if u is already numpy.'''
    return u if isinstance(u, np.ndarray) else np.array(u)
