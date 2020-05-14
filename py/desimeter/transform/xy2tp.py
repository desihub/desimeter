# -*- coding: utf-8 -*-
"""
This module provides the fundamental coordinate transformations between fiber
positioner (theta, phi) angles and an (x, y) cartesian space.

It is a direct copy of static methods from the postransforms module on SVN:
    /code/focalplane/plate_control/trunk/petal/postransforms.py
    SVN revision r131291

That module generally works with python lists and the math module, rather
numpy arrays. That is done for speed, since numpy carries tremendous overhead
in the more atomic operations done on the instrument.

Since desimeter tends to work with numpy arrays, the usage of these functions
is not expected to be done directly. Rather one will use the thin wrappers
in the pos2ptl module.

Another reason for the wrappers is that the pos2ptl module uses some more
specific nomenclature. I.e. there are different kinds of "tp" (internally-
tracked vs externally-measured) and different kinds of xy (petal vs flat vs
local, etc.)
"""

import math
import sys
epsilon = sys.float_info.epsilon

def tp2xy(tp, r):
    """Converts TP angles into XY cartesian coordinates, where arm lengths
    associated with angles theta and phi are respectively r[1] and r[2].
    INPUTS:  tp ... [theta,phi], unit degrees
              r ... [central arm length, eccentric arm length]
    OUTPUT:  xy ... [x,y]
    """
    t = math.radians(tp[0])
    t_plus_p = t + math.radians(tp[1])
    x = r[0] * math.cos(t) + r[1] * math.cos(t_plus_p)
    y = r[0] * math.sin(t) + r[1] * math.sin(t_plus_p)
    return x, y

def xy2tp(xy, r, ranges):
    """Converts XY cartesian coordinates into TP angles, where arm lengths
     associated with angles theta and phi are respectively r[1] and r[2].

    INPUTS:   xy ... [x,y]
               r ... [central arm length, eccentric arm length]
          ranges ... [[min(theta), max(theta)], [min(phi), max(phi)]]

    OUTPUTS:  tp ... [theta,phi], unit degrees
     unreachable ... boolean, True if the requested xy cannot be reached
                     by any tp

    In cases where unreachable == True, the returned tp value will be a
    closest possible approach to the unreachable point requested at xy.
    """
    theta_centralizing_err_tol = 1e-4 # within this much xy error allowance, adjust theta toward center of its range
    n_theta_centralizing_iters = 3 # number of points to try when attempting to centralize theta
    numeric_contraction = epsilon*10 # slight contraction to avoid numeric divide-by-zero type of errors
    x, y, r1, r2 = xy[0], xy[1], r[0], r[1]
    unreachable = False
    
    # adjust targets within reachable annulus
    hypot = (x**2.0 + y**2.0)**0.5
    angle = math.atan2(y, x)
    outer = r[0] + r[1]
    inner = abs(r[0] - r[1])
    if hypot > outer or hypot < inner:
        unreachable = True
    inner += numeric_contraction
    outer -= numeric_contraction
    HYPOT = hypot
    if hypot >= outer:
        HYPOT = outer
    elif hypot <= inner:
        HYPOT = inner
    X = HYPOT*math.cos(angle)
    Y = HYPOT*math.sin(angle)
    
    # transform from cartesian XY to angles TP
    arccos_arg = (X**2.0 + Y**2.0 - (r1**2.0 + r2**2.0)) / (2.0 * r1 * r2)
    arccos_arg = max(arccos_arg, -1.0) # deal with slight numeric errors where arccos_arg comes back like -1.0000000000000002
    arccos_arg = min(arccos_arg, +1.0) # deal with slight numeric errors where arccos_arg comes back like +1.0000000000000002
    P = math.acos(arccos_arg)
    T = angle - math.atan2(r2*math.sin(P), r1 + r2*math.cos(P))
    TP = [math.degrees(T), math.degrees(P)]
    
    # wrap angles into travel ranges
    for i in [0, 1]:
        range_min, range_max = min(ranges[i]), max(ranges[i])
        if TP[i] < range_min: # try +360 phase wrap
            TP[i] += math.floor((range_max - TP[i])/360.0)*360.0
            if TP[i] < range_min:
                TP[i] = range_min
                unreachable = True
        elif TP[i] > range_max: # try -360 phase wrap
            TP[i] -= math.floor((TP[i] - range_min)/360.0)*360.0
            if TP[i] > range_max:
                TP[i] = range_max
                unreachable = True
                
    # centralize theta
    T_ctr = (ranges[0][0] + ranges[0][1])/2.0
    T_options = linspace(TP[0], T_ctr, n_theta_centralizing_iters)
    for T_try in T_options:
        xy_try = tp2xy([T_try, TP[1]], r)
        x_err = xy_try[0] - X
        y_err = xy_try[1] - Y
        vector_err = (x_err**2.0 + y_err**2.0)**0.5
        if vector_err <= theta_centralizing_err_tol:
            TP[0] = T_try
            break

    return tuple(TP), unreachable

def linspace(start,stop,num):
    """Return a list of floats linearly spaced from start to stop (inclusive).
    List has num elements."""
    return [i*(stop-start)/(num-1)+start for i in range(num)]
