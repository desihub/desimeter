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

def xy2tp(xy, r, ranges, t_guess=None, t_guess_tol=20.0):
    """Converts XY cartesian coordinates into TP angles, where arm lengths
     associated with angles theta and phi are respectively r[1] and r[2].

    INPUTS:   xy ... [x,y]
               r ... [central arm length, eccentric arm length]
          ranges ... [[min(theta), max(theta)], [min(phi), max(phi)]]
         t_guess ... optional guess at approx theta expected, unit degrees
     t_guess_tol ... default=20, unit degrees

    OUTPUTS:  tp ... [theta,phi], unit degrees
     unreachable ... boolean, True if the requested xy cannot be reached
                     by any tp

    In cases where unreachable == True, the returned tp value will be a
    closest possible approach to the unreachable point requested at xy.
    
    Guess value for theta (when provided by caller) will be used to disambiguate
    between mathematically possible alternate (t,p) pairs. If not provided, then
    the function will pick a pair for which p is within [0,180].
    """
    numeric_contraction = epsilon*10 # slight contraction to avoid numeric divide-by-zero type of errors
    x, y, r1, r2 = xy[0], xy[1], r[0], r[1]
    unreachable = False

    # adjust targets within reachable annulus
    hypot = (x**2.0 + y**2.0)**0.5
    angle = math.atan2(y, x)
    outer = r1 + r2
    inner = abs(r1 - r2)
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
    
    # check alternate configuration: reflecting the arms across vector (X,Y)
    if t_guess != None:
        tx = r1*math.cos(T)  # i.e., vector to phi fulcrum, which will be reflected
        ty = r1*math.sin(T)
        m = Y/X
        tx_alt = ((1 - m**2)*tx + 2*m*ty)/(m**2 + 1)
        ty_alt = ((m**2 - 1)*ty + 2*m*tx)/(m**2 + 1)
        px_alt = X - tx_alt  # i.e., vector from alternate phi fulcrum to desired X,Y
        py_alt = Y - ty_alt
        T_alt = math.atan2(ty_alt, tx_alt)
        closeness = abs(T - t_guess)
        closeness_alt = abs(T_alt - t_guess)
        if closeness_alt < closeness and closeness_alt <= t_guess_tol:
            T = T_alt
            P = math.atan2(py_alt, px_alt) - T_alt
    
    TP = [math.degrees(T), math.degrees(P)]

    # wrap angles into travel ranges
    for i in [0, 1]:
        range_min, range_max = min(ranges[i]), max(ranges[i])
        TP[i], fail = _wrap_into_range(TP[i], range_min, range_max)
        if fail:
            unreachable = True

    return tuple(TP), unreachable

def _wrap_into_range(angle, range_min, range_max):
    '''Check +/-360 deg phase wraps (as necessary) to put angle within the
    argued range. All units in deg.
    
    INPUT:  angle ... value to be checked
            range_min ... lowest allowed
            range_max ... highest allowed
    
    OUTPUT: new ... new angle
            unreachable ... boolean, True if not not possible to put angle in range
    '''
    new = angle
    unreachable = False
    if new < range_min: # try +360 phase wrap
        new += math.floor((range_max - new)/360.0)*360.0
        if new < range_min:
            new = range_min
            unreachable = True
    elif new > range_max: # try -360 phase wrap
        new -= math.floor((new - range_min)/360.0)*360.0
        if new > range_max:
            new = range_max
            unreachable = True
    return new, unreachable
