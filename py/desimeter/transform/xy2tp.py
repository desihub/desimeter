# -*- coding: utf-8 -*-
"""
This module provides the fundamental coordinate transformations between fiber
positioner (theta, phi) angles and an (x, y) cartesian space.

It is kept manually synchronized with a file of the same name in the online
instrument code:
    /code/focalplane/plate_control/<some_branch>/petal/xy2tp.py

This module generally works with python lists and the math module, rather
numpy arrays. That is done for speed, since numpy carries tremendous overhead
in the small, atomic operations done when operating the instrument.
"""

import math
import sys
try:
    import posconstants as pc
    default_t_guess_tol = pc.default_t_guess_tol
except:
    default_t_guess_tol = 30.0
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

def xy2tp(xy, r, ranges, t_guess=None, t_guess_tol=default_t_guess_tol):
    """Converts XY cartesian coordinates into TP angles, where arm lengths
     associated with angles theta and phi are respectively r[1] and r[2].

    INPUTS:   xy ... [x,y]
               r ... [central arm length, eccentric arm length]
          ranges ... [[min(theta), max(theta)], [min(phi), max(phi)]]
         t_guess ... optional guess at approx theta expected, unit degrees
     t_guess_tol ... ignore t_guess if more than this far from mathematically
                     possible theta options, default=20, unit degrees

    OUTPUTS:  tp ... [theta,phi], unit degrees
     unreachable ... boolean, True if the requested xy cannot be reached by any tp

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
    
    # primary configuration of the arms
    TP = [math.degrees(T), math.degrees(P)]
    TP, range_fail = _wrap_TP_into_ranges(TP, ranges)
    
    # test alternate configurations
    if t_guess != None:
        options = [{'TP':TP, 'range_fail':range_fail}]
        
        # reflecting the arms across vector (X,Y)
        tx = r1*math.cos(T)  # i.e., vector to phi fulcrum, which will be reflected
        ty = r1*math.sin(T)
        m = Y/X
        tx_alt = ((1 - m**2)*tx + 2*m*ty)/(m**2 + 1)
        ty_alt = ((m**2 - 1)*ty + 2*m*tx)/(m**2 + 1)
        px_alt = X - tx_alt  # i.e., vector from alternate phi fulcrum to desired X,Y
        py_alt = Y - ty_alt
        T_alt = math.atan2(ty_alt, tx_alt)
        P_alt = math.atan2(py_alt, px_alt) - T_alt
        TP_alt = [math.degrees(T_alt), math.degrees(P_alt)]
        TP_alt, range_fail_alt = _wrap_TP_into_ranges(TP_alt, ranges)
        options.append({'TP':TP_alt, 'range_fail':range_fail_alt})
        
        # rotating theta by +/-360 deg
        for center_TP in [opt['TP'] for opt in options]:
            TP_alts = [[center_TP[0] + wrap, center_TP[1]] for wrap in [-360.0, 360.0]]
            for TP_alt in TP_alts:
                TP_alt, range_fail_alt = _wrap_TP_into_ranges(TP_alt, ranges)
                options.append({'TP':TP_alt, 'range_fail':range_fail_alt})
        
        # selection
        for opt in options:
            opt['closeness'] = abs(opt['TP'][0] - t_guess)
            opt['closeness_fail'] = opt['closeness'] > t_guess_tol
        all_fail = all(opt['closeness_fail'] for opt in options) or \
                   all(opt['range_fail'] for opt in options)
        if not all_fail:
            options.sort(key=lambda x: x['closeness'])  # tertiary sort
            options.sort(key=lambda x: x['closeness_fail'])  # secondary sort
            options.sort(key=lambda x: x['range_fail'])  # primary sort
            TP = options[0]['TP']
            range_fail = options[0]['range_fail']
            
    unreachable |= range_fail
    return tuple(TP), unreachable

def _wrap_TP_into_ranges(tp, ranges):
    '''Call _wrap_into_range for theta and phi angles.
    
    INPUT:  tp ... [theta, phi], units deg
            ranges ... see xy2tp docstr, units deg
            
    OUTPUT: TP ... wrapped angles
            unreachable ... boolean, true if not possible to put either angle in range
    '''
    TP = [None, None]
    unreachable = False
    for i in [0, 1]:
        range_min, range_max = min(ranges[i]), max(ranges[i])
        TP[i], range_fail = _wrap_into_range(tp[i], range_min, range_max)
        unreachable |= range_fail
    return TP, unreachable

def _wrap_into_range(angle, range_min, range_max):
    '''Check +/-360 deg phase wraps (as necessary) to put angle within the
    argued range. All units in radians.
    
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
