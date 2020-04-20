# -*- coding: utf-8 -*-
"""
This module provides coordinate transformation between petal cartesian (x, y)
and several fiber positioner systems.

    x_ptl, y_ptl   ... Units mm. Cartesian coordinates, relative to petal.
                       Corresponds to ptlXY in postransforms.py.
                     
    x_flat, y_flat ... Units mm. Pseudo-cartesian, relative to petal.
                       Corresponds to flatXY in postransforms.py.
    
    x_loc, y_loc   ... Units mm. Pseudo-cartesian, local to a fiber positioner.
                       On local tangent plane to asphere, to a good approximation.
                       Corresponds to poslocXY in postransforms.py.
                       
    t_ext, p_ext   ... Units deg. Externally-measured angles of positioner's
                       kinematic "arms", in the "flat" or "loc" systems.
                       Corresponds to poslocTP in postransforms.py.
                       
    t_int, p_int   ... Units deg. Internally-tracked angles of positioner's
                       kinematic "arms", in a system where theta (and therefore
                       also phi) depend on the angle at which the positioner
                       was physically mounted in the petal.
                       Corresponds to posintTP in postransforms.py.
                       Also corresponds to POS_T, POS_P in the online database.
                       

TYPICAL INPUTS / OUTPUTS:
-------------------------
    Inputs: Args may be scalars or 1xN lists, tuples, or numpy arrays.
    Output: Tuple of 2 numpy arrays.

"FLAT" SPACE
------------
"Flat" cooridnates are produced by the approximation (Q, S) ~ (Q, R). In other
words:
    
    1. (x_ptl, y_ptl) --> normal polar coordinates, with R = radial position.
    2. Approximate that polar R same as surface coordinate S along asphere.
    3. (Q, S) --> (x_flat, y_flat)
    
The "loc" coordinates are a translation within "flat".
    
The spatial warping in flat xy is < 1 um local to a given positioner. In the
petal or global coordinates, the warping causes up to ~ 500 um variation from
physical reality.

CALIBRATION PARAMETERS:
-----------------------
On the DESI instrument, the calibration parameters (OFFSET_X, OFFSET_Y,
LENGTH_R1, LENGTH_R2, OFFSET_T, OFFSET_P) are all calculated in this pseudo-
cartesian space.

ANGULAR RANGE LIMITS:
---------------------
The physical fiber positioner robot has hard limits on its angular travel
ranges. In usage on the instrument, they are applied to t_int and p_int. The
limit values are derived from a number calibration parameters: PHYSICAL_RANGE_T,
PHYSICAL_RANGE_P, PRINCIPLE_HARDSTOP_CLEARANCE_T, PRINCIPLE_HARDSTOP_CLEARANCE_P,
SECONDARY_HARDSTOP_CLEARANCE_T, SECONDARY_HARDSTOP_CLEARANCE_P, BACKLASH, and
near_full_range_reduced_hardstop_clearance_factor. C.f.
https://desi.lbl.gov/svn/code/focalplane/plate_control/trunk/petal/posmodel.py
https://desi.lbl.gov/svn/code/focalplane/plate_control/trunk/petal/posconstants.py

The use of these parameters to determine range limits varies depending on the
operational context.

This module makes the greatly simplifying assumption that it will not be used
for instrument control, only for target selection or analyis purposes. Therefore
angular range limits are NOT implemented here.
"""

import numpy as np
import xy2qs
import xy2tp

_default_t_int_range = [-180. + xy2tp.epsilon, 180.] # theta min, max
_default_p_int_range = [-20., 200.] # phi min, max

def ptl2flat(x_ptl, y_ptl):
    '''Converts (x_ptl, y_ptl) coordinates to (x_ptl, y_ptl).'''
    q = np.arctan2(y_ptl, x_ptl)
    r = np.hypot(x_ptl, y_ptl)
    s = xy2qs.r2s(r)
    x_flat = s * np.cos(q)
    y_flat = s * np.sin(q)
    return x_flat, y_flat

def flat2ptl(x_flat, y_flat):
    '''Converts (x_flat, y_flat) coordinates to (x_ptl, y_ptl).'''
    q = np.arctan2(y_flat, x_flat)
    s = np.hypot(x_flat, y_flat)
    r = xy2qs.s2r(s)
    x_ptl = r * np.cos(q)
    y_ptl = r * np.sin(q)
    return x_ptl, y_ptl

def flat2loc(u_flat, u_offset):
    '''Converts x_flat or y_flat coordinate to x_loc or y_loc.
    READ THE FINE PRINT: Args here are not  like "(x,y)". More like "(x,xo)".
        
    INPUTS
        u_flat   ... x_flat or y_flat
        u_offset ... OFFSET_X or OFFSET_Y
        
    OUTPUTS:
        u_loc    ... x_loc or y_loc
        
    u_offset may either be a scalar (will be applied to all points), or
    a vector of unique values per point.
    '''
    u_offset = _to_numpy(u_offset) # so that minus sign works in next line
    u_loc = _add_offset(u_flat, -u_offset)
    return u_loc
    
def loc2flat(u_loc, u_offset):
    '''Converts x_loc or y_loc coordinate to x_flat or y_flat.
    READ THE FINE PRINT: Args here are not  like "(x,y)". More like "(x,xo)".
    
    INPUTS:
        u_loc    ... x_loc or y_loc
        u_offset ... OFFSET_X or OFFSET_Y
        
    OUTPUTS:
        u_flat   ... x_flat or y_flat
    
    u_offset may either be a scalar (will be applied to all points), or
    a vector of unique values per point.
    '''
    u_flat = _add_offset(u_loc, u_offset)
    return u_flat

def loc2ext(x_loc, y_loc, r1, r2, t_offset, p_offset):
    '''Converts (x_loc, y_loc) coordinates to (t_ext, p_ext).
    READ THE FINE PRINT: Returned tuple has 3 elements.
    
    INPUTS:
        x_loc, y_loc ... positioner local (x,y), as described in module notes
        r1, r2       ... kinematic arms LENGTH_R1 and LENGTH_R2, scalar or vector
        t_offset, p_offset ... [deg] OFFSET_T and OFFSET_P
        
    OUTPUTS:
        t_ext, p_ext ... externally-measured (theta,phi), see module notes
        unreachable  ... vector of same length, containing booleans
        
    The vector "unreachable" identifies any points where the desired (x,y) can
    not be reached by the positioner. In such cases, the returned t_ext, p_ext
    correspond to a point which is as near as possible to x_loc, y_loc. This
    matches how it works on the instrument.
    
    (The purpose of returning these "nearest possible" points, rather than a
    simple "fail", is to be operationally tolerant of targets that may be just
    a few um from the edge of reachability. In these cases, near-enough can
    still deliver good spectra.)
    '''
    x_loc = _to_list(x_loc)
    y_loc = _to_list(y_loc)
    r1 = _to_list(r1)
    r2 = _to_list(r2)
    n = len(x_loc)
    if len(r1) != n:
        r1 = [r1[0]]*n
    if len(r2) != n:
        r2 = [r2[0]]*n
    t_ext = []
    p_ext = []
    unreachable = []
    for i in range(n):
        ext_ranges = [[float(int2ext(_default_t_int_range[0], t_offset[i])),
                       float(int2ext(_default_t_int_range[1], t_offset[i]))],
                      [float(int2ext(_default_p_int_range[0], p_offset[i])),
                       float(int2ext(_default_p_int_range[1], p_offset[i]))]]
        tp_ext, unreach = xy2tp.xy2tp(xy=[x_loc[i], y_loc[i]],
                                     r=[r1[i], r2[i]],
                                     ranges=ext_ranges)
        t_ext.append(tp_ext[0])
        p_ext.append(tp_ext[1])
        unreachable.append(unreach)
    t_ext = _to_numpy(t_ext)
    p_ext = _to_numpy(p_ext)
    unreachable = _to_numpy(unreachable)
    return t_ext, p_ext, unreachable

def ext2loc(t_ext, p_ext, r1, r2):
    '''Converts (t_ext, p_ext) coordinates to (x_loc, y_loc).
    
    INPUTS:
        t_ext, p_ext ... [deg] externally-measured (theta,phi), see module notes
        r1, r2       ... kinematic arms LENGTH_R1 and LENGTH_R2, scalar or vector
        
    OUTPUTS:
        x_loc, y_loc ... positioner local (x,y), as described in module notes
    '''
    t = np.radians(t_ext)
    t_plus_p = t + np.radians(p_ext)
    r1 = _to_numpy(r1)
    r2 = _to_numpy(r2)
    x_loc = r1 * np.cos(t) + r2 * np.cos(t_plus_p)
    y_loc = r1 * np.sin(t) + r2 * np.sin(t_plus_p)
    return x_loc, y_loc

def ext2int(u_ext, u_offset):
    '''Converts t_ext or p_ext coordinate to t_int or p_int.
    READ THE FINE PRINT: Args here are not  like "(t,p)". More like "(t,to)".
    
    INPUTS:
        u_ext    ... t_ext or p_ext
        u_offset ... OFFSET_T or OFFSET_P
        
    OUTPUTS:
        u_int   ... t_int or p_int
    
    u_offset may either be a scalar (will be applied to all points), or
    a vector of unique values per point.
    '''
    u_offset = _to_numpy(u_offset) # so that minus sign works in next line
    u_int = _add_offset(u_ext, -u_offset)
    return u_int

def int2ext(u_int, u_offset):
    '''Converts t_int or p_int coordinate to t_ext or p_ext.
    READ THE FINE PRINT: Args here are not  like "(t,p)". More like "(t,to)".
        
    INPUTS:
        u_ext    ... t_ext or p_ext
        u_offset ... OFFSET_T or OFFSET_P
        
    OUTPUTS:
        u_int   ... t_int or p_int
    
    u_offset may either be a scalar (will be applied to all points), or
    a vector of unique values per point.
    '''
    u_ext = _add_offset(u_int, u_offset)
    return u_ext

def delta_angle(u_start, u_final, direction=0):
    '''Special function for distance between two angle positions. Allows caller
    to select the direction of rotation. This is important when modeling
    behavior of the instrument, because real positioners have theta travel
    ranges that go a bit beyond +/-180 deg. So specifying "direction" (when
    known) can cover special cases like (b) below:
    
        a. u_start = +179, u_final = -179, direction = 0 or -1 --> delta_u = -358
        b. u_final = +179, u_final = -179, direction = +1      --> delta_u = +2
        
    INPUTS:
        u_start   ... units deg, start t_int or p_int
        u_final   ... units deg, final t_int or p_int
        direction ... direction delta should go around the circle
        
    All of the above inputs can be either a scalar or a vector.        
    
    The direction flags work like this:
        
         1 --> counter-clockwise --> output du >= 0, and abs(du)
        -1 --> clockwise         --> output du <= 0, and abs(du)
         0 --> unspecified       --> sign(du) = sign(u1-u0)
        
    For the special case where abs(du) > 360 deg, while at the same time
    sign(u_final - u_start) != direction, then some decision must be made
    about how to handle the number of revolutions about the circle. In this
    case, the opposite-direction output will be kept within magnitude 360 deg.
    For example,
    
        u_start = 0, u_final = 400, direction = 0 or +1 --> delta_u = +400
        u_start = 0, u_final = 400, direction = -1      --> delta_u = -320
        u_start = 0, u_final = 800, direction = 0 or +1 --> delta_u = +800
        u_start = 0, u_final = 800, direction = -1      --> delta_u = -280
        
    OUTPUTS:
        delta_u ... angular distance (signed)
    '''
    u_start = _to_numpy(u_start)
    u_final = _to_numpy(u_final)
    direction = np.ones(np.shape(u_start)) * _to_numpy(direction)
    natural_delta = u_final - u_start
    forced_terms = np.abs(natural_delta) * np.sign(direction)
    natural_terms = natural_delta * (direction == 0)
    delta_u = forced_terms + natural_terms
    return delta_u
    
def _to_numpy(u):
    '''Internal function to cast values to consistent numpy vectors.'''
    if isinstance(u, list):
        return np.array(u)
    if np.size(u) == 1:
        return np.array([float(u)])
    if isinstance(u, np.ndarray):
        return u
    return np.array(u)

def _to_list(u):
    '''Internal function to cast values to consistent python list vectors.'''
    if isinstance(u, list):
        return u
    if np.size(u) == 1:
        return [float(u)]
    if isinstance(u, np.ndarray):
        return u.tolist()
    return np.array(u).tolist()

def _add_offset(u, offset):
    '''Internal function to apply arithmetic offsets consistently to vectors.'''
    u = _to_numpy(u)
    offset = _to_numpy(offset)
    return u + offset

if __name__ == '__main__':
    '''Tests below compare results from this module to a presumed canonical
    set of coordinate conversions using postransforms.py. See utility script
    desimeter/bin/generate_pos2ptl_testdata, for code to generate the data set.
    '''
    import json
    from pkg_resources import resource_filename
    import os
    
    desimeter_data_dir = resource_filename("desimeter", "data")
    filename = 'pos2ptl_test_data.json'
    path = os.path.join(desimeter_data_dir, filename)
    with open(path, 'r') as file:
        u = json.load(file)
                                                                             
    tests = {'ptl2flat': {'func': ptl2flat, 'inputs': ['x_ptl', 'y_ptl'], 'checks': ['x_flat', 'y_flat']},
             'flat2ptl': {'func': flat2ptl, 'inputs': ['x_flat', 'y_flat'], 'checks': ['x_ptl', 'y_ptl']},
             'flat2loc_x': {'func': flat2loc, 'inputs': ['x_flat', 'OFFSET_X'], 'checks': ['x_loc']},
             'flat2loc_y': {'func': flat2loc, 'inputs': ['y_flat', 'OFFSET_Y'], 'checks': ['y_loc']},
             'loc2flat_x': {'func': loc2flat, 'inputs': ['x_loc', 'OFFSET_X'], 'checks': ['x_flat']},
             'loc2flat_y': {'func': loc2flat, 'inputs': ['y_loc', 'OFFSET_Y'], 'checks': ['y_flat']},
             'ext2loc': {'func': ext2loc, 'inputs': ['t_ext', 'p_ext', 'LENGTH_R1', 'LENGTH_R2'], 'checks': ['x_loc', 'y_loc']},
             'loc2ext': {'func': loc2ext, 'inputs': ['x_loc', 'y_loc', 'LENGTH_R1', 'LENGTH_R2', 'OFFSET_T', 'OFFSET_P'], 'checks': ['t_ext', 'p_ext']},
             'ext2int_t': {'func': ext2int, 'inputs': ['t_ext', 'OFFSET_T'], 'checks': ['t_int']},
             'ext2int_p': {'func': ext2int, 'inputs': ['p_ext', 'OFFSET_P'], 'checks': ['p_int']},
             'int2ext_t': {'func': int2ext, 'inputs': ['t_int', 'OFFSET_T'], 'checks': ['t_ext']},
             'int2ext_p': {'func': int2ext, 'inputs': ['p_int', 'OFFSET_P'], 'checks': ['p_ext']},
             }
    
    # test of delta_angle() function defined separately
    # (there is no exact matching function in postransforms.py)
    delta = {'u_start':   [15, 15,   15, 0,  90,   90,  90,   0,   0,    0,   0,   0,    0,   0,   0,    0,   0,   0,    0,  800,  800,  800],
             'u_final':   [45, 45,   45, 0,  45,   45,  45, 180, 180,  180, 360, 360,  360, 400, 400,  400, 800, 800,  800,    0,    0,    0],
             'direction': [ 0, +1,   -1, 0,   0,   +1,  -1,   0,  +1,   -1,   0,  +1,   -1,   0,  +1,   -1,   0,  +1,   -1,    0,   +1,   -1],
             'delta_u':   [30, 30, -330, 0, -45, +315, -45, 180, 180, -180, 360, 360, -360, 400, 400, -320, 800, 800, -280, -800, +280, -800]}
    u.update(delta)
    delta_test = {'delta_angle': {'func':delta_angle,
                                  'inputs': ['u_start', 'u_final', 'direction'],
                                  'checks':['delta_u']}}
    tests.update(delta_test)
    
    all_out = {}
    for name, args in tests.items():
        inputs = [u[key] for key in args['inputs']]
        checks = [u[key] for key in args['checks']]
        if name == 'ext2loc':
            # only reachable targets are fair tests for this conversion
            reachable = np.array(u['unreachable']) == False
            inputs = [np.array(vec)[reachable].tolist() for vec in inputs]
            checks = [np.array(vec)[reachable].tolist() for vec in checks]
        outputs = args['func'](*inputs)
        if isinstance(outputs, tuple):
            outputs = list(outputs)
        else:
            outputs = [outputs]
        errors = []
        for i in range(len(args['checks'])):
            error = outputs[i] - checks[i]
            errors += [error.tolist()]
        errors = np.array(errors)
        all_out[name] = {'inputs':inputs, 'checks':checks, 'outputs': outputs, 'errors':errors}
        print(f'\n{name} on {len(error)} points:')
        print(f' max(abs(err)) = {np.max(np.abs(errors)):8.6f}')
        print(f'    norm(err)  = {np.linalg.norm(errors):8.6f}')
        worst = np.argmax(np.abs(errors)) % len(errors[0])
        print(f' Worst case at test point {worst}:')
        for i in range(len(errors)):
            print(f'  {args["inputs"][i]}={inputs[i][worst]:9.4f} --> ' +
                  f'{args["checks"][i]}={outputs[i][worst]:9.4f}, ' +
                  f'expected={checks[i][worst]:9.4f} --> ' +
                  f'error={errors[i][worst]:8.2E}')
