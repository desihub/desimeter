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
    u_flat = _to_numpy(u_flat)
    u_offset = _to_numpy(u_offset)
    u_loc = u_flat - u_offset
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
    u_loc = _to_numpy(u_loc)
    u_offset = _to_numpy(u_offset)
    u_flat = u_loc + u_offset
    return u_flat

def loc2ext(x_loc, y_loc, r1, r2):
    '''Converts (x_loc, y_loc) coordinates to (t_ext, p_ext).
    READ THE FINE PRINT: Returned tuple has 3 elements.
    
    INPUTS:
        x_loc, y_loc ... positioner local (x,y), as described in module notes
        r1, r2       ... kinematic arms LENGTH_R1 and LENGTH_R2, scalar or vector
        t_min, t_max ... [deg] range limits on t_int
        
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
    pass

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
    x = r1 * np.cos(t) + r2 * np.cos(t_plus_p)
    y = r1 * np.sin(t) + r2 * np.sin(t_plus_p)
    return x, y

def ext2int(u_ext):
    '''Converts t_ext or p_ext coordinate to t_int or p_int.'''
    pass

def int2ext(u_int):
    '''Converts t_int or p_int coordinate to t_ext or p_ext.'''
    pass

def delta_angle(u0, u1, scale=1.0, direction=None):
    '''Special function for the subtraction operation on angles, like:
        
        output = (t1_int - t0_int) / scale, or
        output = (p1_int - p0_int) / scale
        
    INPUTS:
        u0        ... start t_int or p_int
        
        u1        ... final t_int or p_int
        
        scale     ... SCALE_T or SCALE_P (see notes below)
        
        direction ... specify direction delta should go around the circle
                       1 --> counter-clockwise    --> output du >= 0
                      -1 --> clockwise            --> output du <= 0
                       0 or None --> unspecificed --> sign(du) = sign(u1-u0)
        
        All of the above inputs can be either a scalar or a vector.
        
    OUTPUTS:
        du ... angular distance (signed)

    This function provides consistent application of a factor "scale",
    intended to correspond to SCALE_T and SCALE_P, a.k.a. GEAR_CALIB_T and
    GEAR_CALIB_P.
    
    Note: on the DESI instrument, scale calibrations are applied only when
    converting between number of intended motor rotations and corresponding
    number of expected output shaft rotations. Therefore:
        
        1. They don't work in absolute coordinates. Only relative operations.
        2. The postransforms module knows nothing about them.
        
    When the caller has definite knowledge of the true direction of rotation,
    they may specify value for "direction". This is important when modeling
    behavior of the instrument, because real positioners have theta travel
    ranges that go a bit beyond +/-180 deg. So specifying "direction" (when
    known) can cover special case like (b) below:
    
        a. u0 = +179, u1 = -179, direction = 0 or -1 --> delta_u = -358
        b. u0 = +179, u1 = -179, direction = +1      --> delta_u = +2
    '''
    pass
        
def addto_angle(u0, du, scale=1.0):
    '''Special function for the addition operation on angles, like:
        
        output = t_int + dt_int * scale, or
        output = p_int + dp_int * scale
        
    This is not just "angle1 + angle2", since scale is only applied to the
    "delta" term, du. See notes in delta_angle() regarding scale.
    
    INPUTS:
        u0    ... initial angle
        du    ... delta angle
        scale ... SCALE_T or SCALE_P
    
    OUTPUTS:
        u1    ... result angle
    '''
    pass
    
def _to_numpy(u):
    '''Internal function to cast values to vectors.'''
    return u if isinstance(u, np.ndarray) else np.array(u)

if __name__ == '__main__':
    '''Sample code to generate test values below:
    
    import postransforms # c.f. SVN... /code/focalplane/plate_control/trunk/petal
    import numpy as np
    trans = postransforms.PosTransforms(stateless=True)
    r_max = 420
    minval = -r_max/np.sqrt(2)
    maxval = r_max/np.sqrt(2)
    n_pts = 50
    xy_ptl_fixed = np.array([[0, r_max, 0, -r_max, 0], [0, 0, r_max, 0, -r_max]])
    xy_ptl_rand = (maxval - minval) * np.random.rand(2, n_pts) + minval
    xy_ptl = np.append(xy_ptl_fixed, xy_ptl_rand, axis=1)
    xy_ptl_list = np.transpose(xy_ptl).tolist()
    xy_flat_list = [trans.ptlXY_to_flatXY(xy) for xy in xy_ptl_list]
    xy_flat = np.transpose(xy_flat_list)
    
    u = {'x_ptl': xy_ptl[0],
         'y_ptl': xy_ptl[1],
         'x_flat': xy_flat[0],
         'y_flat': xy_flat[1],
         }
    
    Expected values generated using:
        /code/focalplane/plate_control/trunk/petal/postransforms.py
        SVN revision r131291
    '''
    from numpy import array
    u = {'x_ptl': array([   0.        ,  420.        ,    0.        , -420.        ,
          0.        ,  198.6797643 ,  217.76367415,  -68.96526497,
       -118.77563685,  -11.8485883 ,  -86.0621337 , -145.3024218 ,
       -294.18403665,   65.57687226,  -50.37019779,   73.11852638,
         67.46758988, -182.8670359 , -159.9806036 ,  -26.69621344,
        106.72276593,    9.91677315,  153.84654707,  -66.89887725,
        -14.05876249,   67.44838061,  210.09257336,  248.83002242,
       -122.54501418, -232.75127458,  127.40879077, -270.87779555,
        -60.90895377,   31.09059788, -178.53496939,   67.22125365,
        -60.17374239,  118.6183532 , -230.22221086,  -87.213702  ,
        148.7465539 ,  124.61477444,  233.30921457, -147.97593093,
        -54.21698817,   48.00843646,  -42.71395975,   40.60865468,
        -69.66016231,  -91.81288477,   59.85044594,  213.37756302,
        104.41605221, -119.46251155,  -89.72524086]), 'y_ptl': array([   0.        ,    0.        ,  420.        ,    0.        ,
       -420.        , -207.63569632,   50.97048504, -280.69605571,
       -209.15620472,    9.30864291,  -16.28223159,  183.33541958,
         -8.30041178, -294.79922592, -200.52105041,  269.10930224,
        140.18857537,  132.22992219,  150.31295761, -200.34956021,
       -280.59467187,   89.69159621, -151.04413592,  -74.63020613,
         12.90079222,  279.2330591 ,   52.44596121,  168.50703147,
        103.70798413,   97.82340363,   -0.77061707, -165.38938527,
        -98.84141826,  -88.2088539 ,   48.81445353,   87.97485867,
       -212.43773764, -198.19383034,  109.87821876,  -95.19456588,
        218.12755685, -167.48613175,   50.26264034,  160.95828037,
        -60.24350088, -242.35785372,  119.76292437,  119.51229135,
        239.13288175,   56.17756835, -190.24161258,  241.72042934,
       -249.78610252, -289.93275376, -188.54370551]), 'x_flat': array([ 0.00000000e+00,  4.20645406e+02,  2.57571025e-14, -4.20645406e+02,
        2.57571025e-14,  1.98814837e+02,  2.17851440e+02, -6.90127081e+01,
       -1.18831478e+02, -1.18486067e+01, -8.60667977e+01, -1.45366843e+02,
       -2.94394006e+02,  6.56262105e+01, -5.03873580e+01,  7.31652660e+01,
        6.74800136e+01, -1.82942158e+02, -1.60042568e+02, -2.67048758e+01,
        1.06802088e+02,  9.91734511e+00,  1.53903887e+02, -6.69036868e+01,
       -1.40587969e+01,  6.74941979e+01,  2.10171609e+02,  2.49015361e+02,
       -1.22569168e+02, -2.32872397e+02,  1.27424045e+02, -2.71103344e+02,
       -6.09149325e+01,  3.10925307e+01, -1.78582870e+02,  6.72272188e+01,
       -6.01973397e+01,  1.18669548e+02, -2.30344629e+02, -8.72244436e+01,
        1.48831497e+02,  1.24658119e+02,  2.33417118e+02, -1.48032759e+02,
       -5.42194879e+01,  4.80323144e+01, -4.27190515e+01,  4.06134199e+01,
       -6.96954003e+01, -9.18205586e+01,  5.98693146e+01,  2.13561073e+02,
        1.04478856e+02, -1.19559564e+02, -8.97564643e+01]), 'y_flat': array([ 0.00000000e+00,  0.00000000e+00,  4.20645406e+02,  5.15142050e-14,
       -4.20645406e+02, -2.07776857e+02,  5.09910277e+01, -2.80889154e+02,
       -2.09254537e+02,  9.30865733e+00, -1.62831140e+01,  1.83416703e+02,
       -8.30633606e+00, -2.95021025e+02, -2.00589364e+02,  2.69281325e+02,
        1.40214390e+02,  1.32284242e+02,  1.50371178e+02, -2.00414569e+02,
       -2.80803226e+02,  8.96967693e+01, -1.51100431e+02, -7.46355715e+01,
        1.29008238e+01,  2.79422740e+02,  5.24656911e+01,  1.68632542e+02,
        1.03728425e+02,  9.78743104e+01, -7.70709331e-01, -1.65527098e+02,
       -9.88511204e+01, -8.82143376e+01,  4.88275504e+01,  8.79826655e+01,
       -2.12521046e+02, -1.98279369e+02,  1.09936645e+02, -9.52062904e+01,
        2.18252121e+02, -1.67544388e+02,  5.02858863e+01,  1.61020095e+02,
       -6.02462785e+01, -2.42478395e+02,  1.19777201e+02,  1.19526316e+02,
        2.39253848e+02,  5.61822638e+01, -1.90301589e+02,  2.41928315e+02,
       -2.49936342e+02, -2.90168298e+02, -1.88609317e+02])}
                                                                             
    tests = {ptl2flat: {'inputs': ('x_ptl', 'y_ptl'), 'checks': ('x_flat', 'y_flat')},
             flat2ptl: {'inputs': ('x_flat', 'y_flat'), 'checks': ('x_ptl', 'y_ptl')}}
    
    def calculable(x):
        f = lambda y: not isinstance(y, bool) and np.isfinite(y)
        # try/except here because numpy sometimes gives vectors and sometimes scalars
        try: return f(x[0])
        except: return f(x)    
    for func, args in tests.items():
        inputs = [u[key] for key in args['inputs']]
        checks = [u[key] for key in args['checks']]
        outputs = func(*inputs)
        errors = []
        for i in range(len(outputs)):
            if calculable(outputs[i]):
                error = outputs[i] - checks[i]
                errors += [error.tolist()]
        errors = np.array(errors)
        print(f'\n{func.__name__} on {len(error)} points:')
        print(f' err_max  = {np.max(errors):6.4f}')
        print(f' err_norm = {np.linalg.norm(errors):6.4f}')
        worst = np.argmax(errors) % len(errors[0])
        print(' Worst case:')
        for i in range(len(errors)):
            print(f'  {args["inputs"][i]}={inputs[i][worst]:9.4f} --> ' +
                  f'{args["checks"][i]}={outputs[i][worst]:9.4f}, ' +
                  f'expected={checks[i][worst]:9.4f} --> ' +
                  f'error={errors[i][worst]:6.4f}')
