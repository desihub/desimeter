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
        t_min, t_max ... range limits on t_int
        
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
        t_ext, p_ext ... externally-measured (theta,phi), see module notes
        r1, r2       ... kinematic arms LENGTH_R1 and LENGTH_R2, scalar or vector
        
    OUTPUTS:
        x_loc, y_loc ... positioner local (x,y), as described in module notes
    '''
    pass

def ext2int(u_ext):
    '''Converts t_ext or p_ext coordinate to t_int or p_int.'''
    pass

def int2ext(u_int):
    '''Converts t_int or p_int coordinate to t_ext or p_ext.'''
    pass

# user gives direction (+1, -1) or zero or none if not known
# don't do wraptol here 
def delta_angle(u0_int, u1_int, u_scale=1.0, wraptol=0):
    '''Special function for the subtraction operation on angles, like:
        
        output = (t1_int - t0_int) / scale
        output = (p1_int - p0_int) / scale
        
    For angles crossing the +/-180 deg boundary, the argument wraptol controls
    whether the delta goes the "short way" or the "long way" around. E.g.:
    
        u0 = +179, u1 = -179, wraptol = 0 --> delta_u = -358
        u0 = +179, u1 = -179, wraptol = 3 --> delta_u = +2
    
    This function also provides consistent application of a scale factor,
    intended to correspond to SCALE_T and SCALE_P, a.k.a. GEAR_CALIB_T and
    GEAR_CALIB_P. Understand that on the instrument, these calibrations are
    applied when converting between number of intended motor rotations and
    corresponding number of expected output shaft rotations. Therefore:
        1. They don't work in absolute coordinates. Only relative operations.
        2. The postransforms module knows nothing about them.
    '''
    pass
        
def addto_angle(t_int, dt_int, t_scale=1.0, t_min=nom_t_min, t_max=nom_t_max):
    '''Special function for the addition operation on internal theta:
        
        output = t_int + dt_int * scale
        
    '''
    # note determine order of wrapping operations based on sign of du
    pass
    
def _to_numpy(u):
    '''Internal function to cast values to vectors.'''
    return u if isinstance(u, np.ndarray) else np.array(u)

def _wrap_theta(tp0, dtdp, posintT_range):
    """
    tp0             : initial TP positions
    dtdp            : delta TP movement
    posintT_range   : allowed range of internal theta, tuple or list

    Returns a modified dtdp after appropriately wrapping the delta theta
    to not cross a physical hardstop. Phi angle is untouched.
    """
    ti, dt, dp = tp0[0], dtdp[0], dtdp[1]  # t initial, dt, dp
    wrapped_dt = dt - 360*pc.sign(dt)
    tf = ti + dt
    wrapped_tf = ti + wrapped_dt
    if min(posintT_range) <= wrapped_tf <= max(posintT_range):
        if ((min(posintT_range) > tf or max(posintT_range) < tf)
                or abs(wrapped_dt) < abs(dt)):
            dt = wrapped_dt
    return dt, dp

if __name__ == '__main__':
    '''Sample code to generate test values below:
    
    import postransforms
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
    x_ptl = xy_ptl[0]
    y_ptl = xy_ptl[1]
    x_flat = xy_flat[0]
    y_flat = xy_flat[1]
    
    Expected values generated using:
        /focalplane/plate_control/trunk/petal/postransforms.py
        SVN revision r131291
    '''
    x_ptl = [   0.        ,  420.        ,    0.        , -420.        ,
          0.        ,   49.36016763,  -44.7721312 ,  290.27533596,
       -171.04187504,  149.13131749,   45.22753108, -156.05544102,
        -72.75112396,  -46.182099  ,  177.77252752, -167.04449541,
        242.58666348,   85.49150994, -170.04762666, -198.57262048,
       -223.6451169 , -279.76966979, -134.83332098, -168.87426158,
        182.33588638,   51.40915708, -272.7590442 ,   -0.02863735,
        217.32671465,  167.88475044,  -55.22431018,  130.04988557,
       -186.1852126 , -238.85164632,  288.65000005, -178.90301878,
        240.75019191,  190.34370977,  124.97949325,  209.45026321,
         18.34138867,  150.48441556,  273.58850183,  -26.27090272,
        177.08586715, -257.11588069,  185.39875346,  -42.57734023,
        268.43058242,  195.78150239, -252.75352011, -104.03205623,
          0.77304182,  190.10570586,  -33.63950671]
    
    y_ptl = [   0.        ,    0.        ,  420.        ,    0.        ,
       -420.        ,  141.59941388, -169.51784977,  162.5367851 ,
       -296.4173894 ,   76.33197112,  179.75760781, -248.94473711,
       -232.98087201,  -51.83902311,  122.03529375,   56.76862158,
       -182.48116967,  192.11681348, -103.99198405,  -42.07583055,
        -84.57594401,  -31.69948182,  -85.57340799,  -45.95503519,
        265.30396862,   64.64055929,   89.26975874,  -97.38735631,
        159.23491464, -161.47531464,   89.94517025, -253.01503035,
        236.77534191, -246.95572445, -145.41599946,   58.17770752,
        101.12768602, -250.13650192,  -70.5522265 ,  265.75574848,
        -88.99420182,  282.14807857, -295.06497539,  282.71937891,
        192.47586294,  177.6687368 , -236.94156772,   62.30608525,
        288.77337204, -132.60177287,   70.61486625,  204.5543887 ,
        263.0779443 ,    1.50756598, -170.65328027]
    
    x_flat = [   0.        ,  420.645406  ,    0.        , -420.645406  ,
          0.        ,   49.36856083,  -44.78281368,  290.54160796,
       -171.20827052,  149.16355712,   45.23970442, -156.16644687,
        -72.78639548,  -46.18364487,  177.83880399, -167.08489391,
        242.77110081,   85.52171316, -170.10117328, -198.63764471,
       -223.74895298, -279.95206958, -134.85959493, -168.9143967 ,
        182.49219419,   51.4116223 , -272.94397312,   -0.02863929,
        217.4561269 ,  167.95851428,  -55.22873845,  130.13648573,
       -186.32454269, -239.08591585,  288.89948529, -178.95273908,
        240.88457481,  190.49909697,  124.99882336,  209.64926941,
         18.34246188,  150.61165923,  273.96847279,  -26.28832794,
        177.18504338, -257.32333766,  185.53716173,  -42.57903321,
        268.78620788,  195.87029762, -252.89607019, -104.07634464,
          0.77348003,  190.15975301,  -33.64739524]
    
    y_flat = [   0.        ,    0.        ,  420.645406  ,    0.        ,
       -420.645406  ,  141.62349145, -169.55829617,  162.68588145,
       -296.70575454,   76.34847278,  179.80599096, -249.12181728,
       -233.09382684,  -51.84075834,  122.08079047,   56.78235066,
       -182.61990912,  192.18468624, -104.02473028,  -42.08960863,
        -84.61521173,  -31.72014875,  -85.59008304,  -45.96595699,
        265.531401  ,   64.64365899,   89.33028308,  -97.39394704,
        159.32973478, -161.54626236,   89.95238268, -253.18351294,
        236.95253068, -247.19794259, -145.54168505,   58.19387613,
        101.18413388, -250.340701  ,  -70.56313855,  266.00825253,
        -88.99940914,  282.38665183, -295.47477376,  282.90690378,
        192.58365827,  177.81209091, -237.11845496,   62.3085627 ,
        289.15594828, -132.66191341,   70.65469223,  204.64147138,
        263.22707473,    1.50799459, -170.69329879]
    
    n = len(x_ptl)
    print(f'Test of {n} uniform random points:')
    
    x_flat_test, y_flat_test = ptl2flat(x_ptl, y_ptl)
    flat_error = np.array([x_flat_test - x_flat, y_flat_test - y_flat])
    flat_error_combined = np.hypot(flat_error[0], flat_error[1])
    print('ptl2flat')
    print(f' err_max = {max(flat_error_combined):.4f}')
    print(f' err_norm = {np.linalg.norm(flat_error):.4f}')
    
    print('')
    x_ptl_test, y_ptl_test = flat2ptl(x_flat, y_flat)
    ptl_error = np.array([x_ptl_test - x_ptl, y_ptl_test - y_ptl])
    ptl_error_combined = np.hypot(ptl_error[0], ptl_error[1])
    print('flat2ptl')
    print(f' err_max = {max(ptl_error_combined):.4f}')
    print(f' err_norm = {np.linalg.norm(ptl_error):.4f}')
    
