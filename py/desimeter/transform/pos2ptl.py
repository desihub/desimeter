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
import tp2xy

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
    rand = lambda: 2 * (np.random.rand() - 0.5)
    params = {'LENGTH_R1': 3 + 0.4 * rand(),
              'LENGTH_R2': 3 + 0.4 * rand(),
              'OFFSET_T': 180 * rand(),
              'OFFSET_P': 20 * rand(),
              'OFFSET_X': r_max * rand(),
              'OFFSET_Y': r_max * rand(),
              }
    trans.alt = params
    minval = -r_max/np.sqrt(2)
    maxval = r_max/np.sqrt(2)
    n_pts = 50
    xy_ptl_fixed = np.array([[0, r_max, 0, -r_max, 0], [0, 0, r_max, 0, -r_max]])
    xy_ptl_rand = (maxval - minval) * np.random.rand(2, n_pts) + minval
    xy_ptl = np.append(xy_ptl_fixed, xy_ptl_rand, axis=1)
    xy_ptl_list = np.transpose(xy_ptl).tolist()
    xy_flat_list = [trans.ptlXY_to_flatXY(xy) for xy in xy_ptl_list]
    xy_flat = np.transpose(xy_flat_list)
    xy_loc_list = [trans.flatXY_to_poslocXY(xy) for xy in xy_flat_list]
    xy_loc = np.transpose(xy_loc_list)
    tp_ext_list = [trans.poslocXY_to_poslocTP(xy)[0] for xy in xy_loc_list]
    tp_ext = np.transpose(tp_ext_list)
    tp_int_list = [trans.poslocTP_to_posintTP(tp) for tp in tp_ext_list]
    tp_int = np.transpose(tp_int_list)
    
    u = {'x_ptl': xy_ptl[0], 'y_ptl': xy_ptl[1],
         'x_flat': xy_flat[0], 'y_flat': xy_flat[1],
         'x_loc': xy_loc[0], 'y_loc': xy_loc[1],
         't_ext': tp_ext[0], 'p_ext': tp_ext[1],
         't_int': tp_int[0], 'p_int': tp_int[1],
         }
    u.update(params)
    
    Expected values generated using:
        /code/focalplane/plate_control/trunk/petal/postransforms.py
        SVN revision r131291
    '''
    from numpy import array
    u = {'x_ptl': array([   0.        ,  420.        ,    0.        , -420.        ,
          0.        ,  -49.71170394, -153.27517051,  -24.83618853,
        -87.80039948,  -86.35515488, -225.63124148,   59.51609504,
         77.69144034, -232.51850268,  -88.04976281,  150.42302658,
        146.42099892,    8.5332472 ,  -64.20358526,  -87.92379287,
        200.11358343,    8.90723718, -226.3316552 ,  229.35821985,
       -256.22681083,  -61.28603271,   64.61589362,  282.90313764,
       -257.1243235 ,   77.7868476 , -237.03783857,  249.87872826,
       -199.44972791, -235.96223061,   -0.88383157, -103.93812474,
        238.23929311,  195.50696481,   -6.72547404, -182.95437097,
        -75.25863082,  -96.23397059, -169.9025175 ,  -78.01053368,
       -210.65020817,  115.07601958, -201.22582179, -176.11442907,
        -76.28475675,   22.17648019,  137.41293666,  -53.69131367,
       -256.92635759,  -75.29288387, -104.42993051]), 'y_ptl': array([   0.        ,    0.        ,  420.        ,    0.        ,
       -420.        , -158.9212518 ,  239.41202143,  216.84380098,
        -95.83182443, -226.79285853,   51.13014003,  -97.70543404,
        -30.44307117, -169.02940538,   40.47949605, -102.67597696,
        109.44054966,    4.33082601, -164.79206683, -238.13076979,
        224.5050337 ,  178.81452729,    7.37963012, -240.48724097,
          3.4766854 , -148.87618999,  234.85454852,  161.31358806,
         13.23185039,   97.74780855,   -4.96138959,  -35.09563317,
       -135.10179256, -181.92611178,  100.28895261, -141.25224622,
          8.76859548, -120.07339877,  113.36272233, -273.25323425,
        131.48976626, -143.92578092, -282.8194997 ,  283.40816327,
         19.93069677, -217.94313013,  -49.4418588 ,  277.43925607,
        194.78357439,  140.36232306,  -48.04546695,   77.0365705 ,
       -263.56000458,  191.43510643, -292.19834446]), 'x_flat': array([ 0.00000000e+00,  4.20645406e+02,  2.57571025e-14, -4.20645406e+02,
        2.57571025e-14, -4.97223100e+01, -1.53377084e+02, -2.48456911e+01,
       -8.78113695e+01, -8.63965203e+01, -2.25728957e+02,  5.95217570e+01,
        7.76952465e+01, -2.32676677e+02, -8.80556597e+01,  1.50461999e+02,
        1.46459241e+02,  8.53325267e+00, -6.42191945e+01, -8.79700673e+01,
        2.00262866e+02,  8.90946094e+00, -2.26425322e+02,  2.29568143e+02,
       -2.56364340e+02, -6.12981867e+01,  6.46470883e+01,  2.83151484e+02,
       -2.57263698e+02,  7.77957750e+01, -2.37145949e+02,  2.50008612e+02,
       -1.99543801e+02, -2.36134935e+02, -8.83895203e-01, -1.03962937e+02,
        2.38349223e+02,  1.95590155e+02, -6.72610389e+00, -1.83118225e+02,
       -7.52717156e+01, -9.62563126e+01, -1.70055718e+02, -7.80660761e+01,
       -2.10725627e+02,  1.15132960e+02, -2.01294704e+02, -1.76271931e+02,
       -7.63114080e+01,  2.21798376e+01,  1.37434852e+02, -5.36946792e+01,
       -2.57218229e+02, -7.53182580e+01, -1.04512967e+02]), 'y_flat': array([ 0.00000000e+00,  0.00000000e+00,  4.20645406e+02,  5.15142050e-14,
       -4.20645406e+02, -1.58955158e+02,  2.39571207e+02,  2.16926768e+02,
       -9.58437979e+01, -2.26901496e+02,  5.11522832e+01, -9.77147291e+01,
       -3.04445626e+01, -1.69144391e+02,  4.04822071e+01, -1.02702579e+02,
        1.09469133e+02,  4.33082879e+00, -1.64832131e+02, -2.38256098e+02,
        2.24672512e+02,  1.78859170e+02,  7.38268417e+00, -2.40707350e+02,
        3.47855150e+00, -1.48905714e+02,  2.34967930e+02,  1.61455197e+02,
        1.32390227e+01,  9.77590268e+01, -4.96365242e+00, -3.51138754e+01,
       -1.35165515e+02, -1.82059267e+02,  1.00296173e+02, -1.41285966e+02,
        8.77264152e+00, -1.20124491e+02,  1.13373339e+02, -2.73497960e+02,
        1.31512628e+02, -1.43959195e+02, -2.83074517e+02,  2.83609946e+02,
        1.99378325e+01, -2.18050971e+02, -4.94587834e+01,  2.77687375e+02,
        1.94851625e+02,  1.40383573e+02, -4.80531293e+01,  7.70413994e+01,
       -2.63859412e+02,  1.91499621e+02, -2.92430684e+02]), 'x_loc': array([-405.80862082,   14.83678518, -405.80862082, -826.45402682,
       -405.80862082, -455.53093079, -559.18570435, -430.65431195,
       -493.61999028, -492.20514108, -631.53757768, -346.28686382,
       -328.11337435, -638.48529831, -493.86428055, -255.34662186,
       -259.34938024, -397.27536814, -470.02781527, -493.77868815,
       -205.5457547 , -396.89915987, -632.23394303, -176.24047781,
       -662.1729609 , -467.1068075 , -341.16153249, -122.65713705,
       -663.07231924, -328.01284581, -642.95456943, -155.8000093 ,
       -605.35242198, -641.94355625, -406.69251602, -509.77155754,
       -167.45939827, -210.2184662 , -412.5347247 , -588.92684581,
       -481.08033643, -502.0649334 , -575.86433849, -483.87469687,
       -616.53424802, -290.67566043, -607.10332477, -582.08055195,
       -482.12002882, -383.62878323, -268.37376928, -459.50330001,
       -663.02684978, -481.12687886, -510.3215881 ]), 'y_loc': array([335.46165516, 335.46165516, 756.10706116, 335.46165516,
       -85.18375084, 176.50649739, 575.03286222, 552.38842301,
       239.61785728, 108.56015952, 386.6139384 , 237.7469261 ,
       305.01709258, 166.31726454, 375.94386223, 232.75907638,
       444.93078807, 339.79248395, 170.62952405,  97.20555666,
       560.13416731, 514.32082491, 342.84433933,  94.75430505,
       338.94020666, 186.5559407 , 570.42958474, 496.91685213,
       348.70067791, 433.22068199, 330.49800274, 300.3477798 ,
       200.29613995, 153.40238852, 435.75782809, 194.17568937,
       344.23429668, 215.33716418, 448.83499394,  61.96369518,
       466.97428281, 191.50245996,  52.38713865, 619.07160094,
       355.39948771, 117.4106845 , 286.00287178, 613.14902978,
       530.3132802 , 475.84522827, 287.40852583, 412.50305454,
        71.6022433 , 526.96127642,  43.0309711 ]), 't_ext': array([140.42116437,  87.46757331, 118.2228045 , 157.90752893,
       191.85490133, 158.81993265, 134.19952119, 127.94077078,
       154.10667999, 167.56204659, 148.52591223, 145.52803918,
       137.08919335, 165.39963698, 142.72055237, 137.64952379,
       120.23783302, 139.4593894 , 160.0480858 , 168.86312743,
       110.15101816, 127.65721244, 151.53015103, 151.73565133,
       152.89388358, 158.2289204 , 120.88272482, 103.8655136 ,
       152.26072928, 127.13116678, 152.7954706 , 117.4172136 ,
       161.69189661, 166.5603205 , 133.02402524, 159.14779798,
       115.94149512, 134.3108626 , 132.58683926, 173.99374523,
       135.85243639, 159.12166336, 174.80203326, 128.01166047,
       150.03881692, 158.00506619, 154.77514432, 133.51100261,
       132.27469617, 128.87593059, 133.03846574, 138.08519675,
       173.8363439 , 132.3967438 , 175.18014642]), 'p_ext': array([3.30665111e-06, 3.07832466e-06, 3.07832466e-06, 3.07832466e-06,
       3.07832466e-06, 3.07832466e-06, 3.07832466e-06, 3.07832466e-06,
       3.07832466e-06, 2.83164684e-06, 3.30665111e-06, 3.07832466e-06,
       3.07832466e-06, 3.07832466e-06, 2.41483654e-06, 3.07832466e-06,
       2.83164684e-06, 2.83164684e-06, 3.30665111e-06, 3.07832466e-06,
       3.07832466e-06, 3.07832466e-06, 3.07832466e-06, 3.07832466e-06,
       2.83164684e-06, 3.07832466e-06, 3.07832466e-06, 3.07832466e-06,
       3.07832466e-06, 3.07832466e-06, 2.41483654e-06, 3.07832466e-06,
       2.83164684e-06, 3.07832466e-06, 2.83164684e-06, 2.83164684e-06,
       3.07832466e-06, 3.07832466e-06, 3.07832466e-06, 2.41483654e-06,
       3.07832466e-06, 3.30665111e-06, 3.30665111e-06, 3.07832466e-06,
       3.07832466e-06, 3.07832466e-06, 3.30665111e-06, 3.07832466e-06,
       2.41483654e-06, 3.52019892e-06, 3.30665111e-06, 3.07832466e-06,
       2.83164684e-06, 3.07832466e-06, 3.07832466e-06]), 't_int': array([-12.8908963 , -65.84448736, -35.08925617,   4.59546826,
        38.54284066,   5.50787198, -19.11253947, -25.37128988,
         0.79461932,  14.24998592,  -4.78614844,  -7.78402149,
       -16.22286732,  12.08757631, -10.5915083 , -15.66253688,
       -33.07422765, -13.85267127,   6.73602513,  15.55106677,
       -43.16104251, -25.65484823,  -1.78190964,  -1.57640934,
        -0.41817709,   4.91685973, -32.42933585, -49.44654707,
        -1.05133139, -26.18089389,  -0.51659007, -35.89484707,
         8.37983594,  13.24825983, -20.28803543,   5.83573731,
       -37.37056554, -19.00119807, -20.72522141,  20.68168457,
       -17.45962428,   5.8096027 ,  21.4899726 , -25.3004002 ,
        -3.27324374,   4.69300552,   1.46308365, -19.80105806,
       -21.03736449, -24.43613008, -20.27359493, -15.22686392,
        20.52428323, -20.91531686,  21.86808575]), 'p_int': array([0.56206456, 0.56206433, 0.56206433, 0.56206433, 0.56206433,
       0.56206433, 0.56206433, 0.56206433, 0.56206433, 0.56206409,
       0.56206456, 0.56206433, 0.56206433, 0.56206433, 0.56206367,
       0.56206433, 0.56206409, 0.56206409, 0.56206456, 0.56206433,
       0.56206433, 0.56206433, 0.56206433, 0.56206433, 0.56206409,
       0.56206433, 0.56206433, 0.56206433, 0.56206433, 0.56206433,
       0.56206367, 0.56206433, 0.56206409, 0.56206433, 0.56206409,
       0.56206409, 0.56206433, 0.56206433, 0.56206433, 0.56206367,
       0.56206433, 0.56206456, 0.56206456, 0.56206433, 0.56206433,
       0.56206433, 0.56206456, 0.56206433, 0.56206367, 0.56206478,
       0.56206456, 0.56206433, 0.56206409, 0.56206433, 0.56206433]), 'LENGTH_R1': 2.769707488123567, 'LENGTH_R2': 2.738976719024489, 'OFFSET_T': 153.3120606690524, 'OFFSET_P': -0.5620612559845517, 'OFFSET_X': 405.80862081533746, 'OFFSET_Y': -335.4616551587445}
                                                                             
    tests = {'ptl2flat': {'func': ptl2flat, 'inputs': ['x_ptl', 'y_ptl'], 'checks': ['x_flat', 'y_flat']},
             'flat2ptl': {'func': flat2ptl, 'inputs': ['x_flat', 'y_flat'], 'checks': ['x_ptl', 'y_ptl']},
             'flat2loc_x': {'func': flat2loc, 'inputs': ['x_flat', 'OFFSET_X'], 'checks': ['x_loc']},
             'flat2loc_y': {'func': flat2loc, 'inputs': ['y_flat', 'OFFSET_Y'], 'checks': ['y_loc']},
             'loc2flat_x': {'func': loc2flat, 'inputs': ['x_loc', 'OFFSET_X'], 'checks': ['x_flat']},
             'loc2flat_y': {'func': loc2flat, 'inputs': ['y_loc', 'OFFSET_Y'], 'checks': ['y_flat']},
             }
    
    for name, args in tests.items():
        inputs = [u[key] for key in args['inputs']]
        checks = [u[key] for key in args['checks']]
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
        print(f'\n{name} on {len(error)} points:')
        print(f' max(abs(err)) = {np.max(np.abs(errors)):8.6f}')
        print(f'    norm(err)  = {np.linalg.norm(errors):8.6f}')
        worst = np.argmax(np.abs(errors)) % len(errors[0])
        print(' Worst case:')
        for i in range(len(errors)):
            print(f'  {args["inputs"][i]}={inputs[i][worst]:9.4f} --> ' +
                  f'{args["checks"][i]}={outputs[i][worst]:9.4f}, ' +
                  f'expected={checks[i][worst]:9.4f} --> ' +
                  f'error={errors[i][worst]:8.2E}')
