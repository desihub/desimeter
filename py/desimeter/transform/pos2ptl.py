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

_default_tp_ranges = [[-180. + xy2tp.epsilon, 180.], # theta min, max
                      [-20., 200.]] # phi min, max

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
    u_flat = _add_offset(u_loc, +u_offset)
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
        tp, unreach = xy2tp.xy2tp(xy=[x_loc[i], y_loc[i]],
                                  r=[r1[i], r2[i]],
                                  ranges=_default_tp_ranges)
        t_ext.append(tp[0])
        p_ext.append(tp[1])
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
    x = r1 * np.cos(t) + r2 * np.cos(t_plus_p)
    y = r1 * np.sin(t) + r2 * np.sin(t_plus_p)
    return x, y

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
    u_ext = _add_offset(u_int, +u_offset)
    return u_ext

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
    '''Sample code to generate test values below:
    
    import postransforms # c.f. SVN... /code/focalplane/plate_control/trunk/petal
    import numpy as np
    trans = postransforms.PosTransforms(stateless=True)
    r_max = 420
    n_pts = 50
    def rand(scale=1, offset=0):
        return scale * 2 * (np.random.rand(n_pts) - 0.5) + offset
    params = {'LENGTH_R1': rand(scale=0.4, offset=3),
              'LENGTH_R2': rand(scale=0.4, offset=3),
              'OFFSET_T': rand(scale=180),
              'OFFSET_P': rand(scale=20),
              'OFFSET_X': rand(scale=r_max/np.sqrt(2)),
              'OFFSET_Y': rand(scale=r_max/np.sqrt(2)),
              }
    max_patrol = 6.5 # a bit large, to provide some "unreachable" targs
    x_loc = rand(scale=max_patrol/np.sqrt(2))
    y_loc = rand(scale=max_patrol/np.sqrt(2))
    xy_loc_list = [[x_loc[i], y_loc[i]] for i in range(n_pts)]
    xy_ptl_list, xy_flat_list, tp_ext_list, tp_int_list, unreachable = [], [], [], [], []   
    for i in range(n_pts):
        these_params = {key:params[key][i] for key in params}
        trans.alt = these_params
        xy_flat_list.append(trans.poslocXY_to_flatXY(xy_loc_list[i]))
        xy_ptl_list.append(trans.flatXY_to_ptlXY(xy_flat_list[i]))
        tp_ext_vals, unreach = trans.poslocXY_to_poslocTP(xy_loc_list[i])
        tp_ext_list.append(tp_ext_vals)
        unreachable.append(unreach)
        tp_int_list.append(trans.poslocTP_to_posintTP(tp_ext_list[i]))
    xy_flat = np.transpose(xy_flat_list).tolist()
    xy_ptl = np.transpose(xy_ptl_list).tolist()
    tp_ext = np.transpose(tp_ext_list).tolist()
    tp_int = np.transpose(tp_int_list).tolist()    
    u = {'x_ptl': xy_ptl[0], 'y_ptl': xy_ptl[1],
         'x_flat': xy_flat[0], 'y_flat': xy_flat[1],
         'x_loc': x_loc, 'y_loc': y_loc,
         't_ext': tp_ext[0], 'p_ext': tp_ext[1],
         'unreachable': unreachable,
         't_int': tp_int[0], 'p_int': tp_int[1],
         }
    u.update(params)
    
    Expected values generated using:
        /code/focalplane/plate_control/trunk/petal/postransforms.py
        SVN revision r131291
    '''
    from numpy import array
    u = {'x_ptl': [90.63460879628963, -21.0353545303502, -47.267561963845786, 28.3686292912251, -228.79547014820622, 66.1383727336639, 199.04726006546045, -76.44685203581305, -258.10978039985775, -217.12862492188214, 175.1136041021613, 64.48180930846004, -294.86137899473295, -204.7854791869051, -133.70716660698636, -182.67183923004646, -259.7467064797599, -275.73119564039706, -161.93543445552217, 52.57079523592689, 256.00253288473294, 250.02664724491063, -81.95265251605346, 169.2264394686997, -291.4513766372221, -285.39099825442025, -238.28371976710633, -9.89887449605023, 230.20285298284358, -143.79117987131076, -44.73260485827212, -27.51284107474551, 241.99526553291778, -59.44705324441423, 233.71625764129325, 287.1644610489922, -235.85530745081059, -168.62540367986358, 225.7020855085132, 159.54249763529165, 221.60521934990248, -262.9077776568561, -175.71878315566335, 250.77234820995898, 213.08361035451185, -209.88515625823098, 198.71695660127807, -175.41832100153624, 224.02270236886767, -282.0313503093771], 'y_ptl': [-165.86495978007844, -198.223226937375, 161.49898902360368, 210.56638617537405, -254.75632406661563, -265.41954817389126, 204.19779015294276, 181.36807906388776, 80.58956648306373, -113.24681658438668, -228.31171195146933, 196.39373190717015, -67.7360155166228, 200.3932461842164, -106.87516022105476, 146.18601442824985, -185.08918646421012, 254.0698890383764, 139.6908038054167, -129.940626365478, -177.1083036710536, -96.56552662958845, 123.67204780686092, -221.70995114128633, -218.63700353900734, 151.50406466035284, 150.93958870748685, 52.43772792334824, 50.92201183685703, 42.02316021043662, 239.90646387171918, 181.9142411098302, -32.56731608547286, 127.92265844195882, -219.46803505851815, 113.03794103023415, 109.2335739658909, -176.9632216662303, 255.09567050305708, 137.8618526956137, 27.79841941381922, 184.96438877522945, -256.31487956793205, -87.24479476621836, -212.46846301576696, -258.2098806231467, -273.7099844274206, -236.7283788431389, 228.95534376081483, 285.73157183235094], 'x_flat': [90.66005544065625, -21.04197918698579, -47.27787877359343, 28.378876619772498, -229.01830805647174, 66.17900654019267, 199.18045393989289, -76.47027702246078, -258.2646326317088, -217.23463052367026, 175.23295695088441, 64.50376750032052, -295.0840315131238, -204.92385494315747, -133.7374508568168, -182.75285093128167, -259.96517860230335, -276.0573953193006, -161.99475248460993, 52.578522746147144, 256.20742816627677, 250.17394601540255, -81.9662699092042, 169.33465234564127, -291.77523725573343, -285.6375131384581, -238.43964564147527, -9.899067871542789, 230.30661446033014, -143.81557772159059, -44.7542859051917, -27.520127543619584, 242.1126936194129, -59.45591118018986, 233.91494166574628, 287.3902701057917, -235.98570139105402, -168.70741414798883, 225.91949813502094, 159.59919261561288, 221.69427195303675, -263.1324876703682, -175.85893273615622, 250.91724209682766, 213.24274634794543, -210.07776850292092, 198.90561014729045, -175.54381638884297, 224.21279860184, -282.42099135739875], 'y_flat': [-165.91152815721895, -198.2856533069313, 161.53423844783816, 210.64244705649827, -255.00444683904493, -265.58261548449593, 204.33443054081883, 181.4236542081514, 80.63791596526656, -113.30210545726997, -228.4673232379709, 196.46061047791224, -67.78716360014583, 200.52865406127245, -106.89936710640677, 146.25084531721197, -185.24486438588505, 254.37046263154878, 139.74197347801663, -129.9597266570707, -177.25005483005785, -96.62241649913307, 123.69259736621063, -221.85172491927813, -218.8799528639087, 151.6349307604031, 151.0383591453879, 52.43875229955079, 50.94496438986112, 42.03029051406173, 240.02274199404158, 181.96241906615225, -32.58311935997816, 127.94171961708663, -219.65460655715623, 113.12682734559978, 109.29396436486152, -177.04928720756084, 255.34139716353465, 137.91084324760752, 27.809590277118748, 185.12247976314683, -256.51931088825194, -87.29520398203888, -212.6271396022912, -258.4468406140035, -273.96983320944895, -236.89773583743573, 229.1496256255883, 286.12632492976184], 'x_loc': [-0.2648584821821851, 4.232684252511278, -1.0441013686599954, -3.884184158320891, 1.174474322774752, 1.6598752657757154, 2.9449894543135153, 2.827340887420053, -4.373243522544701, 0.1092304143398498, -4.026918199058748, -4.519956567442952, -0.7971814117061122, 0.3917816872889854, -3.599617619765202, -4.39356165126558, -2.5932581321977035, -3.000816451082693, -0.5817304908246635, 4.065840079733202, 3.5704007998292058, -3.075941675814469, -0.30247961621612324, -1.4224725994093255, 1.1087945076688375, -1.5213739932223451, 2.9562265807583787, -2.2353006929695143, -2.185317003291131, -3.4106117980291066, -0.3991297312171695, 1.821170020666991, -2.31428791468542, -2.22645211601048, -2.588062319840614, -1.8988642566692768, 3.1747622043421324, 3.145081992488124, 2.2810461442318104, -0.22362666500581874, -3.756738151500611, 0.6893888483137933, -3.3052094590021524, 2.4878482871137613, 2.8035355679501457, -0.34542572830483936, 1.9966369212183974, -1.3423001009806392, -4.03898377535665, -4.265174105534376], 'y_loc': [-2.090606319076951, -2.808580046323973, -2.2565415631808525, 4.539613166371122, -1.651508168981695, -0.6468824151792065, 3.7282658140000615, -0.582496116514988, 3.8375589640650447, 0.32275997837194553, -0.7180633075963913, 1.2089057860922006, 1.0206771903390561, -3.9071041753091893, -4.118268823060122, 0.8331579251307693, -0.8010424620569058, -2.6976964408569883, -1.7232525301643788, -2.434747378172805, -0.10211838120791084, -3.709223769618239, -2.74362094927032, -4.313369158842518, 3.846632839318913, -2.5514427051939643, 0.8663681567695172, 0.6291684314334793, -3.6790025670946895, 2.95674096565261, 3.0561409390427277, -4.008007343065422, 2.8673873675863684, 2.3275624041627894, 1.7673937421561687, 3.204230009904061, 0.2842454955401577, -2.9932766684854513, 0.6883926065203279, -1.9727339285413468, -1.1156077753205742, -4.298151031602346, 4.218030110257953, 1.7557971507295564, -4.460309401762153, 4.261680666397903, 3.600897002856988, -0.4642650494114469, -1.662202876834141, 4.024824066884327], 't_ext': [192.36322973693123, -61.09134556489912, 170.58923081625196, 113.429216705318, -123.46120176999416, -93.19741721522685, 8.320602515771075, -77.27736641811114, 138.73276774125924, 251.3028420041431, 148.95856017623817, -238.18070707593634, 55.48515820589272, -137.2756864319763, -156.52430063084714, 130.78255268339495, 135.32961362351148, 174.3613190547112, -184.68886461164345, -59.634817766433564, -57.99158324604043, -165.89996887008138, -154.3524537339689, 218.16216720569747, 31.64016186223705, 186.97712298775411, 310.673911896051, -264.5627606302283, 188.6376588798409, 99.81033697623997, 41.59196710505199, -111.65401392426087, -286.1486220834214, -276.80585623495546, 87.37426087083936, 72.54348521817181, 312.30758293676104, 274.24173546172153, 315.2382921464383, 205.11945266789724, 151.3284931599138, 233.04632938602253, 101.26508234289543, -20.13184757398297, -87.0963438877426, -302.07392696403946, 21.03261763854648, -225.0922005209554, 156.77710895946632, -244.427234056898], 'p_ext': [139.46628226911852, 55.12204960535054, 131.49715408568383, 36.92757900812411, 141.63310629107505, 145.72092883785035, 75.43858644101064, 125.25270243620808, 2.4148365394514667e-06, 179.99999879258172, 91.97021394770464, 79.67515195583277, 156.33459252026879, 94.79926884884728, 46.87319193010012, 79.98983761415178, 120.82948630569112, 87.70479547290991, 148.28952357534735, 57.398910435459875, 102.54996690912598, 67.02128630503006, 119.21087866590187, 66.68991305131948, 99.04194911578335, 123.72237480268426, 117.88123558004625, 129.1444299561245, 92.73511851753386, 85.86622922503327, 123.08813271537956, 83.89711382761676, 110.48158560066322, 111.21910266196883, 119.34104768659817, 97.54619374136313, 116.25478395085221, 84.35042027367021, 126.77330028705198, 143.65302844071718, 102.20361888613401, 91.55142631291095, 52.874924133346894, 113.94772029303297, 55.946772740496, 72.19639957380541, 87.02256711056674, 151.99337226343968, 95.48681385903177, 40.5366026189824], 'unreachable': [False, False, False, False, False, False, False, False, True, True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False], 't_int': [53.985004034628645, 69.94156902867238, 148.88083008570993, -37.34797998546367, -125.8627702041241, -133.03609166091903, -48.234354759322756, -69.08448902537333, 158.83813489467957, 79.25698781119303, 146.27103809475486, -178.8720099066873, 116.70724742506476, -20.666103587012756, -56.456497323465996, 56.71805370557489, 177.45172267359376, 80.17498672532719, -62.11849326044408, -1.2904724333092545, -100.19213710482487, -159.48180470300474, -104.23563389231003, 175.320663647074, 175.26478282439413, 101.12593636723868, 139.74734689806436, -171.32703558033816, 108.40703107181619, 18.63532603987653, 97.15844131263991, -5.0249230628106005, -141.04015493354987, -122.4610936497935, -38.936971724814285, 95.97688881919989, 172.0003963722951, 163.71863447684248, 140.54914548021856, 65.69143313036756, -9.126124806132879, 96.84376888506915, 96.93339739608297, -41.88007849125526, -133.1752491337256, -124.22345555370765, -148.94260776292694, -171.3487196391079, 177.70869157538618, -119.94705522718849], 'p_int': [150.36693175794937, 50.95695229131336, 137.2483077213819, 55.40312990957678, 128.99861998135137, 147.97494674062622, 65.07381756578592, 136.1086688424625, -9.244722170360282, 182.19077053842793, 89.12170495526937, 92.55886225909583, 151.18137235637792, 100.43540705586605, 53.53399833123595, 96.5608759366545, 103.92322183405389, 70.77673887259978, 167.39675528372058, 71.5005112232542, 101.11331074740934, 76.5783954635254, 110.78843209180978, 85.91655691282058, 84.54521062251507, 107.6096529957897, 126.86912858942307, 124.73686968250706, 111.62134291454382, 75.79742228559819, 112.07378487512743, 87.41228031059404, 108.14624918774496, 127.8289032045208, 105.1386267807014, 82.83828050935236, 133.12940046797743, 84.54659075751202, 136.3778620306517, 137.99466201600364, 109.88511145998636, 85.13951847377506, 59.47864088699823, 111.24114925802168, 64.63353088123681, 56.883353354833595, 83.78127957167123, 159.64529105292723, 114.03517156442177, 28.729613768191868], 'LENGTH_R1': array([3.02818875, 2.86841966, 2.78145968, 3.36943122, 3.11815308,
       3.03763295, 2.60595312, 3.04949516, 2.79021019, 2.60485982,
       3.17256341, 2.82676468, 3.20779766, 2.62634876, 2.74720677,
       3.00954226, 2.70914065, 2.60181604, 3.29558948, 2.69961428,
       2.64098502, 2.67918528, 2.76990461, 2.70103857, 3.3904446 ,
       3.38711398, 2.75457162, 2.60089103, 2.8711473 , 3.2881487 ,
       3.39219617, 2.71403644, 3.23878257, 3.01275366, 3.14575201,
       2.85433109, 3.17921174, 2.92933041, 2.70064341, 3.33829667,
       3.36236278, 3.10521114, 2.95239398, 2.84400446, 2.85694003,
       2.60704212, 3.01843041, 3.02251799, 3.3560252 , 3.00433129]), 'LENGTH_R2': array([3.05497838, 2.86156433, 3.2001744 , 2.92760317, 3.04563812,
       3.0065682 , 3.37113115, 3.22017238, 2.6426979 , 3.39650953,
       2.69333712, 3.2560133 , 3.0772613 , 3.14710994, 3.21090109,
       2.82558204, 2.78655314, 2.98187981, 3.35796686, 2.70321738,
       3.04625084, 3.09360213, 2.68366128, 2.73582802, 2.72710083,
       2.82269707, 3.17535102, 2.79234892, 3.31284658, 2.8643998 ,
       3.04430033, 3.18968389, 3.22430575, 2.6674357 , 3.05868541,
       2.79684544, 2.83128597, 2.92929937, 2.61545707, 2.85356186,
       2.84555468, 3.13597064, 3.03208943, 2.74072725, 3.10674313,
       2.6842578 , 2.64797464, 2.72253732, 3.13449725, 3.2464922 ]), 'OFFSET_T': array([ 138.3782257 , -131.03291459,   21.70840073,  150.77719669,
          2.40156843,   39.83867445,   56.55495728,   -8.19287739,
        -20.10536715,  172.04585419,    2.68752208,  -59.30869717,
        -61.22208922, -116.60958284, -100.06780331,   74.06449898,
        -42.12210905,   94.18633233, -122.57037135,  -58.34434533,
         42.20055386,   -6.41816417,  -50.11681984,   42.84150356,
       -143.62462096,   85.85118662,  170.926565  ,  -93.23572505,
         80.23062781,   81.17501094,  -55.56647421, -106.62909086,
       -145.10846715, -154.34476259,  126.3112326 ,  -23.4334036 ,
        140.30718656,  110.52310098,  174.68914667,  139.42801954,
        160.45461797,  136.2025605 ,    4.33168495,   21.74823092,
         46.07890525, -177.85047141,  169.9752254 ,  -53.74348088,
        -20.93158262, -124.48017883]), 'OFFSET_P': array([-10.90064949,   4.16509731,  -5.75115364, -18.4755509 ,
        12.63448631,  -2.2540179 ,  10.36476888, -10.85596641,
         9.24472459,  -2.19077175,   2.84850899, -12.8837103 ,
         5.15322016,  -5.63613821,  -6.6608064 , -16.57103832,
        16.90626447,  16.9280566 , -19.10723171, -14.10160079,
         1.43665616,  -9.55710916,   8.42244657, -19.22664386,
        14.49673849,  16.11272181,  -8.98789301,   4.40756027,
       -18.8862244 ,  10.06880694,  11.01434784,  -3.51516648,
         2.33533641, -16.60980054,  14.20242091,  14.70791323,
       -16.87461652,  -0.19617048,  -9.60456174,   5.65836642,
        -7.68149257,   6.41190784,  -6.60371675,   2.70657104,
        -8.68675814,  15.31304622,   3.24128754,  -7.65191879,
       -18.54835771,  11.80698885]), 'OFFSET_X': array([  90.92491392,  -25.27466344,  -46.2337774 ,   32.26306078,
       -230.19278238,   64.51913127,  196.23546449,  -79.29761791,
       -253.89138911, -217.34386094,  179.25987515,   69.02372407,
       -294.2868501 , -205.31563663, -130.13783324, -178.35928928,
       -257.37192047, -273.05657887, -161.41302199,   48.51268267,
        252.63702737,  253.24988769,  -81.66379029,  170.75712495,
       -292.88403176, -284.11613915, -241.39587222,   -7.66376718,
        232.49193146, -140.40496592,  -44.35515617,  -29.34129756,
        244.42698153,  -57.22945906,  236.50300399,  289.28913436,
       -239.1604636 , -171.85249614,  223.63845199,  159.82281928,
        225.4510101 , -263.82187652, -172.55372328,  248.42939381,
        210.43921078, -209.73234277,  196.90897323, -174.20151629,
        228.25178238, -278.15581725]), 'OFFSET_Y': array([-163.82092184, -195.47707326,  163.79078001,  206.10283389,
       -253.35293867, -264.93573307,  200.60616473,  182.00615032,
         76.800357  , -113.62486544, -227.74925993,  195.25170469,
        -68.80784079,  204.43575824, -102.78109828,  145.41768739,
       -184.44382192,  257.06815907,  141.46522601, -127.52497928,
       -177.14793645,  -92.91319273,  126.43621832, -217.53835576,
       -222.7265857 ,  154.18637347,  150.17199099,   51.80958387,
         54.62396696,   39.07354955,  236.96660105,  185.97042641,
        -35.45050673,  125.61415721, -221.4220003 ,  109.92259734,
        109.00971887, -174.05601054,  254.65300456,  139.88357718,
         28.92519805,  189.42063079, -260.737341  ,  -89.05100113,
       -208.1668302 , -262.70852128, -277.57073021, -236.43347079,
        230.8118285 ,  282.10150086])}
                                                                             
    tests = {'ptl2flat': {'func': ptl2flat, 'inputs': ['x_ptl', 'y_ptl'], 'checks': ['x_flat', 'y_flat']},
             'flat2ptl': {'func': flat2ptl, 'inputs': ['x_flat', 'y_flat'], 'checks': ['x_ptl', 'y_ptl']},
             'flat2loc_x': {'func': flat2loc, 'inputs': ['x_flat', 'OFFSET_X'], 'checks': ['x_loc']},
             'flat2loc_y': {'func': flat2loc, 'inputs': ['y_flat', 'OFFSET_Y'], 'checks': ['y_loc']},
             'loc2flat_x': {'func': loc2flat, 'inputs': ['x_loc', 'OFFSET_X'], 'checks': ['x_flat']},
             'loc2flat_y': {'func': loc2flat, 'inputs': ['y_loc', 'OFFSET_Y'], 'checks': ['y_flat']},
             'ext2loc': {'func': ext2loc, 'inputs': ['t_ext', 'p_ext', 'LENGTH_R1', 'LENGTH_R2'], 'checks': ['x_loc', 'y_loc']},
             'loc2ext': {'func': loc2ext, 'inputs': ['x_loc', 'y_loc', 'LENGTH_R1', 'LENGTH_R2'], 'checks': ['t_ext', 'p_ext']},
             'ext2int_t': {'func': ext2int, 'inputs': ['t_ext', 'OFFSET_T'], 'checks': ['t_int']},
             'ext2int_p': {'func': ext2int, 'inputs': ['p_ext', 'OFFSET_P'], 'checks': ['p_int']},
             'int2ext_t': {'func': int2ext, 'inputs': ['t_int', 'OFFSET_T'], 'checks': ['t_ext']},
             'int2ext_p': {'func': int2ext, 'inputs': ['p_int', 'OFFSET_P'], 'checks': ['p_ext']},

             }
    
    all_out = {}
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
