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
    x_loc = rand(scale=3.5/np.sqrt(2))
    y_loc = rand(scale=3.5/np.sqrt(2))
    xy_loc_list = [[x_loc[i], y_loc[i]] for i in range(n_pts)]
    xy_ptl_list, xy_flat_list, tp_ext_list, tp_int_list = [], [], [], []     
    for i in range(n_pts):
        these_params = {key:params[key][i] for key in params}
        trans.alt = these_params
        xy_flat_list.append(trans.poslocXY_to_flatXY(xy_loc_list[i]))
        xy_ptl_list.append(trans.flatXY_to_ptlXY(xy_flat_list[i]))
        tp_ext_list.append(trans.poslocXY_to_poslocTP(xy_loc_list[i]))
        tp_int_list.append(trans.poslocTP_to_posintTP(tp_ext_list[i]))
    xy_flat = np.transpose(xy_flat_list).tolist()
    xy_ptl = np.transpose(xy_ptl_list).tolist()
    tp_ext = np.transpose(tp_ext_list).tolist()
    tp_int = np.transpose(tp_int_list).tolist()    
    u = {'x_ptl': xy_ptl[0], 'y_ptl': xy_ptl[1],
         'x_flat': xy_flat[0], 'y_flat': xy_flat[1],
         'x_loc': x_loc.tolist(), 'y_loc': y_loc.tolist(),
         't_ext': tp_ext[0], 'p_ext': tp_ext[1],
         't_int': tp_int[0], 'p_int': tp_int[1],
         }
    u.update(params)
    
    Expected values generated using:
        /code/focalplane/plate_control/trunk/petal/postransforms.py
        SVN revision r131291
    '''
    from numpy import array
    u = {'x_ptl': [48.41923375210939, -11.399312169878057, -131.82171577808455, -231.31522071929356, -91.48061600705415, 73.60930295033995, 50.644230013005185, 271.9591041696353, -42.32827202581005, 265.55863086332704, 282.7065229011477, -255.43228703635972, -119.13405762898863, -234.31957415152152, 165.28351396885023, 281.2995529303452, 30.336167560792695, 46.86529958347062, 281.19997255053005, -44.745159815392675, -121.2655496534281, 255.54922781545915, -256.5240683987907, -129.82508034724367, -75.90084484053205, -216.65903790529276, 127.75132452146269, -188.57398790605555, 114.7150375213864, -262.06005570042623, -259.2650276905908, 138.6007027943019, -32.066006765079486, -267.3064453288188, -231.79191317406776, -127.6882746960663, -208.4982183239598, -16.426547470163687, 274.47605732008935, 161.95708714625857, 124.86360054993573, -24.84879900717868, -288.3808925427442, 276.4465889306141, -292.1233649108902, -19.14723524771029, 176.47541607173164, 127.47979460222061, 281.32793798697793, -220.22825905267467], 'y_ptl': [4.03786342028203, -203.4763039865941, 199.82949227632417, -26.221316002453207, 214.473720125598, -157.97884304213397, 30.45131507306872, -244.86423510612764, -279.28576200643556, 150.89728739387473, 149.57675362483047, 23.368066064145605, -236.90034263582487, -32.56829927609965, 99.97973669195238, 147.28566062573046, 224.00280473414003, 130.87023145831867, -136.16291743623245, 180.20719508168762, 121.53535557095137, -209.8510456778891, 24.888100024084512, 195.67763247310515, 249.09530836206076, 241.33300888994015, 62.191205493623265, 143.76277967945964, -106.82245411892788, -284.514265821598, -236.7566163551534, 128.57960162124004, 197.42515980148795, 24.31196784178436, -51.77776307376631, 80.64472663485782, -122.73121773078307, -211.03894618287288, 96.77271519482878, -295.9870574058823, 122.11211369325684, 247.50158151231096, 279.70876040083544, 85.19187688248188, -107.8012930049428, 34.829134164341525, -148.42845254723554, -180.5705804066624, -107.37731099230007, -9.437755169876674], 'x_flat': [48.42001523021687, -11.40307730560457, -131.88307591413775, -231.41672549552035, -91.52089493286341, 73.62663989486916, 50.64544787632607, 272.2642004068767, -42.35605279466145, 265.7630901233028, 282.9456712167316, -255.56966712759632, -119.2027064064606, -234.42595374906594, 165.33215397892116, 281.5340031547051, 30.34867438592314, 46.872065613241354, 281.42671728728277, -44.75724823054199, -121.2931951235005, 255.78072489863183, -256.66341408701385, -129.883105809021, -75.94296822581767, -216.8476224556012, 127.77065989566887, -188.66002066810606, 114.7365075424957, -262.3924591647231, -259.53119437263985, 138.6396366474472, -32.07617914070782, -267.46439090796105, -231.89802207117893, -127.71032213142794, -208.5974533603766, -16.432433797890056, 274.66754619044497, 162.1100505064648, 124.893142218914, -24.861334527244132, -288.7799277967385, 276.63707701634394, -292.3572996289, -19.147440881355614, 176.55131541167867, 127.52990177098535, 281.5383432870875, -220.31431829401686], 'y_flat': [4.037928590703357, -203.54351119087312, 199.9225085503208, -26.23282233050598, 214.56815292937702, -158.01605125262955, 30.45204734878287, -245.13893507248466, -279.46906202853563, 151.01346643732214, 149.70328423442803, 23.38063419750558, -237.03685203732007, -32.5830850769968, 100.00915895726375, 147.4084164417955, 224.09515535495984, 130.88912543509934, -136.272711987738, 180.25588011032596, 121.5630625499308, -210.04114566523967, 24.901619416038933, 195.7650908051423, 249.2335510610587, 241.54307017977195, 62.20061823540936, 143.82836830678417, -106.84244696724518, -284.87515076206796, -236.9996754502798, 128.6157204808567, 197.48778945471216, 24.32633325611677, -51.801465718346975, 80.65865124463087, -122.78963183599173, -211.1145703761196, 96.84022890839212, -296.2666078452987, 122.14100438377896, 247.62643909759615, 280.09579594700784, 85.25057913527455, -107.88762113924966, 34.829508215306525, -148.49228932307923, -180.6415554233277, -107.4576185348284, -9.441443189086325], 'x_loc': [1.7884177758934694, -0.010918374858695826, -3.6855050875494104, 2.899715199352094, -4.363835835613318, -2.1896230265011107, 1.6020430880471055, 4.8750534442281985, -3.3641959331472444, -0.007759353088025527, 1.8765036468986676, 1.3474465084476779, 2.749086369728637, 1.7640894724801288, 3.0150930712858126, 2.990507113469418, 3.879223240071509, 3.720753814466869, 2.6839682338567856, -3.698714159184247, -1.287925077814286, -2.1537724657127377, 0.3779440216678186, -1.8528656455261303, -0.6645599535265143, -1.1260435605508847, 2.860248981672028, -4.099237776499531, -1.9885890338604495, -0.5188998822966033, -2.7259924346406397, 1.163536289313618, 3.5460929885237484, -2.62096350743897, -2.250065632206348, 0.3818620318481923, -3.851419052573654, -0.8384304364852465, -1.5953002814877022, 2.3413231841042945, -4.703898700165436, 2.7463300781349695, 4.92137196783035, 4.902061526994462, -1.4990887674326467, 0.4210972579515027, 2.8274822328086535, -4.369529130524918, -1.2816604344377567, -0.519146941675625], 'y_loc': [3.36564626451834, 2.692647003567939, -3.4536409472893346, 1.522459451414915, -1.3816500492209887, 4.103027104094897, -2.049323788824266, -2.8586672838776197, 1.1207935879545714, -1.0986149726419927, -3.2516310895631455, -2.4704457351513462, 3.831008153418326, 1.1768935023097302, 2.757106099666098, -4.316657584265735, -2.8793531889263306, 0.8020836100059306, -0.5815868307846526, 1.1631884407435458, 2.077263029540494, 4.831061556487216, -4.227983372293183, 2.1742562546016013, 4.694924168628513, -1.9599232400881599, -0.24883208693615838, 3.6486117349958227, -0.03433249933198924, -1.5279306409705236, -4.100480675903912, 2.895934764300035, 2.0708751908613547, 4.3808639959944395, -0.27685938408661115, 0.45549705570227444, -2.7405259020706914, 1.4833411567978676, 2.1529573547996277, 0.6015570712262719, 3.4049857181802925, -1.6664745633556786, 2.579587004085455, -3.338100512403188, 4.325925959759568, -2.4475234608878793, 2.0984395802144795, 3.8563845710625775, 1.0796836796697389, -1.8854376755957616], 't_ext': [(15.74747355256734, 99.14157478296411), (23.033825610182152, 125.7095682147544), (-174.8027363321766, 67.2851098362648), (-28.022114199491476, 114.82594021090986), (163.71580090868594, 70.47939539312411), (-275.4075085624898, 70.03489020543378), (-113.56629359371615, 129.7483536784539), (-30.386801261226896, 3.6222548091772e-06), (113.67609862695032, 105.05983028554677), (-154.28444969369102, 158.7178352457939), (247.20627764560416, 108.98337038912942), (-123.30823521683946, 128.5091329338098), (14.788562219488956, 72.27118821461097), (-30.02528542522458, 137.41095907616702), (-9.332071916187447, 103.09160254249822), (-74.32124189213773, 38.119980791754635), (-70.32167878459704, 62.3487425131605), (-38.52204828507229, 108.25738075743841), (-80.70856524321928, 130.4457944398235), (106.2573628197106, 106.14802713196222), (59.72849616106975, 134.16625897835493), (114.0281544978021, 2.9575586669421963e-06), (236.58696936579233, 88.69196621745813), (64.90818761272509, 128.9632042147988), (64.85956500948458, 63.37412521027217), (169.25120807128025, 138.56250908687443), (-67.02573087583596, 124.34037339408016), (115.74898536295173, 43.71378101142817), (-257.73329290084973, 141.8737959043669), (173.2409215635405, 146.74746664719945), (-165.01043174070043, 81.44042119340772), (16.648143077832803, 114.79063692446815), (-21.721177837812657, 95.49584919876797), (85.01335808940605, 72.71805560199621), (105.24472208632062, 137.37632487127428), (-66.22981995578151, 169.54724302828777), (182.63363491416894, 72.40260102485767), (43.70786927431703, 143.8703836162625), (67.26246079852227, 124.68797008799174), (-48.809398892188035, 132.23914511029628), (119.1433941429993, 52.61685819494954), (-96.08400003900464, 114.29320796388608), (9.711378599140723, 39.90223616766415), (-54.936470022652244, 41.207713602759526), (-290.5894546228789, 86.5611335766062), (213.90267798069655, 134.32735542394144), (342.6978141479232, 107.7980026159942), (112.54203741970586, 49.58448240953362), (80.77320972499133, 147.15340377276164), (-171.92631302753531, 140.57618364194863)], 'p_ext': [False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, False, False, False, False, False, False, True, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False], 't_int': [array([-95.79999935, -12.40589812]), array([ 62.51227736, 165.18801996]), array([-165.44425901,   76.64358716]), array([ 22.97378862, 165.82184303]), array([ 59.79959579, -33.43680973]), array([-176.62422426,  168.81817451]), array([-124.46070196,  118.85394531]), array([-165.60451937, -135.21771449]), array([104.43511115,  95.81884281]), array([-38.26474459, 274.73754035]), array([145.81271603,   7.58980877]), array([-35.71608401, 216.10128414]), array([ 98.21703575, 155.69966175]), array([-16.70004012, 150.73620439]), array([-175.5839347 ,  -63.16026024]), array([-43.50730923,  68.93391345]), array([-165.91890829,  -33.24848699]), array([-56.35293213,  90.42649691]), array([ 52.78002257, 263.93438226]), array([52.72868083, 52.61934515]), array([ 36.01972221, 110.45748503]), array([ 72.19270954, -41.835442  ]), array([ 89.10244944, -58.79255371]), array([108.17997492, 172.23499152]), array([-80.28562794, -81.77106774]), array([ -5.03171956, -35.72041854]), array([ 69.48308197, 260.84918624]), array([ -36.42040247, -108.45560682]), array([-147.14944581,  252.45764299]), array([131.01009848, 104.51664357]), array([-66.76513622, 179.68571671]), array([-121.80689092,  -23.66439708]), array([ 14.37543467, 131.59246171]), array([166.8619168 , 154.56661431]), array([ 68.78943421, 100.921037  ]), array([-135.69905217,  100.07801081]), array([109.43689187,  -0.79414202]), array([-46.41177832,  53.75073602]), array([109.7780033 , 167.20351259]), array([ -2.79712624, 178.25141776]), array([ -56.88004663, -123.40658258]), array([ -5.88051474, 204.49669326]), array([ 98.42032105, 128.61117862]), array([-3.17740457, 92.96677906]), array([-147.33593973,  229.81464847]), array([84.37565361,  4.80033105]), array([172.13696514, -62.76284639]), array([139.04278509,  76.08523008]), array([-52.34710808,  14.03308597]), array([-134.37097085,  178.13152582])], 'p_int': [-17.576782782432193, 16.18377877748909, 0.8447100967219301, -15.078820683399368, -8.835027748566407, -4.442588089864015, 14.567148606864553, -2.4642971925032384, 10.51133948305357, -6.541618407525962, -5.958405359349723, -8.500037081389209, -17.40311585665701, 18.834298640559123, -5.0157658988231235, 3.436060587941867, -7.392902716996446, 11.030482050448231, -17.641800190065933, 10.017901440772356, -15.363072074808732, 7.7143293863898235, -18.105925111829208, -5.042144030718485, -0.7567930381294996, -5.006530742387345, 0.0451230661866342, 1.7062728451385656, -11.302474194541743, -19.330754072919483, 12.074941284872288, -19.820973431145198, -0.48547636929666904, 19.714354245912862, 17.912496481138195, 1.266357555827442, -3.925403355191537, -0.6064387496587464, -13.516503058506974, 4.857680960430089, -1.4807270098155412, -19.74381363901138, 10.902661952366643, 0.18749804381741342, -19.65281408364652, -7.042278571717033, -0.5357872801910801, -0.23718002251590153, -4.86335535677346, 3.3574360429998285], 'LENGTH_R1': array([3.07788448, 2.82781894, 2.68315602, 3.09647306, 2.89737147,
       2.9460441 , 3.14050025, 2.73953597, 3.08531151, 3.01610202,
       3.29917385, 3.30017764, 2.67604725, 3.00731458, 3.27449244,
       2.78150513, 2.61174183, 3.38294325, 3.18517398, 3.08597053,
       3.24228662, 2.65153602, 3.26070588, 3.2860385 , 2.66629794,
       3.15986272, 3.07823386, 2.86321664, 2.87415287, 2.74266612,
       3.20373785, 3.07193761, 2.83926993, 3.20563793, 2.76248943,
       2.62651688, 3.16111511, 2.68134353, 2.96331753, 3.0487321 ,
       3.39254619, 2.67837966, 3.2380708 , 3.15630734, 3.34667475,
       3.22937943, 2.98856799, 3.0592293 , 3.08787324, 2.96080794]), 'LENGTH_R2': array([2.78937964, 3.05700462, 3.36672092, 2.98185894, 2.70535   ,
       2.73065292, 2.9755338 , 2.65836905, 2.72453389, 2.71775764,
       3.16164057, 3.17277817, 3.15208924, 2.81007823, 3.29517647,
       2.77443564, 3.02902737, 3.10097599, 3.35714331, 3.35765113,
       3.01044603, 2.60051586, 2.64439363, 3.3438964 , 2.90417149,
       3.22686266, 3.07163383, 3.04916397, 3.15920936, 2.87852652,
       3.29256784, 2.68899999, 3.25115366, 3.1332975 , 3.31328025,
       2.93821736, 2.68637118, 2.80116558, 2.80136702, 2.9149424 ,
       3.0835704 , 3.18994863, 2.66954599, 3.17965463, 2.92990468,
       3.16804762, 2.98738532, 3.35884899, 2.65158402, 2.82471083]), 'OFFSET_T': array([ 111.5474729 ,  -39.47845175,   -9.35847732,  -50.99590282,
        103.91620512,  -98.7832843 ,   10.89440837,  135.21771811,
          9.24098748, -116.0197051 ,  101.39356162,  -87.59215121,
        -83.42847353,  -13.32524531,  166.25186279,  -30.81393266,
         95.5972295 ,   17.83088385, -133.48858782,   53.52868198,
         23.70877395,   41.83544496,  147.48451992,  -43.27178731,
        145.14519295,  174.28292763, -136.50881284,  152.16938783,
       -110.58384709,   42.23082308,  -98.24529552,  138.455034  ,
        -36.09661251,  -81.84855871,   36.45528787,   69.46923221,
         73.19674304,   90.1196476 ,  -42.5155425 ,  -46.01227265,
        176.02344078,  -90.2034853 ,  -88.70894245,  -51.75906545,
       -143.25351489,  129.52702437,  170.56084901,  -26.50074767,
        133.12031781,  -37.55534218]), 'OFFSET_P': array([ 17.57678278, -16.18377878,  -0.8447101 ,  15.07882068,
         8.83502775,   4.44258809, -14.56714861,   3.46429719,
       -10.51133948,   6.54161841,   5.95840536,   8.50003708,
        17.40311586, -18.83429864,   5.0157659 ,  -3.43606059,
         7.39290272, -11.03048205,  17.64180019, -10.01790144,
        15.36307207,  -6.71432939,  18.10592511,   5.04214403,
         0.75679304,   5.00653074,  -0.04512307,  -1.70627285,
        11.30247419,  19.33075407, -12.07494128,  19.82097343,
         0.48547637, -19.71435425, -17.91249648,  -1.26635756,
         3.92540336,   0.60643875,  13.51650306,  -4.85768096,
         1.48072701,  19.74381364, -10.90266195,  -0.18749804,
        19.65281408,   7.04227857,   0.53578728,   0.23718002,
         4.86335536,  -3.35743604]), 'OFFSET_X': array([  46.63159745,  -11.39215893, -128.19757083, -234.31644069,
        -87.1570591 ,   75.81626292,   49.04340479,  267.38914696,
        -38.99185686,  265.77084948,  281.06916757, -256.91711364,
       -121.95179278, -236.19004322,  162.31706091,  278.54349604,
         26.46945115,   43.1513118 ,  278.74274905,  -41.05853407,
       -120.00527005,  257.93449736, -257.04135811, -128.03024016,
        -75.27840827, -215.7215789 ,  124.91041091, -184.56078289,
        116.72509658, -261.87355928, -256.80520194,  137.47610036,
        -35.62227213, -264.8434274 , -229.64795644, -128.09218416,
       -204.74603431,  -15.59400336,  276.26284647,  159.76872732,
        129.59704092,  -27.60766461, -293.70129976,  271.73501549,
       -290.85821086,  -19.56853814,  173.72383318,  131.8994309 ,
        282.82000372, -219.79517135]), 'OFFSET_Y': array([   0.67228233, -206.23615819,  203.3761495 ,  -27.75528178,
        215.94980298, -162.11907836,   32.50137114, -242.28026779,
       -280.58985562,  152.11208141,  152.95491532,   25.85107993,
       -240.86786019,  -33.75997858,   97.25205286,  151.72507403,
        226.97450854,  130.08704183, -135.69112516,  179.09269167,
        119.48579952, -214.87220722,   29.12960279,  193.59083455,
        244.53862689,  243.50299342,   62.44945032,  140.17975657,
       -106.80811447, -283.34722012, -232.89919477,  125.71978572,
        195.41691426,   19.94546926,  -51.52460633,   80.20315419,
       -120.04910593, -212.59791153,   94.68727155, -296.86816492,
        118.73601867,  249.29291366,  277.51620894,   88.58867965,
       -112.2135471 ,   37.27703168, -150.5907289 , -184.49793999,
       -108.53730221,   -7.55600551])}
                                                                             
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
