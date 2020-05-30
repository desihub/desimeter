# -*- coding: utf-8 -*-
"""
Best-fit calculation module for positioner calibration parameters.
"""

import math
import scipy.optimize
import numpy as np

# imports below require <path to desimeter>/py' to be added to system PYTHONPATH.
import desimeter.transform.pos2ptl as pos2ptl
from desimeter.circles import robust_fit_circle

# default parameter values and bounds
default_values = {'LENGTH_R1': 3.0,
                  'LENGTH_R2': 3.0,
                  'OFFSET_T': 0.0,
                  'OFFSET_P': 0.0,
                  'OFFSET_X': 0.0,
                  'OFFSET_Y': 0.0,
                  'SCALE_T': 1.0,
                  'SCALE_P': 1.0,
                  }

default_bounds = {'LENGTH_R1': (2.5, 3.5),
                  'LENGTH_R2': (2.5, 3.5),
                  'OFFSET_T': (-200.0, 200.0),#(-179.999999, 180.0),
                  'OFFSET_P': (-30.0, 30.0),
                  'OFFSET_X': (-500.0, 500.0),
                  'OFFSET_Y': (-500.0, 500.0),
                  'SCALE_T': (0.0, 1.0),
                  'SCALE_P': (0.0, 1.0),
                  }

# static vs dynamic parameters
'''I use the term "static" in this module to signify the common six geometry
values for each positioner, giving its kinematic arm lengths, center location,
and angle of mounting. To first order, these are reasonably "permanent" values.

Contrast this with the two "dynamic" values, used for effects like effective
output ratio while moving. These ('SCALE_T' and 'SCALE_P)' are mathematically
like 'GEAR_CALIB_T' and 'GEAR_CALIB_P' elsewhere in the positioner control code.
(However, without the misleading "gear" keyword. The need for such calibrations
is generally due to electrical failures, not mechanical gear issues.)

i.e.
SCALE_T = GEAR_CALIB_T = (obsT2 - obsT1) / (posintT2 - posintT1)
        = actual_theta_output_shaft_rotations / (commanded_theta_magnetic_field_rotations * GEAR_RATIO_T)

SCALE_T and SCALE_P are not part of the common coordinates transformation
system. The reason is that they can only be applied to delta motions of the
gearmotor output shafts. They have no meaning in absolute coordinates, and it's
easy to get confused on this fact. Therefore in operation of the instrument,
these parameters are handled very carefully, in only one low-level place:
delta shaft distance calculations done within the PosModel.Axis() class.
'''
static_keys = ['LENGTH_R1', 'LENGTH_R2', 'OFFSET_T', 'OFFSET_P', 'OFFSET_X', 'OFFSET_Y']
dynamic_keys = ['SCALE_T', 'SCALE_P']
all_keys = static_keys + dynamic_keys

def fit_params(posintT, posintP, ptlX, ptlY, gearT, gearP,
               mode='static',
               nominals=None,
               bounds=None,
               keep_fixed=None,
               ):
    '''Best-fit function for parameters used in the transformation between
    internally-tracked (theta,phi) and externally measured (x,y).

    The natural order of usage is in two passes:
        1. mode='static' --> best-fit static parameters
        2. mode='dynamic', feeding in results from (1)

    INPUTS:
        posintT  ... list of theta angles as internally-tracked (a.k.a. POS_T, a.k.a. t_ext)
        posintP  ... list of phi angles as internally-tracked (a.k.a. POS_P, a.k.a. p_ext)

        ptlX     ... list of x as meas with FVC and converted to petal coords (a.k.a. X_PTL)
        ptlY     ... list of y as meas with FVC and converted to petal coords (a.k.a. Y_PTL)

        gearT    ... list of GEAR_CALIB_T at the time the robot made the move
        gearP    ... list of GEAR_CALIB_P at the time the robot made the move

        mode     ... 'static' or 'dynamic'
                     In static mode, all "dynamic" parameters will be held fixed,
                     no matter the value of the keep_fixed argument. And vice-
                     versa when in dynamic mode.

        nominals ... dict with keys = param name and vals = initial values for
                     each parameter

        bounds   ... dict with keys = param name and vals = (min,max) limits to
                     use in the best-fit search

        keep_fixed  ... list of any parameters you want forced to their nominal
                        values, and then never varied, when doing the best-fit.

    OUTPUTS:
        best_params ... dict of best-fit results, keys = param names
        final_err   ... numeric error of the best-fit params
    '''
    if nominals is None:
        nominals = default_values
    if bounds is None:
        bounds = default_bounds
    if keep_fixed is None:
        keep_fixed = []
    # arg checking
    assert len(posintT) == len(posintP) == len(ptlX) == len(ptlY) == len(gearT) == len(gearP)
    assert all(isinstance(arg, list) for arg in (posintT, posintP, ptlX, ptlY, gearT, gearP))
    assert all(math.isfinite(x) for x in posintT + posintP + ptlX + ptlY + gearT + gearP)
    assert mode in {'static', 'dynamic'}
    assert all(key in nominals for key in all_keys)
    assert all(key in bounds for key in all_keys)
    assert all(key in all_keys for key in keep_fixed)

    # selection of which parameters are variable
    if mode == 'static':
        keep_fixed = set(keep_fixed).union(dynamic_keys)
    else:
        keep_fixed = set(keep_fixed).union(static_keys)
    params_to_fit = set(nominals).difference(keep_fixed)
    params_to_fit = sorted(params_to_fit) # because optimizer expects a vector with consistent ordering
    param_idx = {key:params_to_fit.index(key) for key in params_to_fit}

    # some re-names and pre-initializations of numpy arrays
    t_int = np.array(posintT)
    p_int = np.array(posintP)
    x_ptl = np.array(ptlX)
    y_ptl = np.array(ptlY)

    # pre-calculate any one-time coordinate transforms
    n_pts = len(posintT)
    x_flat, y_flat = pos2ptl.ptl2flat(x_ptl, y_ptl)
    if mode == 'dynamic':
        n_pts -= 1 # since in dynamic mode we are using delta values

        # Note how measured_int0 depends on having been fed good "static" param
        # values. Also note extraction of only theta and phi from the tuple
        # returned by loc2int. The 2nd element is the "unreachable" boolean
        # flags, ignored here.
        measured_x_loc_0 = pos2ptl.flat2loc(x_flat, nominals['OFFSET_X'])
        measured_y_loc_0 = pos2ptl.flat2loc(y_flat, nominals['OFFSET_Y'])
        measured_t_int_0, measured_p_int_0, unreachable = pos2ptl.loc2int(
            measured_x_loc_0, measured_y_loc_0,
            nominals['LENGTH_R1'], nominals['LENGTH_R2'],
            nominals['OFFSET_T'], nominals['OFFSET_P']
            )
        del unreachable
        dt_int = pos2ptl.delta_angle(u_start=t_int[:-1], u_final=t_int[1:], direction=0)
        dp_int = pos2ptl.delta_angle(u_start=p_int[:-1], u_final=p_int[1:], direction=0)

        # compensate for any gear calibrations in place at time of move
        # this is somewhat aesthetic, so that SCALE_* correponds 1:1 with GEAR_CALIB_*
        # e.g. with a gearT==0.5 (and let's suppose that calibration is physically
        # accurate), then the operation below will make it look to the fitter
        # like dt_int is 2x bigger than measured distance. hence fitter will
        # output SCALE_T --> 0.5
        dt_int /= gearT[1:]
        dp_int /= gearP[1:]

        # now remove first point, to match up with lists of deltas
        x_flat = x_flat[1:]
        y_flat = y_flat[1:]
        measured_t_int_0 = measured_t_int_0[:-1]
        measured_p_int_0 = measured_p_int_0[:-1]

    # set up the consumer function for the variable parameters vector
    p0 = nominals.copy()
    def expected_xy(params):
        if mode =='static':
            for key, idx in param_idx.items():
                p0[key] = params[idx]
            t_int_expected = t_int
            p_int_expected = p_int
        else:
            dt_int_scaled = dt_int * params[param_idx['SCALE_T']]
            dp_int_scaled = dp_int * params[param_idx['SCALE_P']]
            t_int_expected = measured_t_int_0 + dt_int_scaled
            p_int_expected = measured_p_int_0 + dp_int_scaled
        x_loc, y_loc = pos2ptl.int2loc(t_int=t_int_expected,
                                       p_int=p_int_expected,
                                       r1=p0['LENGTH_R1'],
                                       r2=p0['LENGTH_R2'],
                                       t_offset=p0['OFFSET_T'],
                                       p_offset=p0['OFFSET_P'])
        x_flat_expected = pos2ptl.loc2flat(x_loc, p0['OFFSET_X'])
        y_flat_expected = pos2ptl.loc2flat(y_loc, p0['OFFSET_Y'])
        return x_flat_expected, y_flat_expected

    # define error function and run optimizer
    def err_norm(params):
        x_exp, y_exp = expected_xy(params)
        err_x = x_exp - x_flat
        err_y = y_exp - y_flat
        sumsq = np.sum(err_x**2 + err_y**2)
        return (sumsq / n_pts)**0.5

    if mode=='static' :
        # we initialize OFFSET_X and OFFSET_Y by fitting circles
        # for all unique values of p_int
        unique_p_int = np.unique(p_int)
        offset_x=list()
        offset_y=list()
        for val in unique_p_int :
            selection = np.where(p_int==val)[0]
            if selection.size<3 : continue # no circle to fit here
            xc,yc,_,_ = robust_fit_circle(x_flat[selection],y_flat[selection])
            offset_x.append(xc)
            offset_y.append(yc)
        if len(offset_x)>0 :
            offset_x = np.median(offset_x)
            offset_y = np.median(offset_y)
            if not "OFFSET_X" in keep_fixed : nominals["OFFSET_X"] = offset_x
            if not "OFFSET_Y" in keep_fixed : nominals["OFFSET_Y"] = offset_y
            print("init OFFSET_X={} OFFSET_Y={}".format(offset_x,offset_y))

    initial_params = [nominals[key] for key in params_to_fit]
    bounds_vector = [bounds[key] for key in params_to_fit]
    optimizer_result = scipy.optimize.minimize(fun=err_norm,
                                               x0=initial_params,
                                               bounds=bounds_vector)

    # organize and return results
    best_params = {key: optimizer_result.x[param_idx[key]] for key in params_to_fit}
    fixed_params = {key: nominals[key] for key in keep_fixed}
    best_params.update(fixed_params)
    best_params['OFFSET_T'] = wrap_at_180(best_params['OFFSET_T'])
    err = err_norm(optimizer_result.x)
    return best_params, err

def wrap_at_180(angle):
    angle %= 360
    if angle > 180:
        angle -= 360
    return angle
