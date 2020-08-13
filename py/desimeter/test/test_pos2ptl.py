import unittest
import json
from pkg_resources import resource_filename
import os
import numpy as np
import desimeter.transform.pos2ptl as pos2ptl
import desimeter.transform.xy2tp as xy2tp
        
class TestPos2Ptl(unittest.TestCase):
    '''Tests below compare results from this module to a presumed canonical
    set of coordinate conversions using postransforms.py. See utility script
    desimeter/bin/generate_pos2ptl_testdata, for code to generate the data set.
    '''
    def test_all(self):
        tol = 0.005
        desimeter_data_dir = resource_filename("desimeter", "data")
        filename = 'pos2ptl_test_data.json'
        path = os.path.join(desimeter_data_dir, filename)
        with open(path, 'r') as file:
            u = json.load(file)
                                                                                 
        tests = {'ptl2flat': {'func': pos2ptl.ptl2flat, 'inputs': ['x_ptl', 'y_ptl'], 'checks': ['x_flat', 'y_flat']},
                 'flat2ptl': {'func': pos2ptl.flat2ptl, 'inputs': ['x_flat', 'y_flat'], 'checks': ['x_ptl', 'y_ptl']},
                 'flat2loc_x': {'func': pos2ptl.flat2loc, 'inputs': ['x_flat', 'OFFSET_X'], 'checks': ['x_loc']},
                 'flat2loc_y': {'func': pos2ptl.flat2loc, 'inputs': ['y_flat', 'OFFSET_Y'], 'checks': ['y_loc']},
                 'loc2flat_x': {'func': pos2ptl.loc2flat, 'inputs': ['x_loc', 'OFFSET_X'], 'checks': ['x_flat']},
                 'loc2flat_y': {'func': pos2ptl.loc2flat, 'inputs': ['y_loc', 'OFFSET_Y'], 'checks': ['y_flat']},
                 'ext2loc': {'func': pos2ptl.ext2loc, 'inputs': ['t_ext', 'p_ext', 'LENGTH_R1', 'LENGTH_R2'], 'checks': ['x_loc', 'y_loc']},
                 'loc2ext': {'func': pos2ptl.loc2ext, 'inputs': ['x_loc', 'y_loc', 'LENGTH_R1', 'LENGTH_R2', 'OFFSET_T', 'OFFSET_P'], 'checks': ['t_ext', 'p_ext']},
                 'ext2int_t': {'func': pos2ptl.ext2int, 'inputs': ['t_ext', 'OFFSET_T'], 'checks': ['t_int']},
                 'ext2int_p': {'func': pos2ptl.ext2int, 'inputs': ['p_ext', 'OFFSET_P'], 'checks': ['p_int']},
                 'int2ext_t': {'func': pos2ptl.int2ext, 'inputs': ['t_int', 'OFFSET_T'], 'checks': ['t_ext']},
                 'int2ext_p': {'func': pos2ptl.int2ext, 'inputs': ['p_int', 'OFFSET_P'], 'checks': ['p_ext']},
                 'int2loc': {'func': pos2ptl.int2loc, 'inputs': ['t_int', 'p_int', 'LENGTH_R1', 'LENGTH_R2', 'OFFSET_T', 'OFFSET_P'], 'checks': ['x_loc', 'y_loc']},
                 'loc2int': {'func': pos2ptl.loc2int, 'inputs': ['x_loc', 'y_loc', 'LENGTH_R1', 'LENGTH_R2', 'OFFSET_T', 'OFFSET_P'], 'checks': ['t_int', 'p_int']},
                 }
        
        # test of delta_angle() function defined separately
        # (there is no exact matching function in postransforms.py)
        delta = {'u_start':   [15, 15,   15, 0,  90,   90,  90,   0,   0,    0,   0,   0,    0,   0,   0,    0,   0,   0,    0,  800,  800,  800,    0,  720,  720],
                 'u_final':   [45, 45,   45, 0,  45,   45,  45, 180, 180,  180, 360, 360,  360, 400, 400,  400, 800, 800,  800,    0,    0,    0,  720,    0,    0],
                 'direction': [ 0, +1,   -1, 0,   0,   +1,  -1,   0,  +1,   -1,   0,  +1,   -1,   0,  +1,   -1,   0,  +1,   -1,    0,   +1,   -1,   -1,   -1,   +1],
                 'delta_u':   [30, 30, -330, 0, -45, +315, -45, 180, 180, -180, 360, 360, -360, 400, 400, -320, 800, 800, -280, -800, +280, -800, -720, -720, +720]}
        u.update(delta)
        delta_test = {'delta_angle': {'func':pos2ptl.delta_angle,
                                      'inputs': ['u_start', 'u_final', 'direction'],
                                      'checks':['delta_u']}}
        tests.update(delta_test)
        all_out = {}
        for name, args in tests.items():
            inputs = [u[key] for key in args['inputs']]
            checks = [u[key] for key in args['checks']]
            if name in {'ext2loc', 'int2loc'}:
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
                assert errors[i][worst] < tol

    def test_xy2tp(self):
        tol = 0.001
        r_test = [2.5, 3.5]
        t_test = [-190, 0, 190]
        p_test = [-10, 10, 170, 190]
        ranges = [[-200, 200], [-20, 200]]
        t_guess_err = 3
        t_guess_tol = 10
        R = [[r1, r2] for r1 in r_test for r2 in r_test]
        tp_test = [[t, p] for t in t_test for p in p_test]
        for r in R:
            for tp in tp_test:
                for sign in [-1, 0, +1]:
                    xy_test = xy2tp.tp2xy(tp, r)
                    t_guess = tp[0] + sign*t_guess_err
                    tp_calc, unreachable = xy2tp.xy2tp(xy_test, r, ranges, t_guess, t_guess_tol)
                    assert not unreachable
                    for i in [0, 1]:
                        error = abs(tp_calc[i] - tp[i])
                        if error >= tol:
                            print(f'\nError in xy2tp test:\n error[{i}] = {error:.4f}\n tp_calc = {tp_calc}' +
                                  f'\n tp_nom = {tp}\n r = {r}\n t_guess = {t_guess}')
                            print(f' Sample call:\n  xy2tp.xy2tp({xy_test},{r},{ranges},{t_guess},{t_guess_tol})')
                        assert error < tol 
        
if __name__ == '__main__':
    unittest.main()
    
