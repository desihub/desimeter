# -*- coding: utf-8 -*-
"""
Wrapper module to handle running the best-fit cases. In general this module
is equipped to ingest a data history for a given positioner, do the best fits,
and then save the results to csv report file.
"""

import os
import numpy as np
from astropy.time import Time
from astropy.table import Table
from astropy.table import join

# imports below require <path to desimeter>/py' to be added to system PYTHONPATH.
import desimeter.posparams.fitter as fitter
import desimeter.transform.ptl2fp as ptl2fp
from desimeter.posparams.posmoveselection import posmove_selection

# column headers for processed results of best-fits
# boolean says whether it has "STATIC"/"DYNAMIC" options
output_keys = {'POS_ID': False,
               'NUM_POINTS': False,
               'DATA_START_DATE': False,
               'DATA_END_DATE': False,
               'DATA_START_DATE_SEC': False,
               'DATA_END_DATE_SEC': False,
               'ANALYSIS_DATE': True,
               'FIT_ERROR': True,
               }
output_keys.update({key:True for key in fitter.all_keys})

def run_best_fits(posid, path, period_days, data_window, savedir,
                  static_quantile, printf=print, min_window=10, log_note_selection=None):
    '''Define best-fit analysis cases for one positioner.

    INPUTS:
        posid ... positioner id string
        path ... file path to csv containing measured (x,y) vs internal (t,p) data
        period_days ... minimum number of days per best-fit
        data_window ... number of data points per best-fit
        savedir ... directory to save result files
        static_quantile ... least-error group of static params to median --> overall best
        printf ... print function (so you can control how this module spits any messages)
        min_window ... (optional) minimum number of data points to attempt a fit
        log_note_selection ... (optional) required keywords in LOG_NOTE selection, seperated by '&', like ''arc calibration & use_desimeter=True'

    OUTPUTS:
        logstr ... string describing (very briefly) what was done. suitable for
                   printing to stdout, etc
    '''
    table = _read(posid=posid, path=path, printf=printf, log_note_selection=log_note_selection)
    if not table:
        return f'{posid}: dropped from analysis (invalid table)'
    _make_petal_xy(table, printf=printf)
    _apply_filters(table, printf=printf)
    if not table or len(table) < min_window:
        return f'ERROR: {posid} dropped from analysis because {len(table)}<{min_window} is insufficient data remaining after filters applied'
    _make_sec_from_epoch(table, printf=printf)
    period_sec = period_days * 24 * 60 * 60
    datum_dates = np.arange(min(table['DATE_SEC']), max(table['DATE_SEC']), period_sec)
    datum_dates = datum_dates.tolist()
    cases = _define_cases(table, datum_dates, data_window, printf=printf)

    # FIRST-PASS:  STATIC PARAMETERS
    static_out = _process_cases(table, cases, printf=printf, mode='static',
                                param_nominals=fitter.default_values.copy(),
                                )

    # DECIDE ON BEST STATIC PARAMS
    best_static = fitter.default_values.copy()
    errors = static_out['FIT_ERROR_STATIC']
    quantile = np.percentile(errors, static_quantile * 100)
    selection = static_out[errors <= quantile]
    these_best = {key:np.median(selection[key + _mode_suffix('static')]) for key in best_static}
    best_static.update(these_best)
    printf(f'{posid}: Selected best static params = {_dict_str(best_static)}')

    # SECOND-PASS:  DYNAMIC PARAMETERS
    dynamic_out = _process_cases(table, cases, printf=printf, mode='dynamic',
                                 param_nominals=best_static)

    # MERGED STATIC + DYNAMIC
    same_keys = [key for key, is_variable in output_keys.items() if not is_variable]
    merged = join(static_out, dynamic_out, keys=same_keys)
    merged.sort(['POS_ID','DATA_END_DATE_SEC'])
    merged_filepath = os.path.join(savedir, posid + '_paramfits.csv')
    merged_filepath = os.path.realpath(merged_filepath)
    merged.write(merged_filepath, overwrite=True)
    logstr = f'{posid}: Best-fit results on {len(cases)} data windows written to {merged_filepath}'
    return logstr

def _read(posid, path, printf=print, log_note_selection=None):
    '''Read csv file at path. File should contain positioner measured vs
    commanded data. Data is appropriately type cast as necessary.

    INPUTS:  posid ... expected unique POS_ID string to be found within table
             path ... file path to csv file for one positioner's data
             printf ... print function

    OUTPUT:  table ... astropy table or None if no valid table found
    '''
    required_keys = {'DATE', 'POS_P', 'POS_T', 'X_FP', 'Y_FP', 'CTRL_ENABLED',
                     'TOTAL_MOVE_SEQUENCES', 'GEAR_CALIB_T', 'GEAR_CALIB_P',
                     'POS_ID', 'PETAL_LOC'}
    typecast_funcs = {'CTRL_ENABLED': _str2bool}
    table = Table.read(path)


    table = posmove_selection(table,log_note_selection)

    table.sort('DATE')
    key_search = [key in table.columns for key in required_keys]
    if table and all(key_search) and str(posid) == _get_posid(table):
        printf(f'{_get_posid(table)}: Read from {path}')
        for key,func in typecast_funcs.items():
            table[key] = [func(val) for val in table[key]]
        return table
    return None

def _make_petal_xy(table, printf=print):
    '''Produces new columns for X_PTL and Y_PTL coordinates of measured (x,y)
    data. Input is an astropy table, which is altered in-place.'''
    ptlXYZ = ptl2fp.fp2ptl(petal_loc=table['PETAL_LOC'][0],
                           x_fp=table['X_FP'],
                           y_fp=table['Y_FP'],
                           z_fp=None)
    table['X_PTL'] = ptlXYZ[0]
    table['Y_PTL'] = ptlXYZ[1]
    printf(f'{_get_posid(table)}: (X_FP, Y_FP) converted to (X_PTL, Y_PTL)')

# Data filtering functions
# Typically applied using the _apply_filters() function below
def _ctrl_not_enabled(table):
    '''Returns array of bools whether control was enabled in each row of table.'''
    return table['CTRL_ENABLED'] == False

def _no_move_performed(table):
    '''Returns array of bools whether a move was performed in each row of table.
    Always assumed True for first row.'''
    table.sort('TOTAL_MOVE_SEQUENCES')
    dummy = table['TOTAL_MOVE_SEQUENCES'][0] - 1
    moveseq_diff = _np_diff_with_prepend(table['TOTAL_MOVE_SEQUENCES'], prepend=dummy)
    return moveseq_diff == 0

def _no_cruise_performed(table):
    '''Returns array of bools whether a cruise move was performed in each row
    of table. Always assumed True for first row.'''
    table.sort('TOTAL_MOVE_SEQUENCES')
    dummyT = table['TOTAL_MOVE_SEQUENCES'][0] - 1
    dummyP = table['TOTAL_MOVE_SEQUENCES'][0] - 1
    cruiseT_diff = _np_diff_with_prepend(table['TOTAL_CRUISE_MOVES_T'], prepend=dummyT)
    cruiseP_diff = _np_diff_with_prepend(table['TOTAL_CRUISE_MOVES_P'], prepend=dummyP)
    return cruiseT_diff + cruiseP_diff == 0 # these diffs are always >= 0 for sorted table input

def _apply_filters(table, bad_row_funcs=None, printf=print):
    '''Removes bad rows from table, according to filter functions specified
    in the list bad_row_funcs. Note as of 2020-04-16, default bad_row_funcs
    intentionally deoes not include no_cruise_performed().
    '''
    if bad_row_funcs is None:
        bad_row_funcs = [_ctrl_not_enabled, _no_move_performed]
    initial_len = len(table)
    for func in bad_row_funcs:
        bad_rows = func(table)
        table.remove_rows(bad_rows)
        if len(table) == 0:
            break
    final_len = len(table)
    printf(f'{_get_posid(table)}: dropped {initial_len - final_len} of {initial_len} non-conforming data rows')

def _make_sec_from_epoch(table, printf=print):
    '''Produces new column for seconds-since-epoch version of date. Input is
    an astropy table containing column 'DATE'. The table is altered in-place.'''
    dates = table['DATE']
    table['DATE_SEC'] = Time(dates, format='iso').unix
    printf(f'{_get_posid(table)}: generated seconds-since-epoch column (\'DATE_SEC\')')

def _define_cases(table, datum_dates, data_window, printf=print):
    '''Define best-fit analysis cases for one positioner.

    INPUTS:
        table ... astropy table of measured vs commanded data, for just one positioner
        datum_dates ... keypoint dates at which to evaluate a window of data
        data_window ... mininum number of data points per best-fit
        printf ... print function

    OUTPUT:
        cases ... list of dicts defining the analysis cases, with keys:
                    'start_idx' ... first data index in analyis window
                    'final_idx' ... last data index in analysis window
    '''
    table.sort('DATE_SEC')  # for speed, _row_idx_for_time() ASSUMEs this pre-sort has been done
    cases = []
    widths = []
    start_idxs = []
    final_idxs = []
    posid = _get_posid(table)

    #- Catch special cases of single day and/or short input table
    num_moves = len(table)
    if len(datum_dates) == 1:
        for istart in range(0, num_moves, data_window):
            ifinal = min(len(table)-1, istart+data_window-1)
            #- check if window is smaller than requested
            if (ifinal-istart < data_window):
                if istart == 0:
                    #- even one window is smaller than requested
                    print(f'WARNING: total num_moves {num_moves} < data_window {data_window}')
                else:
                    #- pad backwards, overlapping with previous window
                    istart = ifinal - data_window

            cases.append(dict(start_idx=0, final_idx=ifinal))

        printf(f'{posid}: {len(cases):5d} analysis cases defined')
        return cases

    for j in range(1, len(datum_dates)):
        start_date = datum_dates[j - 1]
        final_date = datum_dates[j]
        start_idxs.append(_row_idx_for_time(table, start_date))
        final_idxs.append(_row_idx_for_time(table, final_date))
        widths.append(final_idxs[-1] - start_idxs[-1])
    for J,final_idx in enumerate(final_idxs):
        backwards_sum = lambda i: sum(w for w in widths[i:J+1])
        satisfies_min_width = lambda i: backwards_sum(i) > data_window
        should_expand_backwards = lambda i: not(satisfies_min_width(i)) and i > 0
        I = J
        while should_expand_backwards(I):
            I -= 1
        can_expand_forwards = J < len(widths) - 1
        if not satisfies_min_width(I) and can_expand_forwards:
            continue  # skip ahead and try the next datum date
        case = {'start_idx': start_idxs[I],
                'final_idx': final_idx}
        if case not in cases:
            cases.append(case)

    printf(f'{posid}: {len(cases):5d} analysis cases defined')
    return cases

def _process_cases(table, cases, mode, param_nominals, printf=print):
    '''Feed analysis cases for a positioner through the best-fit function.

    INPUTS:
        table ... Astropy table of measured vs commanded data, for just one positioner
        cases ... List of analysis cases as generated by _define_cases()
        mode ... 'static' or 'dynamic'
        param_nominals ... Dict of nominal parameter values, structured like fitter.default_values
        printf ... print function

    OUTPUT:
        output ... Astropy table of results, headers given by output_keys.
                   However any key 'X' with output_key['X'] == True gets a
                   suffix for the mode ('_STATIC' or '_DYNAMIC') appended.
    '''
    table.sort('DATE_SEC') # to ensure consistent ordering as in _define_cases
    posid = _get_posid(table)
    output = {col:[] for col in output_keys}
    for case in cases:
        m = case['start_idx']
        n = case['final_idx']
        xytp_data, subtable = _select_by_index(table, start=m, final=n+1)
        printf(f'{posid}: fitting {mode.upper()} params over data period:')
        printf(f'  start idx = {m:5d}, date = {subtable["DATE"][0]}')
        printf(f'  final idx = {n:5d}, date = {subtable["DATE"][-1]}')
        printf(f'  num points = {n-m+1:5d}')
        params, err = fitter.fit_params(posintT=xytp_data['posintT'],
                                        posintP=xytp_data['posintP'],
                                        ptlX=xytp_data['ptlX'],
                                        ptlY=xytp_data['ptlY'],
                                        gearT=xytp_data['gearT'],
                                        gearP=xytp_data['gearP'],
                                        mode=mode,
                                        nominals=param_nominals,
                                        bounds=fitter.default_bounds,
                                        keep_fixed=[],
                                        )

        output['ANALYSIS_DATE'].append(Time.now().iso)
        output['POS_ID'].append(posid)
        for suffix in {'', '_SEC'}:
            d = f'DATE{suffix}'
            output[f'DATA_START_{d}'].append(table[d][m])
            output[f'DATA_END_{d}'].append(table[d][n])
        output['NUM_POINTS'].append(n - m + 1)
        output['FIT_ERROR'].append(err)
        for key in fitter.static_keys:
            output[key].append(params[key])
        for key in fitter.dynamic_keys:
            output[key].append(params[key])
        printf(f'{posid}: best params = {_dict_str(params)}')
        printf(f'  fit error = {err:.3f}')
    for key in list(output):
        if not output[key]:
            output[key] = [None]*len(cases)
        if output_keys[key]:
            output[key + _mode_suffix(mode)] = output.pop(key)
    return Table(output)

def _row_idx_for_time(presorted_table, t):
    '''Returns index of first row in table which matches the argued time t.
    Value t should be in same scale as 'DATE_SEC' column. I.e. the date in
    seconds-since-epoch. If t is outside the time range of the table, the
    function will return either index 0 or max index.

    Note: For speed, this function assumes the provided table has already been
    pre-sorted, ascending by time!
    '''
    if t < presorted_table['DATE_SEC'][0]:
        return 0
    if t > presorted_table['DATE_SEC'][-1]:
        return len(presorted_table) - 1
    return int(np.argwhere(presorted_table['DATE_SEC'] >= t)[0])

def _select_by_index(table, start=0, final=-1):
    '''Returns a subset of data formatted for the best-fitting function.'''
    data = {}
    subtable = table[start:final]
    data['posintT'] = subtable['POS_T'].tolist()
    data['posintP'] = subtable['POS_P'].tolist()
    data['ptlX'] = subtable['X_PTL'].tolist()
    data['ptlY'] = subtable['Y_PTL'].tolist()
    data['gearT'] = subtable['GEAR_CALIB_T'].tolist()
    data['gearP'] = subtable['GEAR_CALIB_P'].tolist()
    return data, subtable

def _np_diff_with_prepend(array, prepend):
    '''For numpy users v1.16 and earlier, where no prepend arg available in the
    diff function. Can perhaps be replaced in future by direct usage of
    np.diff(array, prepend=prepend). This wrapper function only works for scalar
    values of prepend.'''
    prepend_array = np.array([prepend])
    prepended = np.concatenate((prepend_array, np.array(array)))
    return np.diff(prepended)

def _dict_str(d):
    '''Return a string displaying formatted values in dict d.'''
    s = '{'
    for key,val in d.items():
        s += f"'{key}':{val:.3f}, "
    s = s[:-2]
    s += '}'
    return s

def _mode_suffix(mode):
    '''Return a common suffix string for mode.'''
    assert mode in {'static', 'STATIC', 'dynamic', 'DYNAMIC'}
    return '_' + mode.upper()

def _str2bool(s):
    '''Converts common string representations of booleans to bool.'''
    return s in {'True', 'true', 'TRUE', True, 1, 'Yes', 'yes', 'YES'}

def _get_posid(table):
    '''Grabs the positioner id from an astropy table. Returns set if more than
    one or zero posids found in table.'''
    posid = set(table['POS_ID'])
    if len(posid) == 1:
        posid = posid.pop()
    return posid
