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

def run_best_fits(table, datum_dates, data_window, savedir, static_quantile):
    '''Define best-fit analysis cases for one positioner.
    
    INPUTS:
        table ... astropy table of measured vs commanded data, for just one positioner
        datum_dates ... keypoint dates at which to evaluate a window of data
        data_window ... mininum number of data points per best-fit
        savedir ... directory to save result files
        static_quantile ... least-error group of static params to median --> overall best
        
    OUTPUTS:
        logstr ... string describing (very briefly) what was done. suitable for
                   printing to stdout, etc
    '''
    posid = table['POS_ID'][0]
    table.sort('DATE_SEC')  # for speed, this pre-sort is ASSUMED by some sub-functions
    cases = _define_cases(table, datum_dates, data_window)
    
    # FIRST-PASS:  STATIC PARAMETERS
    static_out = _process_cases(table, cases, mode='static',
                                param_nominals=fitter.default_values.copy())
	
    # DECIDE ON BEST STATIC PARAMS
    best_static = fitter.default_values.copy()
    errors = static_out['FIT_ERROR_STATIC']
    quantile = np.percentile(errors, static_quantile * 100)
    selection = static_out[errors <= quantile]
    these_best = {key:np.median(selection[key + _mode_suffix('static')]) for key in best_static}
    best_static.update(these_best)
    print(f'{posid}: Selected best static params = {_dict_str(best_static)}')
	
    # SECOND-PASS:  DYNAMIC PARAMETERS
    dynamic_out = _process_cases(table, cases, mode='dynamic',
                                 param_nominals=best_static)

    # MERGED STATIC + DYNAMIC
    #static_out.sort(['POS_ID','DATA_END_DATE'])
    #dynamic_out.sort(['POS_ID','DATA_END_DATE'])
    # why non-unique rows in output?
    # I think because i lost my dict of cases --> unique cases only logic
    merged = join(static_out, dynamic_out, keys=['POS_ID', 'DATA_END_DATE_SEC'])
    merged.sort(['POS_ID','DATA_END_DATE_SEC'])
    merged_filepath = os.path.join(savedir, posid + '_paramfits.csv')
    merged_filepath = os.path.realpath(merged_filepath)
    merged.write(merged_filepath, overwrite=True)
    logstr = f'{posid}: Best-fit results on {len(cases)} data windows written to {merged_filepath}'
    return logstr

def _define_cases(table, datum_dates, data_window):
    '''Define best-fit analysis cases for one positioner.
    
    INPUTS:
        table ... astropy table of measured vs commanded data, for just one positioner
        datum_dates ... keypoint dates at which to evaluate a window of data
        data_window ... mininum number of data points per best-fit
        
    OUTPUT:
        cases ... list of dicts defining the analysis cases, with keys:
                    'start_idx' ... first data index in analyis window
                    'final_idx' ... last data index in analysis window
    '''
    cases = []
    widths = []
    start_idxs = []
    final_idxs = []
    for j in range(1, len(datum_dates)):
        start_date = datum_dates[j - 1]
        final_date = datum_dates[j]
        start_idxs.append(_row_idx_for_time(table, start_date))
        final_idxs.append(_row_idx_for_time(table, final_date))
        widths.append(final_idxs[-1] - start_idxs[-1])
    for J in range(len(final_idxs)):
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
                'final_idx': final_idxs[J]}
        cases.append(case)
    posid = table['POS_ID'][0]
    print(f'{posid}: {len(cases):5d} analysis cases defined')
    return cases

def _process_cases(table, cases, mode, param_nominals):
    '''Feed analysis cases for a positioner through the best-fit function.
    
    INPUTS:
        table ... Astropy table of measured vs commanded data, for just one positioner
        cases ... List of analysis cases as generated by _define_cases()
        mode ... 'static' or 'dynamic'
        param_nominals ... Dict of nominal parameter values, structured like fitter.default_values

    OUTPUT:
        output ... Astropy table of results, headers given by output_keys.
                   However any key 'X' with output_key['X'] == True gets a
                   suffix for the mode ('_STATIC' or '_DYNAMIC') appended.
    '''
    posid = table['POS_ID'][0]
    output = {col:[] for col in output_keys}
    for case in cases:
        m = case['start_idx']
        n = case['final_idx']                   
        xytp_data, subtable = _select_by_index(table, start=m, final=n+1)
        print(f'{posid}: fitting {mode.upper()} params for over data period:')
        print(f'  start idx = {m:5d}, date = {subtable["DATE"][0]}')
        print(f'  final idx = {n:5d}, date = {subtable["DATE"][-1]}')
        print(f'  num points = {n-m+1:5d}')
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
        print(f'{posid}: best params = {_dict_str(params)}')
        print(f'  fit error = {err:.3f}')
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