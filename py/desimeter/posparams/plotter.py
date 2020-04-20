# -*- coding: utf-8 -*-
"""
Plot time series of positioner parameters.
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
import astropy.stats as stats

# common options
img_ext = '.png'
error_series_filename = 'fiterror' + img_ext
tick_period_days = 7
day_in_sec = 24*60*60
error_bins = [0, 0.01, 0.05, 0.1, np.inf] # for error_series plot, bin edges in which to categorize positioners

# plot definitions for parameter subplots
param_subplot_defs = [
    {'keys': ['FIT_ERROR_STATIC', 'FIT_ERROR_DYNAMIC'],
     'units': 'um RMS',
     'mult': 1000,
     'subplot': 1,
     'logscale': True,
     'equal_scales': True},
            
    {'keys': ['NUM_POINTS'],
     'units': 'USED IN FIT',
     'mult': 1,
     'subplot': 4,
     'logscale': False},
            
    {'keys': ['LENGTH_R1', 'LENGTH_R2'],
     'units': 'mm',
     'mult': 1,
     'subplot': 2,
     'logscale': False,
     'equal_scales': True},
            
    {'keys': ['OFFSET_X', 'OFFSET_Y'],
     'units': 'mm',
     'mult': 1,
     'subplot': 5,
     'logscale': False,
     'equal_scales': False},
            
    {'keys': ['OFFSET_T', 'OFFSET_P'],
     'units': 'deg',
     'mult': 1,
     'subplot': 6,
     'logscale': False,
     'equal_scales': False},
            
    {'keys': ['SCALE_T', 'SCALE_P'],
     'units': '',
     'mult': 1,
     'subplot': 3,
     'logscale': False,
     'equal_scales': True},
    ]

def plot_params(table, savepath, statics_during_dynamic):
    '''Plot time series of positioner parameters for a single positioner.
    
    Inputs:
        table    ... Astropy table as generated by fit_params, then reduced to
                     just the rows for a single POS_ID.
        
        savepath ... Where to save output plot file. Extension determines image
                     format.
        
        statics_during_dynamic ... Dict of static params used during the
                     dynamic params best-fit
        
    Outputs:
        The plot image file is saved to savepath. Also the savepath string is
        returned unaltered (for ease of tracking progress during multi-
        processing).
    '''
    fig = _init_plot()
    posid = table['POS_ID'][0]
    fig.subplots_adjust(wspace=.3, hspace=.3)
    times = table['DATE_SEC']
    tick_values = np.arange(times[0], times[-1]+day_in_sec, tick_period_days*day_in_sec)
    tick_labels = [Time(t, format='unix', out_subfmt='date').iso for t in tick_values]
    n_pts = len(table)
    marker = ''
    for p in param_subplot_defs:
        plt.subplot(2, 3, p['subplot'])
        for key in p['keys']:
            ax_right = None
            if p['keys'].index(key) == 1:
                ax_right = plt.twinx()
                color = 'red'
                linestyle = '--'
                if n_pts == 1:
                    marker = '^'
            else:
                color = 'blue'
                linestyle = '-'
                if n_pts == 1:
                    marker = 'v'
            y = [val * p['mult'] for val in table[key]]#.tolist()]
            plt.plot(times, y, color=color, linestyle=linestyle, marker=marker)
            if not ax_right:
                ax_left = plt.gca()
            units = f' ({p["units"]})' if p['units'] else ''
            plt.ylabel(key + units, color=color)
            if p['logscale']:
                plt.yscale('log')
            if 'ylims' in p:
                plt.ylim(p['ylims'])
            if ax_right and p['equal_scales']:
                min_y = min(ax_left.get_ylim()[0], ax_right.get_ylim()[0])
                max_y = max(ax_left.get_ylim()[1], ax_right.get_ylim()[1])
                ax_left.set_ylim((min_y, max_y))
                ax_right.set_ylim((min_y, max_y))
            plt.xticks(tick_values, tick_labels, rotation=90, horizontalalignment='center', fontsize=8)
            plt.yticks(fontsize=8)
            if key == 'SCALE_P':
                s = statics_during_dynamic
                plt.text(min(plt.xlim()), min(plt.ylim()),
                         f' Using static params:\n'
                         f' LENGTH_R1 = {s["LENGTH_R1"]:>5.3f}, LENGTH_R2 = {s["LENGTH_R2"]:>5.3f}\n'
                         f' OFFSET_X = {s["OFFSET_X"]:>8.3f}, OFFSET_Y = {s["OFFSET_Y"]:>8.3f}\n'
                         f' OFFSET_T = {s["OFFSET_T"]:>8.3f}, OFFSET_P = {s["OFFSET_P"]:>8.3f}\n',
                         verticalalignment='bottom')
    analysis_date = table['ANALYSIS_DATE_DYNAMIC'][-1]
    title = f'{posid}'
    title += f'\nbest-fits to historical data'
    title += f'\nanalysis date: {analysis_date}'
    plt.suptitle(title)
    _save_and_close_plot(fig, savepath)
    return savepath

def plot_errors(table, savepath):
    '''Plot time series of posparam fit errors, summing up at each time point
    for all postioners represented in table.
    
    Inputs:
        table    ... Astropy table as generated by fit_params, containing data
                     rows for multiple positioners.
        
        savepath ... Where to save output plot file. Extension determines image
                     format.
        
    Outputs:
        The plot image file is saved to savepath.
    '''
    fig = _init_plot()

    # This is very much in-progress code. Idea is:
    #
    #   for any pos_id with multiple rows on the same day:
    #       select row with lowest FIT_ERROR value and delete the others
    #   for each threshold level: (like "ceiling" in proto code below)
    #       select just the rows within that threshold zone
    #       set up daily bins
    #       histogram the dates-in-seconds into those bins
    #
    #   This yields a count of how many positioners were measured during each
    #   day as meeting the threshold conditions.
    #
    #   Finally, cumulatively sum the bins.
    #
    #   Merits: easy / quick to write
    #   Demerits: awkward to capture positioners failing out of the cumulative
    #             counts, because losing track of posids.
    #
    # Alternate method:
    #   
    #   A bit more ploddingly, make a dict carrying an empty set for each day.
    #   These are the bins. Now just march through the goddamned table, chucking
    #   posids into these sets if they meet the threshold criteria.
    #
    #   Merits: keeps track of POS_IDs
    #   Demerits: not fancy (probably a merit)
    ceiling = 0.02
    time_min = min(table['DATE_SEC'])
    time_max = max(table['DATE_SEC'])
    nbins = int(np.ceil((time_max - time_min) / day_in_sec))
    time_range = (time_min, time_min + nbins * day_in_sec)
    d = table[table['FIT_ERROR_DYNAMIC'] <= ceiling]
    h = stats.histogram(d['DATE_SEC'], bins=nbins, range=time_range)
        
    _save_and_close_plot(fig, savepath)
    

def _init_plot():
    '''Internal common plot initialization function. Returns figure handle.'''
    plt.ioff()
    fig = plt.figure(figsize=(20,10), dpi=150)
    plt.clf()
    return fig
    
def _save_and_close_plot(fig, savepath):
    '''Internal common plot saving and closing function. Argue the figure
    handle to close, and the path where to save the image. Extension determines
    image format.'''
    plt.savefig(savepath, bbox_inches='tight')
    plt.close(fig)    
