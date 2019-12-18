"""
Utility functions to fit and apply coordinates transformation from PL (petal local) to FP (focal plane ~ CS5)
"""

import numpy as np
from desiutil.log import get_logger
from astropy.table import Table,Column
from pkg_resources import resource_filename


# %% rotation matrices
def Rx(angle):  # all in radians
    Rx = np.array([
        [1.0,           0.0,            0.0],
        [0.0,           np.cos(angle),  -np.sin(angle)],
	[0.0,           np.sin(angle),  np.cos(angle)]
    ])
    return Rx

	
def Ry(angle):  # all in radians
    Ry = np.array([
	[np.cos(angle),  0.0,            np.sin(angle)],
	[0.0,            1.0,            0.0],
	[-np.sin(angle), 0.0,            np.cos(angle)]
    ])
    return Ry
	
	
def Rz(angle):  # all in radians
    Rz = np.array([
	[np.cos(angle), -np.sin(angle), 0.0],
    [np.sin(angle), np.cos(angle),  0.0],
	[0.0,           0.0,            1.0]
    ])
    return Rz

	
def Rxyz(alpha, beta, gamma):  # yaw-pitch-roll system, all in radians
    return Rz(gamma) @ Ry(beta) @ Rx(alpha)  # @ is matrix multiplication
	
def apply_pl2fp(spots,petal_alignment_dict) :

    log = get_logger()
    
    nspot = spots['Petal Loc ID'].size

    # local petal coordinates 'pl'
    xyzpl = np.zeros((3,nspot))

    # rename the columns if needed
    if 'X FCL' in spots.dtype.names :
        spots.rename_column('X FCL', 'XPL')
        log.warning("rename_column('X FCL', 'XPL')")
    if 'Y FCL' in spots.dtype.names :
        spots.rename_column('Y FCL', 'YPL')
        log.warning("rename_column('Y FCL', 'YPL')")
    if 'Z FCL' in spots.dtype.names :
        spots.rename_column('Z FCL', 'ZPL')
        log.warning("rename_column('Z FCL', 'ZPL')")
    
        
    xyzpl[0] = spots['XPL']
    xyzpl[1] = spots['YPL']
    xyzpl[2] = spots['ZPL']
    
    # global focal plane coordinates 'fp'
    xyzfp = np.zeros((3,nspot))
    
    for petal in np.unique(spots['Petal Loc ID']) :
        ii = np.where(spots['Petal Loc ID']==petal)[0]
        params = petal_alignment_dict[petal]
        Rotation = Rxyz(params["alpha"],params["beta"],params["gamma"])
        Translation = np.array([params["Tx"],params["Ty"],params["Tz"]])
        xyzfp[:,ii] = Rotation.dot(xyzpl[:,ii]) + Translation[:,None]
    
    spots['XFP'] = xyzfp[0]
    spots['YFP'] = xyzfp[1]
    spots['ZFP'] = xyzfp[2]

    return spots
    
    
