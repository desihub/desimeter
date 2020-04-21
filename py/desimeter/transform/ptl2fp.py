"""
Utility functions to fit and apply coordinates transformation from PTL (petal local) to FP (focal plane ~ CS5)
"""

import yaml
import numpy as np
from desimeter.log import get_logger
from desimeter.transform import rszn_lookups
from pkg_resources import resource_filename

petal_alignment_dict = None

# rotation matrices
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

def get_petal_alignment_data() :
    global petal_alignment_dict
    if petal_alignment_dict is None :
        filename = resource_filename('desimeter',"data/petal-alignments.yaml")
        ifile=open(filename)
        petal_alignment_dict = yaml.safe_load(ifile)
        ifile.close()
    return petal_alignment_dict

def apply_ptl2fp(spots) :
    
    log = get_logger()

    petal_alignment_dict = get_petal_alignment_data()
        
    
    nspot = spots['PETAL_LOC'].size

    # local petal coordinates 'PTL'
    xyzptl = np.zeros((3,nspot))
    xyzptl[0] = spots['X_PTL']
    xyzptl[1] = spots['Y_PTL']
    xyzptl[2] = spots['Z_PTL']
    
    # global focal plane coordinates 'FP'
    xyzfp = np.zeros((3,nspot))
    
    for petal in np.unique(spots['PETAL_LOC']) :
        ii = np.where(spots['PETAL_LOC']==petal)[0]
        params = petal_alignment_dict[petal]
        Rotation = Rxyz(params["alpha"],params["beta"],params["gamma"])
        Translation = np.array([params["Tx"],params["Ty"],params["Tz"]])
        xyzfp[:,ii] = Rotation.dot(xyzptl[:,ii]) + Translation[:,None]
    
    spots['X_FP'] = xyzfp[0]
    spots['Y_FP'] = xyzfp[1]
    spots['Z_FP'] = xyzfp[2]

    return spots
    
    
def ptl2fp(petal_loc, xptl, yptl, zptl=None) :
    
    if zptl is None:
        radius = np.hypot(xptl, yptl)
        zptl = rszn_lookups.r2z(radius) # estimate as approx nominal echo22
    xyzptl = np.vstack([xptl,yptl,zptl])

    # global focal plane coordinates 'FP'
    params = get_petal_alignment_data()[petal_loc]
    Rotation = Rxyz(params["alpha"],params["beta"],params["gamma"])
    Translation = np.array([params["Tx"],params["Ty"],params["Tz"]])
    xyzfp = Rotation.dot(xyzptl) + Translation[:,None]
    
    return xyzfp[0],xyzfp[1],xyzfp[2]

