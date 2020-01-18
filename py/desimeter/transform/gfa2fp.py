"""
Transforms between GFA pixel coordinates and FP mm
"""

import numpy as np
from scipy.optimize import minimize
from desimeter import io
from desimeter.log import get_logger

#- Cached GFA pix -> GFA FP metrology scale, rotation, offsets per petal
_gfa_transforms = None

def gfa2fp(petal_loc, xgfa, ygfa):
    global _gfa_transforms
    if _gfa_transforms is None:
        metrology = io.load_metrology()
        _gfa_transforms = fit_gfa2fp(metrology)
    
    log = get_logger()
    if petal_loc not in _gfa_transforms:
        log.error('PETAL_LOC {} GFA metrology missing'.format(petal_loc))
    
    params = _gfa_transforms[petal_loc]
    xfp, yfp = apply_scale_rotation_offset(xgfa, ygfa, **params)
    
    return xfp, yfp

def fp2gfa(petal_loc, xfp, yfp):
    global _gfa_transforms
    if _gfa_transforms is None:
        metrology = io.load_metrology()
        _gfa_transforms = fit_gfa2fp(metrology)

    log = get_logger()
    if petal_loc not in _gfa_transforms:
        log.error('PETAL_LOC {} GFA metrology missing'.format(petal_loc))
    
    params = _gfa_transforms[petal_loc]
    xgfa, ygfa = undo_scale_rotation_offset(xfp, yfp, **params)
    
    return xgfa, ygfa


#- TODO: fvc2fp has similar transforms; merge that code with this
def apply_scale_rotation_offset(x, y, scale=1.0, rotation=0.0, xoffset=0.0, yoffset=0.0, ):
    """
    Convenience wrapper to return x,y transformed by offset, scale, rotation
    
    Args:
        x, y (float, 1D array, or list): inputs to transform
        xoffset, yoffset (float): offsets in same units as x,y
        scale (float): scale factor
        rotation (float): counter-clockwise rotation in degrees

    returns xx, yy where
        xx = scale ( x cos(rot) - y sin(rot) ) + xoffset
        yy = scale ( x sin(rot) + y cos(rot) ) + xoffset
    """
    #- support lists, arrays, and scalar float -> float
    if not np.isscalar(x):
        x = np.asarray(x)
    if not np.isscalar(y):
        y = np.asarray(y)

    sinrot = np.sin(np.radians(rotation))
    cosrot = np.cos(np.radians(rotation))
    
    xx = scale*(x*cosrot - y*sinrot) + xoffset
    yy = scale*(x*sinrot + y*cosrot) + yoffset
    
    return xx, yy

def undo_scale_rotation_offset(x, y, scale=1.0, rotation=0.0, xoffset=0.0, yoffset=0.0):
    """
    Applies opposite transform of `apply_scale_rotation_offset``
    
    Args:
        x, y (float, 1D array, or list): inputs to transform
        xoffset, yoffset (float): offsets in same units as x,y
        scale (float): scale factor
        rotation (float): clockwise rotation in degrees

    returns reverse-transformed x, y
    """
    #- support lists, arrays, and scalar float -> float
    if not np.isscalar(x):
        x = np.asarray(x)
    if not np.isscalar(y):
        y = np.asarray(y)
    
    x1 = (x - xoffset)/scale
    y1 = (y - yoffset)/scale
    sinrot = np.sin(np.radians(-rotation))
    cosrot = np.cos(np.radians(-rotation))
    
    xx = (x1*cosrot - y1*sinrot)
    yy = (x1*sinrot + y1*cosrot)
    
    return xx, yy

def fit_scale_rotation_offset(x1, y1, x2, y2, p0=None):
    """
    TODO: document
    
    p0: starting guess (scale, rotation, xoffset, yoffset)
    """
    def fn(params, x1, y1, x2, y2):
        scale, rot, xoff, yoff = params
        xx, yy = apply_scale_rotation_offset(x1, y1, scale, rot, xoff, yoff)
        chi2 = np.sum((xx-x2)**2 + (yy-y2)**2)
        return chi2

    x1 = np.asarray(x1)
    x2 = np.asarray(x2)
    y1 = np.asarray(y1)
    y2 = np.asarray(y2)
    
    dx = np.mean(x1) - np.mean(x2)
    dy = np.mean(y1) - np.mean(y2)

    if p0 is None:
        p0 = (0.015, 0.0, dx, dy)
    results = minimize(fn, p0, args=(x1, y1, x2, y2))
    scale, rot, xoff, yoff = results.x
    return dict(scale=scale, rotation=rot, xoffset=xoff, yoffset=yoff)

def fit_gfa2fp(metrology):
    """
    TODO: document
    """
    #- HARDCODE: GFA pixel dimensions
    nx, ny = 2048, 1032

    #- Trim to just GFA entries without altering input table
    gfarows = (metrology['DEVICE_TYPE'] == 'GFA')
    metrology = metrology[gfarows]
    metrology.sort(['PETAL_LOC', 'PINHOLE_ID'])

    #- Metrology corners start at (0,0) for middle of pixel
    #- Thankfully this is consistent with gfa_reduce, desimeter fvc spots,
    #- and the Unified Metrology Table in DESI-5421
    xgfa = np.array([0, nx-1, nx-1, 0])
    ygfa = np.array([0, 0, ny-1, ny-1])
    
    gfa_transforms = dict()
    
    for p in range(10):
        ii = (metrology['PETAL_LOC'] == p)
        if np.count_nonzero(ii) > 0:
            xfp = np.asarray(metrology['X_FP'][ii])
            yfp = np.asarray(metrology['Y_FP'][ii])
            
            #- minimization starting parameters
            theta = np.radians(p*36 - 77)
            rmax = 407
            p0 = (0.015, 18+p*36.0, rmax*np.cos(theta), rmax*np.sin(theta))
        
            gfa_transforms[p] = fit_scale_rotation_offset(xgfa, ygfa, xfp, yfp, p0=p0)
    
    return gfa_transforms
        
        




