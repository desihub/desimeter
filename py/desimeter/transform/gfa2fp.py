"""
Transforms between GFA pixel coordinates and FP mm
"""

import numpy as np
from desimeter import io
from desimeter.log import get_logger
from desimeter.simplecorr import SimpleCorr

#- Cached GFA pix -> GFA FP metrology scale, rotation, offsets per petal
_gfa_transforms = None
def _get_transform(petal_loc):
    global _gfa_transforms
    if _gfa_transforms is None:
        metrology = io.load_metrology()
        _gfa_transforms = fit_gfa2fp(metrology)

    log = get_logger()
    if petal_loc not in _gfa_transforms:
        log.error('PETAL_LOC {} GFA metrology missing'.format(petal_loc))

    return _gfa_transforms[petal_loc]

def gfa2fp(petal_loc, xgfa, ygfa):
    """
    Transforms from GFA pixel coordinates to focal plane mm

    Args:
        petal_loc (int): Petal location 0-9
        xgfa, ygfa: GFA pixel coordinates, (0,0) is corner pixel center

    Returns CS5 xfp, yfp in mm
    """
    trans = _get_transform(petal_loc)
    xfp, yfp = trans.apply(xgfa, ygfa)
    return xfp, yfp

def gfa2fp_has_metrology(petal_loc):
    trans = _get_transform(petal_loc)
    return trans.from_metrology

def fp2gfa(petal_loc, xfp, yfp):
    """
    Transforms from focal plane mm to GFA pixel coordinates

    Args:
        petal_loc (int): Petal location 0-9
        xfp, yfp: CS5 focal plane mm

    Returns xgfa, ygfa pixel coordinates with (0,0) as center of corner pixel
    """
    global _gfa_transforms
    if _gfa_transforms is None:
        metrology = io.load_metrology()
        _gfa_transforms = fit_gfa2fp(metrology)

    log = get_logger()
    if petal_loc not in _gfa_transforms:
        log.error('PETAL_LOC {} GFA metrology missing'.format(petal_loc))

    xgfa, ygfa = _gfa_transforms[petal_loc].apply_inverse(xfp, yfp)

    return xgfa, ygfa

def fit_gfa2fp(metrology):
    """
    Fit GFA pix -> FP mm scale, rotation, xyoffsets for each GFA

    Returns dict keyed by PETAL_LOC, with dictionaries of transform coeffs.
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

    # Fill in missing petal metrology with averages of measured ones.
    xx,yy,zz = [],[],[]
    for p in range(10):
        ii = (metrology['PETAL_LOC'] == p)
        if not np.any(ii):
            continue
        xx.append(np.asarray(metrology['X_PTL'][ii]))
        yy.append(np.asarray(metrology['Y_PTL'][ii]))
        zz.append(np.asarray(metrology['Z_PTL'][ii]))
    xx = np.mean(np.vstack(xx), axis=0)
    yy = np.mean(np.vstack(yy), axis=0)
    zz = np.mean(np.vstack(zz), axis=0)

    for p in range(10):
        ii = (metrology['PETAL_LOC'] == p)
        if np.any(ii):
            xfp = np.asarray(metrology['X_FP'][ii])
            yfp = np.asarray(metrology['Y_FP'][ii])
            zfp = np.asarray(metrology['Z_FP'][ii])
            from_metrology = True
        else:
            # Approximate from average of measured ones!
            # XYZ_PTL are for petal_loc = 3
            rot = np.deg2rad(36. * (p - 3))
            R = np.array([[np.cos(rot), -np.sin(rot)],
                          [np.sin(rot), np.cos(rot)]])
            xyr = np.dot(R, np.vstack((xx,yy)))
            xfp = xyr[0,:]
            yfp = xyr[1,:]
            zfp = zz
            from_metrology = False

        #- fit transform
        corr = SimpleCorr()
        corr.fit(xgfa, ygfa, xfp, yfp)

        #- measure norm of plane
        x01 =  np.array( [ xfp[1]-xfp[0], yfp[1]-yfp[0], zfp[1]-zfp[0] ] )
        x01 /= np.sqrt(np.sum(x01**2))
        x12 =  np.array( [ xfp[2]-xfp[1], yfp[2]-yfp[1], zfp[2]-zfp[1] ] )
        x12 /= np.sqrt(np.sum(x12**2))
        norm_vector= np.cross(x01,x12)
        # I checked the sign of all components

        #- compute correction to apply
        delta_z = 2.23 # mm
        delta_x = delta_z*norm_vector[0]/norm_vector[2]
        delta_y = delta_z*norm_vector[1]/norm_vector[2]

        #- apply correction to offsets
        corr.dx += delta_x
        corr.dy += delta_y
        corr.from_metrology = from_metrology
        gfa_transforms[p] = corr

    return gfa_transforms
