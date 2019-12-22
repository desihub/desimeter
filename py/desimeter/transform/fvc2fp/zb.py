"""
Utility functions to fit and apply coordinates transformation from FVC to FP
"""

import os
import json
from pkg_resources import resource_filename

import numpy as np
from astropy.table import Table,Column
from scipy.optimize import minimize

# from desimodel.focalplane.geometry import qs2xy,xy2qs
from desiutil.log import get_logger

from .base import FVC2FP_Base
from desimeter.transform.zhaoburge import getZhaoBurgeXY, getZhaoBurgeTerm

#- Utility transforms to/from reduced [-1,1] coordinates
def _reduce_xyfp(x, y):
    """
    Rescale FP xy coordinates [-420,420] -> [-1,1] and flip x axis
    """
    a = 420.0
    return -x/a, y/a

def _expand_xyfp(x, y):
    """
    Undo _redux_xyfp() transform
    """
    a = 420.0
    return -x*a, y*a

def _reduce_xyfvc(x, y):
    """
    Rescale FVC xy pix coords [0,6000] -> [-1,1]
    """
    c = 3000.0
    return (x-c)/c, (y-c)/c

def _expand_xyfvc(x, y):
    """
    Undo _redux_xyfvc() transform
    """
    c = 3000.0
    return (x+1)*c, (y+1)*c

#- Tranformation in reduced coordinates space, to be used by minimizer
def transform(x, y, scale, rotation, offset_x, offset_y, zbcoeffs=None):
    """
    TODO: document
    """
    xx = (x*np.cos(rotation) - y*np.sin(rotation))*scale + offset_x
    yy = (x*np.sin(rotation) + y*np.cos(rotation))*scale + offset_y

    if zbcoeffs is not None:
        dx, dy = getZhaoBurgeXY(zbcoeffs, xx, yy)
        xx += dx
        yy += dy

    return xx, yy

def fit_scale_rotation_offset(x, y, xp, yp, fitzb=False):

    def func(params, x, y, xp, yp, fitzb):
        scale, rotation, offset_x, offset_y = params[0:4]
        xx, yy = transform(x, y, scale, rotation, offset_x, offset_y)
        
        if fitzb:
            c, zbx, zby = fitZhaoBurge(xx, yy, xp, yp)
            xx += zbx
            yy += zby
        
        dr2 = np.sum((xx-xp)**2 + (yy-yp)**2)
        return dr2

    p0 = np.array([1.0, 0.0, 0.0, 0.0])
    p = minimize(func, p0, args=(x, y, xp, yp, fitzb)) #, method='Nelder-Mead')

    scale, rotation, offset_x, offset_y = p.x
    if fitzb:
        xx, yy = transform(x, y, scale, rotation, offset_x, offset_y)
        zbcoeffs, zbx, zby = fitZhaoBurge(xx, yy, xp, yp)
        return scale, rotation, offset_x, offset_y, zbcoeffs
    else:
        return scale, rotation, offset_x, offset_y

def fitZhaoBurge(x, y, xp, yp):
    dx = xp-x
    dy = yp-y

    nx = len(x)
    nzb = 8
    H = np.zeros((2*nx, nzb))
    for i in range(nzb):
        zbx, zby, name = getZhaoBurgeTerm(i, x, y)
        H[0:nx, i] = zbx
        H[nx:, i] = zby

    A = H.T.dot(H)
    b = H.T.dot(np.concatenate([dx, dy]))
    c = np.linalg.solve(A, b)

    zbx, zby = getZhaoBurgeXY(c, x, y)
    
    return c, zbx, zby

#-------------------------------------------------------------------------

class FVCFP_ZhaoBurge(FVC2FP_Base):
    def tojson(self):
        params = dict()
        params['method'] = 'Zhao-Burge'
        params['scale'] = self.scale
        params['rotation'] = self.rotation
        params['offset_x'] = self.offset_x
        params['offset_y'] = self.offset_y
        params['zbcoeffs'] = list(self.zbcoeffs)
        return json.dumps(params)
    
    @classmethod
    def fromjson(cls, jsonstring):
        tx = cls()
        params = json.loads(jsonstring)
        assert params['method'] == 'Zhao-Burge'
        tx.scale = params['scale']
        tx.rotation = params['rotation']
        tx.offset_x = params['offset_x']
        tx.offset_y = params['offset_y']
        tx.zbcoeffs = np.asarray(params['zbcoeffs'])
        return tx
    
    def fit(self, spots, metrology=None, update_spots=False):
        """TODO: document"""
        log = get_logger()
        if metrology is not None:
            self.metrology = metrology
        else:
            filename = resource_filename('desimeter',"data/fp-metrology.csv")
            if not os.path.isfile(filename) :
                log.error("cannot find {}".format(filename))
                raise IOError("cannot find {}".format(filename))
            log.info("reading fiducials metrology in {}".format(filename))
            self.metrology = Table.read(filename,format="csv")

        #- Trim spots to just fiducial spots (not posioners, not unmatchs spots)
        ii = (spots['LOCATION']>0) & (spots['PINHOLE_ID']>0)
        fidspots = spots[ii]

        #- trim metrology to just the ones that have spots
        fidspots_pinloc = fidspots['LOCATION']*10 + fidspots['PINHOLE_ID']
        metro_pinloc = self.metrology['LOCATION']*10 + self.metrology['PINHOLE_ID']
        jj = np.in1d(metro_pinloc, fidspots_pinloc)
        metrology = self.metrology[jj]
        
        #- Sort so that they match each other
        fidspots.sort(keys=('LOCATION', 'PINHOLE_ID'))
        metrology.sort(keys=('LOCATION', 'PINHOLE_ID'))
        assert np.all(fidspots['LOCATION'] == metrology['LOCATION'])
        assert np.all(fidspots['PINHOLE_ID'] == metrology['PINHOLE_ID'])
        
        #- Get reduced coordinates
        rxpix, rypix = _reduce_xyfvc(fidspots['XPIX'], fidspots['YPIX'])
        rxfp, ryfp = _reduce_xyfvc(metrology['X_FP'], metrology['Y_FP'])
        
        #- Perform fit        
        scale, rotation, offset_x, offset_y, zbcoeffs = \
            fit_scale_rotation_offset(rxpix, rypix, rxfp, ryfp, fitzb=True)
                
        self.scale = scale
        self.rotation = rotation
        self.offset_x = offset_x
        self.offset_y = offset_y      
        self.zbcoeffs = zbcoeffs      

        #- Goodness of fit
        xfp_fidmeas, yfp_fidmeas = self.fvc2fp(fidspots['XPIX'], fidspots['YPIX'])
        dx = (metrology['X_FP'] - xfp_fidmeas)
        dy = (metrology['Y_FP'] - yfp_fidmeas)
        dr = np.sqrt(dx**2 + dy**2)
        log.info('Mean and median distance = {:.1f}, {:.1f} um'.format(
            np.mean(dr), np.median(dr)))

        if update_spots:
            xfp_meas, yfp_meas = self.fvc2fp(spots['XPIX'], spots['YPIX'])
            spots["X_FP"] = xfp_meas
            spots["Y_FP"] = yfp_meas
            # spots["X_FP_METRO"] = np.zeros(xfp_meas.size)
            # spots["Y_FP_METRO"] = np.zeros(xfp_meas.size)
            # spots["X_FP_METRO"][selection] = xfp
            # spots["Y_FP_METRO"][selection] = yfp

    def fvc2fp(self, xpix, ypix, xerr=None, yerr=None):
        """
        Converts fiber view camera pixel x,y -> focal plane x,y
        """
        rx, ry = _reduce_xyfvc(xpix, ypix)
        rxfp, ryfp = transform(rx, ry, self.scale, self.rotation,
            self.offset_x, self.offset_y, self.zbcoeffs)
        xfp, yfp = _expand_xyfp(rxfp, ryfp)
        return xfp, yfp

    def fp2fvc(self, xfp, yfp):
        """
        Converts focal plane x,y -> fiber view camera pixel x,y
        """
        rxfp, ryfp = _reduce_xyfp(xfp, yfp)
        
        #- first undo Zhao-Burge terms
        dx, dy = getZhaoBurgeXY(self.zbcoeffs, rxfp, ryfp)
        rxfp -= dx
        ryfp -= dy

        #- Then apply inverse scale, roation, offset
        rxfp -= self.offset_x
        ryfp -= self.offset_y
        
        rxfp /= self.scale
        ryfp /= self.scale

        xx = (rxfp*np.cos(-rotation) - ryfp*np.sin(-self.rotation))
        yy = (rxfp*np.sin(-rotation) + ryfp*np.cos(-self.rotation))
        
        xpix, ypix = _expand_xyfvc(xx, yy)

        return xpix, ypix