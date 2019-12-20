"""
Utility functions to fit and apply coordinates transformation from FVC to FP
"""

import os
import json
from pkg_resources import resource_filename

import numpy as np
from astropy.table import Table,Column

# from desimodel.focalplane.geometry import qs2xy,xy2qs
from desiutil.log import get_logger

def _vander2d(rx, ry, degree=3):
    """
    Return 2D equivalent of a Vandermonde matrix (loosely speaking)
    
    TODO: document
    """
    assert len(rx) == len(ry)
    npar = ((degree+1) * (degree+2))//2
    H = np.zeros((npar,len(rx)))
    p=0
    for dx in range(0,degree+1):
        for dy in range(0,degree+1-dx):
            H[p] = (rx**dx) * (ry**dy)
            p+=1
    
    return H

def _polyfit2d(rx, ry, px, py, degree=3):
    """TODO: document"""

    assert len(rx) == len(ry)
    
    #- de-astropyify
    rx = np.asarray(rx)
    ry = np.asarray(ry)
    px = np.asarray(px)
    py = np.asarray(py)

    # Setup matrix to solve
    # TODO: probably better to use np.linalg.solve instead of .inv
    H = _vander2d(rx, ry, degree)
    A  = H.dot(H.T)
    Ai = np.linalg.inv(A)

    # solve for x
    Bx  = ((px)*H).sum(axis=1)
    Px = Ai.dot(Bx)

    # solve for y
    By  = (py*H).sum(axis=1)
    Py = Ai.dot(By)
    
    return Px, Py

class FVC2FP_Base(object):
    """
    Base class for transforms between FVC and FP coordinates

    Subclasses implement this interface with specific parameterizations
    """
    def __init__(self):
        pass
    
    def tojson(self):
        raise NotImplementedError
    
    @classmethod
    def fromjson(cls, jsonstring):
        raise NotImplementedError
    
    @classmethod
    def read_jsonfile(cls, filename):
        with open(filename) as fx:
            s = fx.read()
        return cls.fromjson(s)

    def write_jsonfile(self, filename):
        with open(filename, 'w') as fx:
            fx.write(self.tojson())

    def fit(spots, metrology=None):
        raise NotImplementedError

    def fvc2fp(xpix, ypix, xerr=None, yerr=None):
        """
        Converts fiber view camera pixel x,y -> focal plane x,y
        """
        raise NotImplmentedError

    def fpxy2fvc(xfp, yfp):
        """
        Converts focal plane x,y -> fiber view camera pixel x,y
        """
        raise NotImplmentedError

class FVCFP_Polynomial(FVC2FP_Base):
    """
    Naive 2D polynomial transform between FVC and FP, mainly to demonstrate
    interfaces rather than quality of actual transform.
    """

    def tojson(self):
        params = dict()
        params['method'] = 'xy polynomial'
        params['degree'] = self.degree
        params['fvc2fp_x_coeff'] = list(self.Px)
        params['fvc2fp_y_coeff'] = list(self.Py)
        params['fp2fvc_x_coeff'] = list(self.Qx)
        params['fp2fvc_y_coeff'] = list(self.Qy)
        return json.dumps(params)
    
    @classmethod
    def fromjson(cls, jsonstring):
        tx = cls()
        params = json.loads(jsonstring)
        assert params['method'] == 'xy polynomial'
        tx.degree = params['degree']
        tx.Px = np.asarray(params['fvc2fp_x_coeff'])
        tx.Py = np.asarray(params['fvc2fp_y_coeff'])
        tx.Qx = np.asarray(params['fp2fvc_x_coeff'])
        tx.Qy = np.asarray(params['fp2fvc_y_coeff'])
        return tx
    
    def _fvc_reduced_coords(self,xpix,ypix,xerr=None,yerr=None) :
        """
        convert FVC pixel range [0,6000] -> [-1,1]
        """
        a = 1./3000.
        
        xx, yy = -(a*xpix-1) , a*ypix-1
        if xerr is not None:
            assert yerr is not None
            xxerr, yyerr = a*xerr, a*yerr
            return xx, yy, xxerr, yyerr
        else:
            return xx, yy

    def _fp_reduced_coords(self,xpix,ypix,xerr=None,yerr=None) :
        """
        convert FP mm range [-500,500] -> [-1,1]
        """
        a = 1./500
        xx, yy = a*xpix , a*ypix
        if xerr is not None:
            assert yerr is not None
            xxerr, yyerr = a*xerr, a*yerr
            return xx, yy, xxerr, yyerr
        else:
            return xx, yy

    def fit(self, spots, degree=3, metrology=None, update_spots=False):
    
        log = get_logger()
        self.degree = degree
        if metrology is not None:
            self.metrology_table = metrology
        else:
            filename = resource_filename('desimeter',"data/fp-metrology.csv")
            if not os.path.isfile(filename) :
                log.error("cannot find {}".format(filename))
                raise IOError("cannot find {}".format(filename))
            log.info("reading fiducials metrology in {}".format(filename))
            self.metrology_table = Table.read(filename,format="csv")
    
        selection = np.where((spots["LOCATION"]>0))[0]
        if len(selection)<4 : 
            log.error("Only {} fiducials were matched, I cannot fit a transform".format(len(selection)))
            raise RuntimError("Only {} fiducials were matched, I cannot fit a transform".format(len(selection)))
    
        xpix = spots["XPIX"][selection]
        ypix = spots["YPIX"][selection]
        xerr = spots["XERR"][selection]
        yerr = spots["YERR"][selection]

        # now match

        spots_identifier = np.array(spots["LOCATION"][selection]*100 + spots["DOTID"][selection])
        metro_identifier = np.array(self.metrology_table["LOCATION"]*100 + self.metrology_table["DOTID"])
    
        tmpid = { f : i for i,f in enumerate(metro_identifier) } # dictionary of indices
    
        spots_indices = []
        metro_indices = []

        for i,f in enumerate(spots_identifier) :
            if f in tmpid :
               spots_indices.append(i)
               metro_indices.append(tmpid[f])
            else :
                log.warning("cannot find metrology for LOCATION={} DOTID={}".format(int(f//100),int(f%100)))
    
        xfp = self.metrology_table["XFP"][metro_indices]
        yfp = self.metrology_table["YFP"][metro_indices]
        xpix = xpix[spots_indices]
        ypix = ypix[spots_indices]
        xerr = xerr[spots_indices]
        yerr = yerr[spots_indices]
        selection = selection[spots_indices]  
        log.warning("Trivial transformation fit with polynomials (pretty bad residuals)")

    
        #---
        # solve linear system for FVC [pix] -> FP [mm]
        
        rx,ry = self._fvc_reduced_coords(xpix,ypix)
        self.Px, self.Py = _polyfit2d(rx, ry, xfp, yfp, self.degree)
        
        log.info('FVC -> FP: {} parameters'.format(len(self.Px)))
        log.debug("Px={}".format(self.Px))
        log.debug("Py={}".format(self.Py))
    
        #---
        # Solve inverse linear system for FP [mm] -> FVC [pix]
        # Use higher degree and denser sampling of transform

        rpix = np.linspace(0, 6000, num=25)
        xpix, ypix = np.meshgrid(rpix, rpix)
        xpix = xpix.ravel()
        ypix = ypix.ravel()
        xfp, yfp = self.fvc2fp(xpix, ypix)
        
        rx, ry = self._fp_reduced_coords(xfp, yfp)
        log.debug('FP minmax(rx) = {:.3f}, {:.3f}'.format(np.min(rx), np.max(rx)))
        log.debug('FP minmax(ry) = {:.3f}, {:.3f}'.format(np.min(ry), np.max(ry)))
        self.Qx, self.Qy = _polyfit2d(rx, ry, xpix, ypix, self.degree+5)
        log.info('FP -> FVC: {} parameters'.format(len(self.Qx)))

        #---
        # check goodness of fit metrics
        
        xfp = self.metrology_table["XFP"][metro_indices]
        yfp = self.metrology_table["YFP"][metro_indices]        
        xfp_meas, yfp_meas = self.fvc2fp(spots["XPIX"], spots["YPIX"])
        
        dist = np.sqrt((xfp_meas[selection]-xfp)**2 + (yfp_meas[selection]-yfp)**2)
        mdist = np.mean(dist)
        log.info("Mean and median distance = {:.1f}, {:.1f} um".format(
            1000*np.mean(dist),1000*np.median(dist)))

        if update_spots:
            spots["XFP"] = xfp_meas
            spots["YFP"] = yfp_meas
            spots["XMETRO"] = np.zeros(xfp_meas.size)
            spots["YMETRO"] = np.zeros(xfp_meas.size)
            spots["XMETRO"][selection] = xfp
            spots["YMETRO"][selection] = yfp


    def fvc2fp(self, xpix, ypix):
    
        xpix = np.asarray(xpix)
        ypix = np.asarray(ypix)
    
        rx,ry = self._fvc_reduced_coords(xpix,ypix)
        H = _vander2d(rx, ry, self.degree)
        xfp = self.Px.dot(H)
        yfp = self.Py.dot(H)

        return xfp, yfp

    def fp2fvc(self, xfp, yfp):
    
        xfp = np.asarray(xfp)
        yfp = np.asarray(yfp)
    
        rx,ry = self._fp_reduced_coords(xfp,yfp)
        H = _vander2d(rx, ry, self.degree+5)
        xfvc = self.Qx.dot(H)
        yfvc = self.Qy.dot(H)

        return xfvc, yfvc


        
        
        

