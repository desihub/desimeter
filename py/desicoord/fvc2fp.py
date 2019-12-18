"""
Utility functions to fit and apply coordinates transformation from FVC to FP
"""

import os
import numpy as np
from desimodel.focalplane.geometry import qs2xy,xy2qs
from desiutil.log import get_logger
from astropy.table import Table,Column
from pkg_resources import resource_filename


fiducials_fpqs_table = None


def _reduced_coords(xpix,ypix,xerr,yerr) :
    a= 1./3000.
    return -(a*xpix-1) , a*ypix-1, a*xerr , a*yerr
    
    
def fit_fvc2fp(spots) :

    # spots are supposed to be matched

    global fiducials_fpqs_table


    log = get_logger()

    
    if fiducials_fpqs_table is None :
        filename = resource_filename('desicoord',"data/fp-metrology.csv")
        if not os.path.isfile(filename) :
            log.error("cannot find {}".format(filename))
            raise IOError("cannot find {}".format(filename))
        log.info("reading fiducials metrology in {}".format(filename))
        fiducials_fpqs_table = Table.read(filename,format="csv")
    
    selection = np.where((spots["LOCATION"]>0)&(spots["DOTID"]==4))[0]
    if len(selection)<4 : 
        log.error("Only {} fiducials were matched, I cannot fit a transfo".format(len(selection)))
        raise RuntimError("Only {} fiducials were matched, I cannot fit a transfo".format(len(selection)))
    
    
    xpix = spots["XPIX"][selection]
    ypix = spots["YPIX"][selection]
    xerr = spots["XERR"][selection]
    yerr = spots["YERR"][selection]

    # now match Q and S
    fid_id = { f : i for i,f in enumerate(fiducials_fpqs_table["LOCATION"]) } # dictionary of indices
    indices = [ fid_id[f] for f in spots["LOCATION"][selection] ]
    xfp = fiducials_fpqs_table["XFP"][indices]
    yfp = fiducials_fpqs_table["YFP"][indices]

    log.warning("Trivial transformation fit with polynomials (pretty bad residuals)")

    rx,ry,ex,ey = _reduced_coords(xpix,ypix,xerr,yerr)
    
    # solve linear system
    nspot=xpix.size
    deg=4 # that includes translation,rotation,dilatation

    npar=0
    for dx in range(0,deg+1) :
        for dy in range(0,deg+1-dx) :
            npar += 1
    log.info("Number of parameters for the transformation = {}".format(npar))
    H=np.zeros((npar,nspot))
    p=0
    for dx in range(0,deg+1) :
        for dy in range(0,deg+1-dx) :
            H[p] = (rx**dx) * (ry**dy)
            #print("{} x^{} * y^{}".format(p,dx,dy))
            p+=1
        
    # solve for x
    A  = H.dot(H.T)
    B  = ((xfp)*H).sum(axis=1)
    Ai = np.linalg.inv(A)
    Px = Ai.dot(B)
    log.debug("Px={}".format(Px))
    
    # solve for y
    A  = H.dot(H.T)
    B  = ((yfp)*H).sum(axis=1)
    Ai = np.linalg.inv(A)
    Py = Ai.dot(B)
    log.debug("Py={}".format(Py))
    
    # apply to all spots
    xpix = spots["XPIX"]
    ypix = spots["YPIX"]
    xerr = spots["XERR"]
    yerr = spots["YERR"]
    rx,ry,ex,ey = _reduced_coords(xpix,ypix,xerr,yerr)

    nspot=xpix.size
    H=np.zeros((npar,nspot))
    p=0
    for dx in range(0,deg+1) :
        for dy in range(0,deg+1-dx) :
            H[p] = (rx**dx) * (ry**dy)
            p+=1
    xfp_meas = Px.dot(H)
    yfp_meas = Py.dot(H)
    
    dist = np.sqrt((xfp_meas[selection]-xfp)**2 + (yfp_meas[selection]-yfp)**2)
    mdist = np.mean(dist)
    log.info("Mean distance = {} mm".format(mdist))
    
    spots["XFP"] = xfp_meas
    spots["YFP"] = yfp_meas

    if True : # add matched metrology data
        spots["XMETRO"] = np.zeros(xfp_meas.size)
        spots["YMETRO"] = np.zeros(xfp_meas.size)
        spots["XMETRO"][selection] = xfp
        spots["YMETRO"][selection] = yfp
    
    return spots


