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
        filename = resource_filename('desicoord',"data/fiducials-fpqs.csv")
        if not os.path.isfile(filename) :
            log.error("cannot find {}".format(filename))
            raise IOError("cannot find {}".format(filename))
        log.info("reading fiducials metrology in {}".format(filename))
        fiducials_fpqs_table = Table.read(filename,format="csv")
    

    # select matched spots and central pinhole identified with PIN_ID=1027 for now (will change)
    selection = np.where((spots["FID_ID"]>0)&(spots["PIN_ID"]==1027))[0]
    if len(selection)<4 : 
        log.error("Only {} fiducials were matched, I cannot fit a transfo".format(len(selection)))
        raise RuntimError("Only {} fiducials were matched, I cannot fit a transfo".format(len(selection)))
    
    
    xpix = spots["XPIX"][selection]
    ypix = spots["YPIX"][selection]
    xerr = spots["XERR"][selection]
    yerr = spots["YERR"][selection]

    # now match Q and S
    fid_id = { f : i for i,f in enumerate(fiducials_fpqs_table["FID_ID"]) } # dictionary of indices
    indices = [ fid_id[f] for f in spots["FID_ID"][selection] ]
    q = fiducials_fpqs_table["Q"][indices]
    s = fiducials_fpqs_table["S"][indices]

    # convert to X and Y
    xfp,yfp = qs2xy(q,s)

    log.warning("Trivial transformation fit with polynomials (pretty bad residuals)")

    rx,ry,ex,ey = _reduced_coords(xpix,ypix,xerr,yerr)
    
    """
    # first trivial steps
    xtmp = xpix - np.mean(xpix)
    ytmp = ypix - np.mean(ypix)
    xtmp *= -1
    a = np.std(xfp)/np.std(xtmp)
    xtmp *= a
    xerr *= a
    a = np.std(yfp)/np.std(ytmp)
    ytmp *= a
    yerr *= a

    xerr[xerr<=0] = 100. # deweight
    yerr[yerr<=0] = 100. # deweight
    
    xtmp += np.median(xfp-xtmp)
    ytmp += np.median(yfp-ytmp)
    """
    
    # solve linear system
    nspot=xpix.size
    deg=3 # that includes translation,rotation,dilatation

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
    
    if 'XFP' not in spots.dtype.names : spots.add_column(Column(np.zeros(nspot,dtype=float)),name='XFP')
    if 'YFP' not in spots.dtype.names : spots.add_column(Column(np.zeros(nspot,dtype=float)),name='YFP')
    spots["XFP"] = xfp_meas
    spots["YFP"] = yfp_meas

    if True : # add matched metrology data
        if 'XMETRO' not in spots.dtype.names : spots.add_column(Column(np.zeros(nspot,dtype=float)),name='XMETRO')
        if 'YMETRO' not in spots.dtype.names : spots.add_column(Column(np.zeros(nspot,dtype=float)),name='YMETRO')
        if 'QMETRO' not in spots.dtype.names : spots.add_column(Column(np.zeros(nspot,dtype=float)),name='QMETRO')
        if 'SMETRO' not in spots.dtype.names : spots.add_column(Column(np.zeros(nspot,dtype=float)),name='SMETRO')
        spots["XMETRO"][selection] = xfp
        spots["YMETRO"][selection] = yfp
        spots["QMETRO"][selection] = q
        spots["SMETRO"][selection] = s
    
    if False : # cannot run without error: "A value in x_new is above the interpolation "
        if 'QFP' not in spots.dtype.names : spots.add_column(Column(np.zeros(nspot,dtype=float)),name='QFP')
        if 'SFP' not in spots.dtype.names : spots.add_column(Column(np.zeros(nspot,dtype=float)),name='SFP')
        qfp_meas , sfp_meas = xy2qs(xfp_meas,yfp_meas)
        spots["QFP"] = qfp_meas
        spots["SFP"] = sfp_meas
    
    #import matplotlib.pyplot as plt
    #plt.plot(xmeas_fp,ymeas_fp,"o")
    #plt.plot(xfp,yfp,"x")
    #plt.figure()
    #plt.hist(dist,bins=40)
    #plt.show()

    return spots


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

    nspot = spots['Petal Loc ID'].size

    # local petal coordinates 'pl'
    xyzpl = np.zeros((3,nspot))
    xyzpl[0] = spots['X FCL']
    xyzpl[1] = spots['Y FCL']
    xyzpl[2] = spots['Z FCL']

    # global focal plane coordinates 'fp'
    xyzfp = np.zeros((3,nspot))
    
    for petal in np.unique(spots['Petal Loc ID']) :
        ii = np.where(spots['Petal Loc ID']==petal)[0]
        params = petal_alignment_dict[petal]
        Rotation = Rxyz(params["alpha"],params["beta"],params["gamma"])
        Translation = np.array([params["Tx"],params["Ty"],params["Tz"]])
        xyzfp[:,ii] = Rotation.dot(xyzpl) + Translation
    
    import matplotlib.pyplot as plt
    plt.plot(xyzfp[0],xyzfp[1],"o")
    plt.show()
    
