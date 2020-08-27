import numpy as np
from skimage.feature import canny
from skimage.transform import rotate
import fitsio
import multiprocessing


def detect_phi_arm(x,y,image,template,ang_step=1.) :
    """
    detection of phi arm angle

    Args:
     x: float, x pixel coordinate of fiber tip
     y: float, y pixel coordinate of fiber tip
     image: 2D np.array : bright image
     template: 2D np.array : template
     ang_step: float, step in degree of azimuthal angle scan

    returns angle, ccfval with angle in radians and ccfval the cross
    correlation coefficient
    """

    assert(template.shape[0] == template.shape[1])
    pad=template.shape[0]//2


    D = image

    curx = int(np.round(x))
    cury = int(np.round(y))
    curD = D[cury-pad : cury+pad+1, curx-pad : curx+pad+1 ]
    curD -= np.median(curD)
    curD /= np.median(np.abs(curD))
    can  =  canny(curD, sigma=2)

    # apply a circular mask
    xgrid, ygrid = np.mgrid[-pad:pad+1,-pad:pad+1]
    can *= ((xgrid**2 + ygrid**2) < pad**2)

    angles = np.arange(0, 360, ang_step)
    imas = np.zeros((len(angles), ) + can.shape, dtype=can.dtype)
    for i, ang in enumerate(angles):
        imas[i] = rotate(can, ang)
    imas = imas / (np.sum(imas**2, axis=1).sum(axis=1)**.5)[:, None, None]
    ccf = (template[None, :, :] * imas).sum(axis=1).sum(axis=1)
    pos = np.argmax(ccf)
    angle  = np.deg2rad(angles[pos]-90)
    ccfval = ccf[pos]

    return angle, ccfval

def detect_phi_arm_with_index(index,x,y,image,template,ang_step) :
    angle, ccfval = detect_phi_arm(x,y,image,template,ang_step)
    return index,angle,ccfval

def _func(arg) :
    """ Used for multiprocessing.Pool """
    index, angle, ccfval = detect_phi_arm_with_index(**arg)
    print("{} angle={:4.1f} ccf={:4.3f}".format(index,angle,ccfval))
    return index, angle, ccfval

def detect_phi_arms(spots,image_filename,template_filename,ang_step,nproc=1) :

    template          = fitsio.read(template_filename).astype(float)
    template /= np.sqrt(np.sum(template**2))

    image  = fitsio.read(image_filename).astype(float)

    det    = spots
    det["ANGLE"] = np.nan
    det["CCF"] = np.nan

    ndet=len(det)

    if nproc > 1 :
        pool = multiprocessing.Pool(nproc)
        func_args = []
        for i in range(ndet) :
            arguments={"index":i,"x":det["XPIX"][i],"y":det["YPIX"][i],"image":image,"template":template,"ang_step":ang_step}
            func_args.append( arguments )
        results  =  pool.map(_func, func_args)
        pool.close()
        pool.join()
        for result in results :
            i=result[0]
            det['ANGLE'][i]=result[1]
            det['CCF'][i]=result[2]
    else :
        for i in range(ndet) :
            det['ANGLE'][i], det['CCF'][i] = detect_phi_arm(det["XPIX"][i],det["YPIX"][i],image,template,ang_step=ang_step)
            print("{}/{} angle={:4.1f} ccf={:4.3f}".format(i,ndet,det['ANGLE'][i], det['CCF'][i]))
