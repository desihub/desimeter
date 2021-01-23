import numpy as np
from skimage.feature import canny
from skimage.transform import rotate
import fitsio
import multiprocessing


def detect_phi_arm(x,y,image,template,ang_step=1.,plot=False) :
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

    curx = int(np.round(x))
    cury = int(np.round(y))
    stamp = image[cury-pad : cury+pad+1, curx-pad : curx+pad+1 ].copy()
    mean_image = np.mean(stamp)

    stamp -= np.median(stamp)
    stamp /= np.std(stamp)

    xgrid, ygrid = np.mgrid[-pad:pad+1,-pad:pad+1]
    mask = ((xgrid**2 + ygrid**2) < pad**2)

    template *= mask
    template /=np.sqrt(np.sum(template**2))

    best_angle = []
    best_ccf = []
    mean_ccf = []
    rms_ccf = []
    best_rtemp = []
    xx=[]
    yy=[]

    coarse_step = 5.
    coarse_angles = np.arange(0, 360, coarse_step)

    for canny_sigma in [0.1,0.2,0.5,1.,1.5] :
        edges = canny(stamp, sigma=canny_sigma).astype(float)
        edges *= mask
        norm = np.sqrt(np.sum(edges**2))
        edges /= norm

        coarse_step = 5.
        coarse_angles = np.arange(0, 360, coarse_step)
        templates = np.zeros((len(coarse_angles), ) + template.shape, dtype=template.dtype)
        for i, ang in enumerate(coarse_angles):
            templates[i] = rotate(template, -ang)
        coarse_ccf = (templates * edges[None, :, :]).sum(axis=1).sum(axis=1)
        pos = np.argmax(coarse_ccf)

        xx.append(coarse_angles)
        yy.append(coarse_ccf)

        # rerun with finer grid
        angle = coarse_angles[pos]
        angles = np.arange(angle-2*coarse_step, angle+2*coarse_step+ang_step, ang_step)
        templates = np.zeros((len(angles), ) + template.shape, dtype=template.dtype)
        for i, ang in enumerate(angles):
            templates[i] = rotate(template, -ang)
        ccf = (templates * edges[None, :, :]).sum(axis=1).sum(axis=1)
        pos = np.argmax(ccf)
        best_ccf.append(ccf[pos])
        mean_ccf.append(np.mean(ccf))
        rms_ccf.append(np.std(ccf))
        best_angle.append(angles[pos])
        best_rtemp.append(templates[pos])

    bb = np.argmax(best_ccf)
    best_edges = best_rtemp[bb]
    angle   = np.deg2rad(best_angle[bb]-90)
    ccfval  = best_ccf[bb]
    rmsccf  = rms_ccf[bb]
    meanccf = mean_ccf[bb]
    chi2 = 0.

    if plot :
        import matplotlib.pyplot as plt
        plt.figure()
        plt.subplot(231,title="x={} y={}".format(curx,cury))
        plt.imshow(stamp,origin=0,vmax=3*np.std(stamp))
        plt.subplot(232)
        plt.imshow(edges,origin=0)
        plt.subplot(233)
        plt.imshow(best_edges,origin=0)
        plt.subplot(212)
        for x,y in zip(xx,yy) :
            plt.plot(x,y)
        plt.axvline(best_angle[bb],linestyle="--")
        plt.xlabel("angle (deg)")
        plt.grid()
        print("angle= {:3.1f} deg ccf={:3.2f}".format(best_angle[bb],ccfval))
        plt.show()

    return angle, ccfval, meanccf, rmsccf, mean_image, norm, chi2

def detect_phi_arm_with_index(index,x,y,image,template,ang_step,plot=False) :
    angle, ccfval, meanccf, rmsccf, mean_image, norm, chi2 = detect_phi_arm(x,y,image,template,ang_step,plot=plot)
    return index,angle,ccfval, meanccf, rmsccf, mean_image, norm, chi2

def _func(arg) :
    """ Used for multiprocessing.Pool """
    index, angle, ccfval, meanccf, rmsccf, mean_image, norm, chi2= detect_phi_arm_with_index(**arg)
    print("{} angle={:4.1f} ccf={:4.3f}".format(index,angle,ccfval))
    return index, angle, ccfval, meanccf, rmsccf, mean_image, norm, chi2

def detect_phi_arms(spots,image_filename,template_filename,ang_step=1.,nproc=1,plot=False) :

    template = fitsio.read(template_filename).astype(float)

    image = None
    fits=fitsio.FITS(image_filename)
    for hdu in fits :
        if hdu.get_exttype() == 'IMAGE_HDU' and len(hdu.get_dims())==2 :
            image = hdu.read()
            break
    if image is None :
        print("error reading",image_filename)
        sys.exit(12)

    image=image.astype(float)


    det    = spots
    det["ANGLE"] = np.nan
    det["CCF"] = np.nan
    det['MCCF']=np.nan
    det['RMSCCF']=np.nan
    det['MIMAGE']=np.nan
    det['NORM']=np.nan
    det['CHI2']=np.nan

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
            det['MCCF'][i]=result[3]
            det['RMSCCF'][i]=result[4]
            det['MIMAGE'][i]=result[5]
            det['NORM'][i]=result[6]
            det['CHI2'][i]=result[7]


    else :
        for i in range(ndet) :
            det['ANGLE'][i], det['CCF'][i], det['MCCF'][i], det['RMSCCF'][i], det['MIMAGE'][i], det['NORM'][i], det['CHI2'][i] = detect_phi_arm(det["XPIX"][i],det["YPIX"][i],image,template,ang_step=ang_step,plot=plot)
            print("{}/{} angle={:4.1f} ccf={:4.3f}".format(i,ndet,det['ANGLE'][i], det['CCF'][i]))
