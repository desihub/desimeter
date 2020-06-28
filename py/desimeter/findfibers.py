"""
Code to find fibers and rotations of fiber arms
"""

import astropy.io.fits as pyfits
import numpy as np
import sep
import matplotlib.pyplot as plt
import astropy.table as atpy
from skimage.feature import canny
from skimage.transform import rotate
import scipy.stats
import scipy.signal

def norm(x):
    """ Normalize the image by median and MAD"""
    med = np.median(x)
    std = np.median(np.abs(x - med))
    return (x - med) / std


def find_angle(ima, ref, ang_step=1):
    """ Find the best orientation by matching the rotated image to
    a template
    Parameters:
    ima: ndarray Input image
    ref: ndarray Reference image
    ang_step: float in degrees, step of the search; The smaller the step the slower the code.
    TODO: I should redo this by computing once the rotated reference image
    for all angles instead of rotating actual images all the time
    """
    angles = np.arange(0, 360, ang_step)
    ima1 = ima * 1.
    imas = np.zeros((len(angles), ) + ima1.shape, dtype=ima1.dtype)
    for i, ang in enumerate(angles):
        imas[i] = rotate(ima1, ang)
    imas = imas / (np.sum(imas**2, axis=1).sum(axis=1)**.5)[:, None, None]
    ref = ref / np.sum(ref**2)**.5
    ccf = (ref[None, :, :] * imas).sum(axis=1).sum(axis=1)
    pos = np.argmax(ccf)
    val = ccf[pos]
    return angles[pos], val


def perc(X, per):
    """ Percentile (for plotting)"""
    return scipy.stats.scoreatpercentile(X, per)


def find_fibers(fname_front,
                fname_back,
                ref_fname='data/fiber_arm_outline.fits',
                ang_step=1,
                doplot=False):
    """ Process the front/back iluminated images
    ang_step step in degrees; The smaller the step the slower the code.
    Returns: the tuple with
    1) table of detections with x,y,angle columns
    2) the extracted rotated canny images of all detected fibers
    """
    thresh = 120  # threshold for detecting fibers
    peakv = 10000  # minimum value of the peak corresponding to illuminated fiber
    pad = 40  # padding around the fiber to analyze
    # changing padding would require regenerating the template image
    minccf = 0.2  # minimum cross-correlation between the template and data to remove non-fibers/fiducials

    D = pyfits.getdata(fname_front)
    dat = pyfits.getdata(fname_back)
    dat = dat - np.median(dat)
    det = sep.extract(dat, thresh)
    cans = []
    ref = pyfits.getdata(ref_fname)
    det = atpy.Table(det)
    det['angle'] = np.nan
    det['ccf'] = np.nan
    xgrid, ygrid = np.mgrid[-pad:pad + 1, -pad:pad + 1]
    mask = np.sqrt(xgrid**2 + ygrid**2) < pad
    for curi in range(len(det['x'])):
        if det['peak'][curi] < peakv:
            continue
        curx, cury = int(np.round(det['y'][curi])), int(
            np.round(det['x'][curi]))
        curD = D[curx - pad:curx + pad + 1, cury - pad:cury + pad + 1]
        curd = dat[curx - pad:curx + pad + 1, cury - pad:cury + pad + 1]
        curD = norm(curD)
        can = canny(curD, sigma=2) * mask
        #cans.append(can)
        angle, ccfval = find_angle(can, ref, ang_step=ang_step)
        cans.append(rotate(can * 1., angle))
        angle = np.deg2rad(angle - 90)
        det['angle'][curi] = angle
        det['ccf'][curi] = ccfval
        if doplot:
            plt.clf()
            plt.subplot(131)
            Rs = np.r_[0, 30]
            plt.imshow(curD, vmax=perc(curD, .95))
            plt.plot(pad + np.cos(angle) * Rs,
                     pad + np.sin(angle) * Rs,
                     color='red')
            plt.subplot(132)
            plt.imshow(curd, vmax=perc(curd, .99))
            plt.subplot(133)
            plt.imshow(can)
    det = det[np.isfinite(det['angle']) & (det['ccf'] > minccf)]
    cans = np.array(cans)
    return det, cans


def find_fibers_front(fname_front,
                      ref_fname='data/fiber_arm_outline.fits',
                      ang_step=1,
                      doplot=False):
    """ (experimental) 
    Process the front iluminated images. return the matched filtered image
    ang_step step in degrees; The smaller the step the slower the code.
    """

    pad = 40  # padding around the fiber to analyze
    # changing padding would require regenerating the template image
    minccf = 0.2  # minimum cross-correlation between the template and data to remove non-fibers/fiducials

    D = pyfits.getdata(fname_front)
    nx, ny = D.shape
    block = 300
    D1 = D * 0.
    sigma = 2
    xgrid, ygrid = np.mgrid[-pad:pad + 1, -pad:pad + 1]
    mask = (np.sqrt(xgrid**2 + ygrid**2) < pad).astype(int)
    #mask = mask/(mask**2).sum()**.5
    for i in np.arange(0, nx, block):
        for j in np.arange(0, ny, block):
            curd = D[i:i+block, j :j + block]
            curd = norm(curd)
            D1[i:i + block,
               j:j + block] = canny(curd, sigma=sigma)
    ref = pyfits.getdata(ref_fname).astype('<f8')
    refs = []
    for angle in np.arange(0, 360, ang_step):
        curr = rotate(ref, angle)
        curr = curr/(curr**2).sum()**.5
        refs.append(curr)

    dats=[]
    allconv = None
    im0 = scipy.signal.fftconvolve(D1**2, mask, mode='same')
    for r in refs:
        curt = scipy.signal.convolve(D1, r, mode='same')
        if allconv is None:
            allconv = curt
        else:
            allconv = np.maximum(allconv, curt)
    
    return allconv/im0**.5
