
"""
Utility functions to detect spots in fiber view camera image
"""

import numpy as np
from desiutil.log import get_logger
from astropy.table import Table
from scipy.signal import fftconvolve

def gaussian_convolve(image,sigma=1.) :
    """FFT convolution of an image with a gaussian 2d kernel
     Args:
        image : 2D numpy array

     Optional:
        sigma : float, default=1

     returns convolved image, same shape as input
    """
    
    hw=int(3*sigma)
    x2d=np.tile(np.arange(-hw,hw+1),(2*hw+1,1))
    kernel = np.exp(-(x2d**2+x2d.T**2)/2./sigma**2)
    kernel /= np.sum(kernel)
    return fftconvolve(image,kernel,"same")

def fitcentroid(stamp) :
    """
    fit the centroid of a 2D image square stamp of size (2*hw+1,2*hw+1)
    in a coordinate system centered on the stamp, i.e. the central pixel has coordinates = (0,0)
    """
    hw=stamp.shape[0]//2
    x1d=np.arange(-hw,hw+1)
    
    # for the moment it's a simplistic barycenter
    # one can do much better than this
    norm=np.sum(stamp)
    if norm <= 0 :
        raise RuntimeError("sum of flux in stamp = {}".format(norm))
    
    xc=np.sum(x1d*np.sum(stamp,axis=0))/norm
    yc=np.sum(x1d*np.sum(stamp,axis=1))/norm
    return xc,yc,norm

def detectspots(fvcimage,threshold=None) :
    """
    Detect spots in a fiber view image and measure their centroids and flux
    Args:
        fvcimage : 2D numpy array

    Optional:
        threshold : float, use this threshold in counts/pixel to do a first
                    detection of peaks
    returns astropy.Table with spots, columns are xpix,ypix,xerr,yerr,counts
    """

    log = get_logger()

    
    n0=fvcimage.shape[0]
    n1=fvcimage.shape[1]
    
    

    # find peaks = local maximum above threshold
    log.info("gaussian convolve...")
    convolved_image=gaussian_convolve(fvcimage)
    log.info("done")

    
    # measure pedestal and rms
    # look the values of a random subsample of the image
    nrand=20000
    ii0=(np.random.uniform(size=nrand)*n0).astype(int)
    ii1=(np.random.uniform(size=nrand)*n1).astype(int)
    vals=convolved_image[ii0,ii1].ravel()
    mval=np.median(vals)
    #- normalized median absolute deviation as robust version of RMS
    #- see https://en.wikipedia.org/wiki/Median_absolute_deviation
    rms=1.4826*np.median(np.abs(vals-mval)) 
    ok=np.abs(vals-mval)<4*rms
    mval=np.mean(vals[ok])
    rms=np.std(vals[ok])
    log.info("pedestal={:4.2f} rms={:4.2f}".format(mval,rms))
    # remove mean
    fvcimage -= mval
    convolved_image      -= mval
    
    #import matplotlib.pyplot as plt
    #plt.hist(vals[ok]-mval,bins=100)
    #plt.show()
    
    if threshold is None :
        threshold=7*rms
    
    
    peaks=np.zeros((n0,n1))
    peaks[1:-1,1:-1]=(convolved_image[1:-1,1:-1]>convolved_image[:-2,1:-1])*(convolved_image[1:-1,1:-1]>convolved_image[2:,1:-1\
])*(convolved_image[1:-1,1:-1]>convolved_image[1:-1,:-2])*(convolved_image[1:-1,1:-1]>convolved_image[1:-1,2:])*(convolved_image[1:-1,1:-1]>threshold)

    # loop on peaks
    peakindices=np.where(peaks.ravel()>0)[0]
    npeak=len(peakindices)
    if npeak == 0 :
        log.error("no spot found")
        raise RuntimeError("no spot found")
    else :
        log.info("found {} peak".format(npeak))
                 
    xpix=np.zeros(npeak)
    ypix=np.zeros(npeak)
    counts=np.zeros(npeak)
    hw=3
    for j,index in enumerate(peakindices) :
        i0=index//n1
        i1=index%n1
        try :
            x,y,c=fitcentroid(fvcimage[i0-hw:i0+hw+1,i1-hw:i1+hw+1])
            x += i1 # x is along axis=1 in python
            y += i0 # y is along axis=1 in python
            xpix[j] = x
            ypix[j] = y
            counts[j] = c
        except Exception as e:
            log.error("failed to fit a centroid {}".format(e))
            
        log.debug("{} x={} y={} counts={}".format(j,xpix[j],ypix[j],counts[j]))
    
    log.warning("NOT IMPLEMENTED: would need some cleaning here for multiple detections of same spot")

    table = Table([xpix,ypix,counts],names=("XPIX","YPIX","COUNTS"))
    return table

    
    
    
