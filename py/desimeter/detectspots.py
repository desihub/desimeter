
"""
Utility functions to detect spots in fiber view camera image
"""
import sys
import numpy as np
from desiutil.log import get_logger
from astropy.table import Table
from scipy.signal import fftconvolve
from scipy.special import erf

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

def fitcentroid_barycenter(stamp) :
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

    # uncertainties, made up for now ... 
    xe=1.
    ye=1.
    
    return xc,yc,xe,ye,norm


def psf(i0,i1,xc,yc,sigma) :
    a=1/(np.sqrt(2)*sigma)
    return 0.25*(erf(a*(i1+0.5-xc))-erf(a*(i1-0.5-xc)))*(erf(a*(i0+0.5-yc))-erf(a*(i0-0.5-yc)))

def dpsfdxc(i0,i1,xc,yc,sigma) :
    a=1/(np.sqrt(2)*sigma)
    return -a*0.25*2/np.sqrt(np.pi)*(np.exp(-(a*(i1+0.5-xc))**2)-np.exp(-(a*(i1-0.5-xc))**2))*(erf(a*(i0+0.5-yc))-erf(a*(i0-0.5-yc)))

def dpsfdyc(i0,i1,xc,yc,sigma) :
    a=1/(np.sqrt(2)*sigma)
    return -a*0.25*2/np.sqrt(np.pi)*(erf(a*(i1+0.5-xc))-erf(a*(i1-0.5-xc)))*(np.exp(-(a*(i0+0.5-yc))**2)-np.exp(-(a*(i0-0.5-yc))**2))

def fitcentroid_gaussian(stamp,sigma=1.,noise=10.) :
    """
    fit the centroid of a 2D image square stamp of size (2*hw+1,2*hw+1)
    in a coordinate system centered on the stamp, i.e. the central pixel has coordinates = (0,0)
    """
    # iterative gauss-newton fit
    n0=stamp.shape[0]
    n1=stamp.shape[1]
    xc=float(n1//2)
    yc=float(n0//2)
    flux=np.sum(stamp)
    #print(flux,xc,yc)
    mod    = np.zeros(stamp.shape)
    dmoddx = np.zeros(stamp.shape)
    dmoddy = np.zeros(stamp.shape)
    # pixel indices (in 2D)
    ii0 = np.tile(np.arange(n0),(n1,1)).T
    ii1 = np.tile(np.arange(n1),(n0,1))
    
    for loop in range(5) :
        mod    = psf(ii0,ii1,xc,yc,sigma)
        dmoddx = dpsfdxc(ii0,ii1,xc,yc,sigma)
        dmoddy = dpsfdyc(ii0,ii1,xc,yc,sigma)
        H=np.array([mod,flux*dmoddx,flux*dmoddy]).reshape(3,n0*n1)
        B=((stamp-flux*mod).reshape(n0*n1)*H).sum(axis=1)
        A=H.dot(H.T)
        Ai=np.linalg.inv(A)
        delta=Ai.dot(B)
        val=np.max(np.abs(delta[1:]))
        if val>0.2 :
            delta *= (0.2/val) # limiting range
            
        flux += delta[0]
        xc += delta[1]
        yc += delta[2]
        #print(loop,delta,[flux,xc,yc])
        if np.abs(delta[1])<0.001 and np.abs(delta[2])<0.001 : break

    # coordinates should be zero at center of stamp
    xc -= float(n1//2)
    yc -= float(n0//2)
    
    xe = noise * np.sqrt(Ai[1,1])
    ye = noise * np.sqrt(Ai[2,2])

    return xc,yc,xe,ye,flux
    #sys.exit(12)

    
def fitcentroid(stamp,noise=10.) :
    """
    fit the centroid of a 2D image square stamp of size (2*hw+1,2*hw+1)
    in a coordinate system centered on the stamp, i.e. the central pixel has coordinates = (0,0)
    """

    return fitcentroid_gaussian(stamp,sigma=1,noise=noise)
    #return fitcentroid_barycenter(stamp)


def detectspots(fvcimage,threshold=500,nsig=7,psf_sigma=1.) :
    """
    Detect spots in a fiber view image and measure their centroids and flux
    Args:
        fvcimage : 2D numpy array

    Optional:
       threshold : float, use max of this threshold in counts/pixel and nsig*rms to do a first
                    detection of peaks 
    returns astropy.Table with spots, columns are xpix,ypix,xerr,yerr,counts
    """

    log = get_logger()

    
    n0=fvcimage.shape[0]
    n1=fvcimage.shape[1]
    
    

    # find peaks = local maximum above threshold
    log.info("gaussian convolve with sigma = {:2.1f} pixels".format(psf_sigma))
    convolved_image=gaussian_convolve(fvcimage,psf_sigma)
    
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
        threshold=nsig*rms
    else :
        threhold=max(threshold,nsig*rms)
    
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
        log.info("found {} peaks".format(npeak))
                 
    xpix=np.zeros(npeak)
    ypix=np.zeros(npeak)
    xerr=np.zeros(npeak)
    yerr=np.zeros(npeak)
    counts=np.zeros(npeak)
    hw=3
    
    xoffset=0. # would need offsets of 1 to match with POS file. not sure what is the best choice here.
    yoffset=0.
    if xoffset !=0 or yoffset !=0 :
        log.warning("Applying offsets x += {} and y += {} (to match with others, like the POS files)".format(xoffset,yoffset))
    if xoffset == 0 and yoffset == 0 :
        log.warning("Here center of first pixel has coord=(0,0); so we expect offsets of 1 with respect to coordinates in .pos files.")
    
    for j,index in enumerate(peakindices) :
        i0=index//n1
        i1=index%n1
        if 1 : #try :
            x,y,ex,ey,c=fitcentroid(fvcimage[i0-hw:i0+hw+1,i1-hw:i1+hw+1],noise=rms)
            xpix[j] = x + i1 + xoffset #  x is along axis=1 in python , also adding offset (definition of pixel coordinates)
            ypix[j] = y + i0 + yoffset #  y is along axis=0 in python , also adding offset (definition of pixel coordinates)
            xerr[j] = ex
            yerr[j] = ey
            counts[j] = c
        #except Exception as e:
        #    log.error("failed to fit a centroid {}".format(e))
            
        log.debug("{} x={} y={} counts={}".format(j,xpix[j],ypix[j],counts[j]))
    
    #log.warning("Would need some cleaning here for multiple detections of same spot")

    table = Table([xpix,ypix,xerr,yerr,counts],names=("XPIX","YPIX","XERR","YERR","COUNTS"))
    return table

    
    
    
