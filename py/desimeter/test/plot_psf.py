#!/usr/bin/env python


import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf


def psf(i0,i1,xc,yc,sigma) :
    a=1/(np.sqrt(2)*sigma)
    return 0.25*(erf(a*(i1+0.5-xc))-erf(a*(i1-0.5-xc)))*(erf(a*(i0+0.5-yc))-erf(a*(i0-0.5-yc)))
def dpsfdxc(i0,i1,xc,yc,sigma) :
    a=1/(np.sqrt(2)*sigma)
    return -a*0.25*2/np.sqrt(np.pi)*(np.exp(-(a*(i1+0.5-xc))**2)-np.exp(-(a*(i1-0.5-xc))**2))*(erf(a*(i0+0.5-yc))-erf(a*(i0-0.5-yc)))
def dpsfdyc(i0,i1,xc,yc,sigma) :
    a=1/(np.sqrt(2)*sigma)
    return -a*0.25*2/np.sqrt(np.pi)*(erf(a*(i1+0.5-xc))-erf(a*(i1-0.5-xc)))*(np.exp(-(a*(i0+0.5-yc))**2)-np.exp(-(a*(i0-0.5-yc))**2))


sigma=0.8
xc=14.2
yc=5.7


if 1 :

    img=np.zeros((30,30))
    img2=np.zeros((30,30))
    
    ii0 = np.tile(np.arange(img.shape[0]),(img.shape[1],1)).T
    ii1 = np.tile(np.arange(img.shape[1]),(img.shape[0],1))
    eps=0.0001
    img  = (psf(ii0,ii1,xc+eps/2,yc,sigma)-psf(ii0,ii1,xc-eps/2,yc,sigma))/eps
    img2 = dpsfdxc(ii0,ii1,xc,yc,sigma)

    plt.figure("dpsfdxc")
    j=np.argmax(np.sum(np.abs(img),axis=1))
    plt.plot(img[j,:],label="numeric derivative")
    plt.plot(img2[j,:],"--",label="analytic derivative")
    plt.legend()
    
    #plt.figure("dpsfdxc1")
    #plt.imshow(img,origin=0)
    #plt.figure("dpsfdxc2")
    #plt.imshow(img2,origin=0)
     
if 1 :
    ii1 = np.tile(np.arange(img.shape[1]),(img.shape[0],1))
    ii0 = np.tile(np.arange(img.shape[0]),(img.shape[1],1)).T
    img = psf(ii0,ii1,xc,yc,sigma)
    plt.figure("img")
    plt.imshow(img,origin=0)
    
plt.show()
    
