#!/usr/bin/env python

import numpy as np
#import matplotlib.pyplot as plt
from astropy.table import Table,Column # for the IO
import RT185v2    # needs DESI-ADC-2.OPT.CSV  DESI-ADC-2.RAY.CSV  DESI-ADC-2.MED.CSV
# These files give us up to 680 rays, 8 wavelengths.
from desimeter.trig import sind

nthread = 19

refwave = np.array([3650,   4360,   4860,   5880,   6380,   6560,   8520,  10140]) #A

def incoming_rays(adc1=0.,adc2=0., dangle=0.2) :
    # w = wavelength number 0 through 7 corresponding to the arrays in line 976:
    # u = incoming ray direction cosines: u=+x component of unit length vector = eastward ray from western target
    # v = incoming ray direction cosine:   v = +y component of  unit length vector southward from northern target
    # adc1 angle in degrees
    # adc2 angle in degrees

    wuv12 = [] # list of arrays with 5 components [w,u,v,adc1,adc2]

    ntheta=int(1.6/dangle)+1
    for theta_deg in np.linspace(0,1.6,ntheta) :
        st=sind(theta_deg)
        nphi = int((2.*np.pi*theta_deg)/dangle)+1
        phi_rad = np.arange(nphi)*(2.*np.pi/nphi)
        uu = st*np.cos(phi_rad)
        vv = st*np.sin(phi_rad)
        w = 2 # I like 486nm
        for u,v in zip(uu,vv) :
            wuv12.append([w,u,v,adc1,adc2])
    wuv12=np.vstack(wuv12)
    return wuv12

def trace(wuv12):                # TASK 10
    resultsList = []
    arrayNine = RT185v2.getNine(wuv12)          # one star, but either monochromatic or polychromatic
    resultsList.append(arrayNine)
    return resultsList  # adc1, adc2, ngood, xave, yave, zave, xrms, yrms, zrms


def work(adc2):
    result = []
    adc1 = -adc2
    rays=incoming_rays(adc1=adc1,adc2=adc2,dangle=0.1)
    xfp = np.zeros((rays.shape[0]))
    yfp = np.zeros((rays.shape[0]))
    zfp = np.zeros((rays.shape[0]))
    for i,ray in enumerate(rays) :
        adc1, adc2, ngood, xfp[i] , yfp[i], zfp[i], xrms, yrms, _ = RT185v2.getNine(ray)

    wave = refwave[ rays[:,0].astype(int) ]
    xtan = rays[:,1]
    ytan = rays[:,2]
    for i in range(len(rays)) :
        result.append([adc1,adc2,xtan[i],ytan[i],xfp[i],yfp[i],zfp[i],wave[i]])
    return result


def main() :


    result=[]
    from multiprocessing.pool import Pool
    pool = Pool(nthread)
    result = pool.map(work, np.linspace(0., 180., 19))
    result = np.concatenate(result)
    print(result.shape)

    # make a table
    nray=len(result)
    table = Table()
    table.add_column(Column(name="ADC1",dtype=float,unit="deg",length=nray))
    table.add_column(Column(name="ADC2",dtype=float,unit="deg",length=nray))
    table.add_column(Column(name="X_TAN",dtype=float,unit="none",length=nray))
    table.add_column(Column(name="Y_TAN",dtype=float,unit="none",length=nray))
    table.add_column(Column(name="X_FP",dtype=float,unit="mm",length=nray))
    table.add_column(Column(name="Y_FP",dtype=float,unit="mm",length=nray))
    table.add_column(Column(name="Z_FP",dtype=float,unit="mm",length=nray))
    table.add_column(Column(name="WAVELENGTH",dtype=float,unit="A",length=nray))

    table["ADC1"] = result[:,0]
    table["ADC2"] = result[:,1]
    table["X_TAN"] = result[:,2]
    table["Y_TAN"] = result[:,3]
    table["X_FP"] = result[:,4]
    table["Y_FP"] = result[:,5]
    table["Z_FP"] = result[:,6]
    table["WAVELENGTH"] = result[:,7]

    table.write("raytrace-tan2fp-4957-v18.csv")

main()
