import pylab as plt
import numpy as np
import os
from glob import glob

from astrometry.util.fits import fits_table, merge_tables

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('date', type=str, help='MMDD date string')
opts = parser.parse_args()
date = opts.date
assert(len(date) == 4)
MM = date[:2]
DD = date[2:]

# 2020-02-24
datestr = '2020-'+MM+'-'+DD
plotdir = MM+DD+'/'
fns = glob(plotdir + 'fm*.fits')

# 2020-03-15
# datestr = '2020-03-15'
# plotdir = '0315/'
# fns = glob('0315/fm*.fits')

fns.sort()
print(len(fns), 'files')

TT = []
for fn in fns:
    T = fits_table(fn)
    TT.append(T)
T = merge_tables(TT)

# FIXME
cosdec = np.cos(np.deg2rad(T.dec_gaia))
T.dist = np.hypot((T.ra_gaia - T.fit_ra)*cosdec, T.dec_gaia - T.fit_dec)
T.dist *= 3600.

cams = np.unique(T.camera)
for cam in cams:
    I = np.flatnonzero(T.camera == cam)
    print(len(I), 'in', cam)
    #I = np.flatnonzero((T.camera == cam) * (T.dist < 0.5))
    #print(len(I), 'in', cam, 'w match < 0.5"')
    TI = T[I]
    plt.clf()
    plt.quiver(TI.xcentroid, TI.ycentroid,
               TI.fit_gaia_x - TI.xcentroid,
               TI.fit_gaia_y - TI.ycentroid,
               TI.mjd_obs,
               label='Measurements - Gaia',
               angles='xy', scale_units='xy', scale=0.01)
    # GFA physical CCD size
    H,W = 1032,2048
    xx = np.array([1, 1, W, W, 1])
    yy = np.array([1, H, H, 1, 1])
    plt.plot(xx, yy, 'k-')
    plt.axis('equal')
    m = 20
    plt.axis([1-m, W+m, 1-m, H+m])
    plt.title('%s: Camera %s' % (datestr, cam))
    #plt.title('2020-03-15: Camera %s' % cam)
    plt.savefig(plotdir + 'summary-1-%s.png' % cam)


cams = np.unique(T.camera)
expids = np.unique(T.expid)

camradecs = {}
for cam in cams:
    I = np.flatnonzero((T.camera == cam))
    TI = T[I]
    ## FIXME RA wrap-around
    camra = np.mean((TI.fit_ra - TI.target_ra) * np.cos(np.deg2rad(TI.target_dec)))
    camdec = np.mean((TI.fit_dec - TI.target_dec))
    camradecs[cam] = (camra, camdec)

dr,dd = [],[]
cr,cd = [],[]
ee = []

for expid in expids:
    for cam in cams:
        I = np.flatnonzero((T.camera == cam) * (T.expid == expid))
        print(len(I), 'in', cam, 'and expid', expid)

        TI = T[I]

        dra  = (TI.fit_ra  - TI.ra_gaia) * np.cos(np.deg2rad(TI.fit_dec))
        ddec = (TI.fit_dec - TI.dec_gaia)
        
        dr.append(np.mean(dra))
        dd.append(np.mean(ddec))
        r,d = camradecs[cam]
        cr.append(r)
        cd.append(d)
        ee.append(expid)

plt.clf()
theta = np.linspace(0, 2.*np.pi, 100)
r = 1.6
plt.plot(r*np.cos(theta), r*np.sin(theta), 'k-')
plt.quiver(cr, cd, dr, dd, ee,
           label='Measurements - Gaia',
           angles='xy', scale_units='xy', scale=0.0001)
#plt.axis('equal')
#plt.xlim(-2.5, 2.5)
plt.axis('square')
plt.axis([-2.5, 2.5, -2.5, 2.5])
plt.title('%s' % datestr)
plt.savefig(plotdir + 'summary-2.png')
