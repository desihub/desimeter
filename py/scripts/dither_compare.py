import os

import pylab as plt
import numpy as np

import fitsio
import astropy
from astropy.table import Table

from desimeter.desimeter import Desimeter
from desimeter.fieldmodel import FieldModel
from desimeter.log import get_logger
from desimeter.transform.xy2qs import qs2xy
from desimeter.match import match_same_system
from desimeter.transform.radec2tan import hadec2xy
from desimeter.simplecorr import SimpleCorr

def main():
    plt.figure(figsize=(6,6))
    log = get_logger()

    dm = Desimeter(desimeter_dir='dm-fid-sys', proc_data_dir='proc-dither')

    dither_file = 'dither20200315-63224-B.fits'

    tile = 63225
    expnum = 55670
    frame = 1

    max_match_distance = 7.

    fn = dm.find_file('fieldmodel', expnum=expnum, frame=frame)
    if not os.path.exists(fn):
        fm = dm.fit_guide_star_coords(expnum, frame)
        with open(fn, 'w') as f:
            f.write(fm.tojson())
    else:
        with open(fn) as f:
            fm = FieldModel.fromjson(f.read())

    fn = dm.find_file('fvc-spots', expnum=expnum, frame=frame)
    if os.path.exists(fn):
        spots = astropy.table.Table.read(fn)
    else:
        spots = dm.measure_spots(expnum, frame)
        spots.write(fn, overwrite=True)

    spots["RA"] = np.zeros(len(spots),dtype=float)
    spots["DEC"] = np.zeros(len(spots),dtype=float)
    ii=(spots["X_FP"]!=0)&(spots["Y_FP"]!=0)
    ra,dec = fm.fp2radec(spots["X_FP"][ii],spots["Y_FP"][ii])
    spots["RA"][ii]  = ra
    spots["DEC"][ii] = dec

    coordfn = dm.find_file('coordinates', expnum=expnum)

    expected_pos = Table.read(coordfn) # auto-magically guess format
    if not "X_FP_EXP" in expected_pos.keys() :
        if "EXP_Q_0" in  expected_pos.keys() :
            log.info("EXP_Q_0,EXP_S_0 -> X_FP,Y_FP")
            x,y = qs2xy(q=expected_pos["EXP_Q_0"],s=expected_pos["EXP_S_0"])
            expected_pos["X_FP"] = x
            expected_pos["Y_FP"] = y
            bad_expected_pos = (np.isnan(x)|np.isnan(y))
            if np.sum(bad_expected_pos) > 0 :
                expected_pos["X_FP"][bad_expected_pos]=-99999.
                expected_pos["Y_FP"][bad_expected_pos]=-99999.
        else :
            log.error("No EXP_Q_0 nor X_FP in expected positions file {}".format(coordfn))

    if not "LOCATION" in expected_pos.keys() :
        # add useful location keyword
        expected_pos["LOCATION"] = (np.array(expected_pos["PETAL_LOC"])*1000 +
                                    np.array(expected_pos["DEVICE_LOC"]))

    if "PINHOLE_ID" in expected_pos.dtype.names :
        # exclude pinhole because here we want to match fibers
        ii = np.where(expected_pos["PINHOLE_ID"]==0)[0]
        expected_pos = expected_pos[:][ii]

    # select spots that are not already matched
    selection  = (spots["LOCATION"]==-1)

    # match
    indices_of_expected_pos,distances = match_same_system(spots["X_FP"][selection],spots["Y_FP"][selection],expected_pos["X_FP"],expected_pos["Y_FP"])

    is_matched = (distances<max_match_distance)&(indices_of_expected_pos>=0)
    ii=np.where(selection)[0]
    selection[ii]          &=  is_matched
    indices_of_expected_pos = indices_of_expected_pos[is_matched]
    distances               = distances[is_matched]

    # add columns after matching fibers
    for k1,k2 in zip(["X_FP","Y_FP"],["X_FP_EXP","Y_FP_EXP"]) :
        if k2 not in spots.keys() : spots[k2] = np.zeros(len(spots))
        spots[k2][selection]=expected_pos[k1][indices_of_expected_pos]
    for k in ["EXP_Q_0","EXP_S_0","PETAL_LOC","DEVICE_LOC","DEVICE_ID","DEVICE_TYPE","LOCATION"] :
        if k in expected_pos.keys() :
            if k not in spots.keys() :
                if k in ["DEVICE_ID","DEVICE_TYPE"] :
                    spots[k] = np.repeat("None           ",len(spots))
                else :
                    spots[k] = np.zeros(len(spots))
            spots[k][selection]=expected_pos[k][indices_of_expected_pos]

    # for spots with metrology X_FP_EXP=X_FP_METRO
    selection = (spots["X_FP_METRO"]!=0)
    spots["X_FP_EXP"][selection]=spots["X_FP_METRO"][selection]
    selection = (spots["Y_FP_METRO"]!=0)
    spots["Y_FP_EXP"][selection]=spots["Y_FP_METRO"][selection]

    from fiducial_systematics import plot_fiducial_offsets

    plt.clf()
    plot_fiducial_offsets(spots, expnum=expnum, frame=frame)
    plt.savefig('fvc-dither-%08i-F%04i.png' % (expnum, frame))


    # plot_dither_vs_desimeter

    fiberassign_file = dm.find_file('fiberassign', expnum=expnum, tile=tile)

    fvc = spots
    ii=(fvc["RA"]!=0)
    fvc=fvc[ii]

    # Eddie's dithers
    #######################################################################
    t=fitsio.read(dither_file)

    expid = expnum
    ind = np.flatnonzero(t['expid'][0, :] == expid)
    if len(ind) > 0:
        i = ind[0]
    else:
        raise ValueError(f'could not find expid {expid} in dither file.')

    err=np.sqrt(t['dxfiboff'][:,i]**2+t['dyfiboff'][:,i]**2)
    roff=np.sqrt(t['xfiboff'][:,i]**2+t['yfiboff'][:,i]**2)
    jj=np.where((err<0.5)&(roff<4.))[0]

    dither_ra    = t['fiber_ditherfit_ra'][jj,i]
    dither_dec   = t['fiber_ditherfit_dec'][jj,i]
    dither_fiber = t['fiber'][jj,i].astype(int)
    print("number of valid targets from dither analysis = {}".format(len(dither_ra)))

    # fiberassign
    #######################################################################
    fiberassign=Table.read(fiberassign_file,1)
    head = fitsio.read_header(fiberassign_file)
    dico = { loc : i for i,loc in enumerate(fiberassign["LOCATION"])}

    fvc_indices=[]
    fiberassign_indices=[]
    for i,loc in enumerate(fvc["LOCATION"]) :
        if loc in dico :
            fvc_indices.append(i)
            fiberassign_indices.append(dico[loc])

    fvc=fvc[fvc_indices]
    fiberassign=fiberassign[fiberassign_indices]

    # re-address the dither results to map desimeter table which has location
    #######################################################################
    dico = { fiber : i for i,fiber in enumerate(dither_fiber) }

    fvc_indices=[]
    dither_indices=[]
    for i,fiber in enumerate(fiberassign["FIBER"]):
        if fiber in dico :
            fvc_indices.append(i)
            dither_indices.append(dico[fiber])

    fvc=fvc[fvc_indices]
    fiberassign=fiberassign[fvc_indices]
    dither_ra=dither_ra[dither_indices]
    dither_dec=dither_dec[dither_indices]

    tel_ra=head["REQRA"]
    tel_dec=head["REQDEC"]
    desimeter_x,desimeter_y=hadec2xy(-fvc["RA"]+tel_ra,fvc["DEC"],0,tel_dec)
    target_x,target_y=hadec2xy(-fiberassign["TARGET_RA"]+tel_ra,fiberassign["TARGET_DEC"],0,tel_dec)
    dither_x,dither_y=hadec2xy(-dither_ra+tel_ra,dither_dec,0,tel_dec)

    dd = np.hypot(desimeter_x - target_x, desimeter_y - target_y)
    print('Desimeter-target:', len(desimeter_x), 'median:', np.median(dd), 'max', np.max(dd))
    print('99th pct:', np.percentile(dd, 99))
    # ~ 1"
    I = np.flatnonzero(dd < 5e-6)
    desimeter_x = desimeter_x[I]
    desimeter_y = desimeter_y[I]
    target_x = target_x[I]
    target_y = target_y[I]
    dither_x = dither_x[I]
    dither_y = dither_y[I]

    # fit a transform to adjust desimeter to the dither results
    corr=SimpleCorr()
    corr.fit(desimeter_x,desimeter_y,dither_x,dither_y)
    desimeter_x_bis , desimeter_y_bis = corr.apply(desimeter_x,desimeter_y)

    residual_ha_offset_rad     = corr.dx
    residual_ha_offset_arcsec  = corr.dx*180*3600./np.pi
    residual_dec_offset_rad    = corr.dy
    residual_dec_offset_arcsec = corr.dy*180*3600./np.pi
    residual_rotation_deg      = corr.rot_deg -360*(corr.rot_deg>180)
    residual_rotation_arcsec   =  residual_rotation_deg*3600.

    print("residual scale along HA       = 1{:+4.2g}".format(corr.sxx-1.))
    print("residual scale along Dec      = 1{:+4.2g}".format(corr.syy-1.))
    print("residual pointing offset dHA  = {:4.2f} arcsec".format(residual_ha_offset_arcsec))
    print("residual pointing offset dDec = {:4.2f} arcsec".format(residual_dec_offset_arcsec))
    print("residual rotation             = {:4.2f} arcsec".format(residual_rotation_arcsec))

    # from radian to degree
    scale=180/np.pi
    desimeter_x *= scale
    desimeter_y *= scale
    target_x *= scale
    target_y *= scale
    dither_x *= scale
    dither_y *= scale
    desimeter_x_bis *= scale
    desimeter_y_bis *= scale

    def myquiver(x1,y1,x2,y2,title=None):
        x0=-1.5 ; y0=1.7
        dx=(x2-x1)
        dy=(y2-y1)
        Q = plt.quiver(x1,y1,dx,dy)
        dr=np.sqrt(dx[:-1]**2+dy[:-1]**2)
        rms2d=np.sqrt(np.mean(dr**2))
        darcsec = 0.5
        plt.quiverkey(Q, x0, y0, darcsec, '%f arcsec' % darcsec)
        text="rms$_{2D}$ = %3.2f''"%(rms2d*3600.)
        plt.text(1.6,y0,text,fontsize=10,
                 horizontalalignment="right",verticalalignment="center")
        text="rms 2D = %3.2f arcsec"%(rms2d*3600.)
        if title:
            print("%30s %s"%(title,text))
            plt.title('%s %s' % (title, text))
        else:
            print(text)

    title="desimeter(guide+fvc)-target"
    plt.clf()
    myquiver(target_x,target_y,desimeter_x,desimeter_y,title=title)
    plt.savefig('dither-%08i-1.png' % expnum)

    title="dither-target"
    plt.clf()
    myquiver(target_x,target_y,dither_x,dither_y,title=title)
    plt.savefig('dither-%08i-2.png' % expnum)

    title="desimeter(guide+fvc)-dither"
    plt.clf()
    myquiver(dither_x,dither_y,desimeter_x,desimeter_y,title=title)
    plt.savefig('dither-%08i-3.png' % expnum)

    title="transformed(desimeter) -dither"
    plt.clf()
    myquiver(dither_x,dither_y,desimeter_x_bis,desimeter_y_bis,title=title)
    plt.savefig('dither-%08i-4.png' % expnum)

if __name__ == '__main__':
    main()
