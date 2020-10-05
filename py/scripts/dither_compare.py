import os
import pylab as plt
import numpy as np

import astropy
from astropy.table import Table

from desimeter.desimeter import Desimeter
from desimeter.fieldmodel import FieldModel
from desimeter.log import get_logger
from desimeter.transform.xy2qs import qs2xy
from desimeter.match import match_same_system

def main():
    plt.figure(figsize=(10,10))
    log = get_logger()

    dm = Desimeter(desimeter_dir='dm-fid-sys', proc_data_dir='proc-dither')

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


if __name__ == '__main__':
    main()
