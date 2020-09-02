"""
Test spotmatch interface with desimeter
"""

import os
import unittest
import shutil
import numpy as np
from astropy.table import Table

from desimeter.spotmatch import spotmatch
from desimeter.io import load_metrology,fvc2fp_filename
from desimeter.transform.fvc2fp import FVC2FP

class TestSpotmatch(unittest.TestCase):

    def test_spotmatch(self):

        if shutil.which("match_positions") is None :
            print("cannot test spotmatch because match_positions is not in PATH.")
            return

        measured_spots=None
        expected_spots=None
        print("testing spotmatch")

        fvc2fp = FVC2FP.read(fvc2fp_filename())
        metrology = load_metrology()
        spots = metrology[(metrology["DEVICE_TYPE"]=="POS")|(metrology["DEVICE_TYPE"]=="FIF")|(metrology["DEVICE_TYPE"]=="GIF")]
        xfp = spots["X_FP"]
        yfp = spots["Y_FP"]

        selection = (spots["DEVICE_TYPE"]=="POS")

        # pos move error
        xfp[selection] += 0.01 * np.random.normal(size=xfp[selection].size)
        yfp[selection] += 0.01 * np.random.normal(size=yfp[selection].size)

        # measurement error
        xfp += 0.005 * np.random.normal(size=xfp.size)
        yfp += 0.005 * np.random.normal(size=yfp.size)

        xpix,ypix = fvc2fp.fp2fvc(xfp,yfp)

        print("test: number of spots sent to spotmatch=",len(xpix))

        res = spotmatch(xpix,ypix,expected_x_fp=None,expected_y_fp=None,expected_location=None)
        #res.write("res.csv",overwrite=True)

if __name__ == '__main__':
    unittest.main()
