import unittest
from pkg_resources import resource_filename
from collections import Counter

import numpy as np
from desimeter.fiberassign import fiberassign_radec2xy_cs5 , fiberassign_cs5_xy2radec, fiberassign_radec2xy_flat , fiberassign_flat_xy2radec
from desimeter.transform.radec2tan import xy2hadec,hadec2xy
from desimeter.trig import cosd

class TestFA(unittest.TestCase):

     def test_it(self):
        print("Testing fiberassign_radec2xy_cs5")
        nrand = 500
        square_side_in_degrees =  1.5
        xtan = np.pi/180. * square_side_in_degrees * (np.random.uniform(size=nrand)-0.5)
        ytan = np.pi/180. * square_side_in_degrees * (np.random.uniform(size=nrand)-0.5)
        tile_ra = 12.
        tile_dec = 40.
        tile_ha  = 1
        tile_mjd = 59344.4
        tile_fieldrot = -130/3600. # deg
        lst = tile_ra+tile_ha

        ha,dec = xy2hadec(xtan,ytan,tile_ha,tile_dec)
        ra=lst-ha

        x,y = fiberassign_radec2xy_cs5(ra,dec,tile_ra,tile_dec,tile_mjd,tile_ha,tile_fieldrot)

        print("Testing fiberassign_cs5_xy2radec")
        ra2,dec2 = fiberassign_cs5_xy2radec(x,y,tile_ra,tile_dec,tile_mjd,tile_ha,tile_fieldrot)

        drad=cosd(dec)*(ra2-ra)
        ddec=(dec2-dec)
        distdeg=np.sqrt(drad**2+ddec**2)
        rms=(np.sqrt(np.mean(distdeg**2)))
        print("rms={} arcsec".format(rms*3600.))
        assert(rms*3600.<0.1)


        print("Testing fiberassign_radec2xy_flat")
        x,y = fiberassign_radec2xy_flat(ra,dec,tile_ra,tile_dec,tile_mjd,tile_ha,tile_fieldrot)

        print("Testing fiberassign_flat_xy2radec")
        ra2,dec2 = fiberassign_flat_xy2radec(x,y,tile_ra,tile_dec,tile_mjd,tile_ha,tile_fieldrot)

        drad=cosd(dec)*(ra2-ra)
        ddec=(dec2-dec)
        distdeg=np.sqrt(drad**2+ddec**2)
        rms=(np.sqrt(np.mean(distdeg**2)))
        print("rms={} arcsec".format(rms*3600.))
        assert(rms*3600.<0.1)

if __name__ == '__main__':
    unittest.main()
