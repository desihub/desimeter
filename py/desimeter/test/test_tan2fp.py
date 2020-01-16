import unittest
from pkg_resources import resource_filename

import numpy as np

from desimeter.transform.tan2fp.echo22 import (
    radius2theta, theta2radius, tan2fp, fp2tan)

class TestTan2FP(unittest.TestCase):

    def test_signs(self):
        #- phi=0 is aligned with +RA = +xtan = -xfp
        #- phi=90 is aligned with +dec = +ytan = +yfp
        xfp, yfp = tan2fp(0.01, 0.0)
        self.assertLess(xfp, 0.0)
        self.assertAlmostEqual(yfp, 0.0)

        xfp, yfp = tan2fp(0.00, 0.01)
        self.assertGreater(yfp, 0.0)
        self.assertAlmostEqual(xfp, 0.0)

        xtan, ytan = fp2tan(100, 0)
        self.assertLess(xtan, 0.0)
        self.assertAlmostEqual(ytan, 0.0)

        xtan, ytan = fp2tan(0, 100)
        self.assertGreater(ytan, 0.0)
        self.assertAlmostEqual(xtan, 0.0)

    def test_tan2fp(self):
        a = np.arcsin(np.radians(1.63))
        xtan1, ytan1 = np.random.uniform(-a, a, size=(2, 1000))
        keep = (xtan1**2 + ytan1**2) < a**2
        xtan1 = xtan1[keep]
        ytan1 = ytan1[keep]

        xfp, yfp = tan2fp(xtan1, ytan1)
        rfp = np.sqrt(xfp**2 + yfp**2)
        self.assertLess(np.max(rfp), 415)
        
        xtan2, ytan2 = fp2tan(xfp, yfp)
        dx = xtan2 - xtan1
        dy = ytan2 - ytan1
        dr = np.sqrt(dx**2 + dy**2)
        self.assertLess(np.max(dr), a*2e-6)
        self.assertLess(np.sqrt(np.mean(dr**2)), a*1e-6)
    
    def test_fp2tan(self):
        rmax = 414.0
        a = np.arcsin(np.radians(1.63))
        
        xfp1, yfp1 = np.random.uniform(-rmax, rmax, size=(2, 1000))
        keep = (xfp1**2 + yfp1**2) < rmax**2
        xfp1, yfp1 = xfp1[keep], yfp1[keep]
        
        xtan, ytan = fp2tan(xfp1, yfp1)
        rtan = np.sqrt(xtan**2 + ytan**2)
        self.assertLess(np.max(rtan), a)
        
        xfp2, yfp2 = tan2fp(xtan, ytan)
        dx, dy = xfp2-xfp1, yfp2-yfp1
        dr = np.sqrt(dx**2 + dy**2)

        self.assertLess(np.max(dr), 1e-3)   #- 1 um
        self.assertLess(np.sqrt(np.mean(dr**2)), 0.5e-3)   #- 0.5 um
        
    
    def test_radius2theta(self):
        r1 = np.arange(0,415.1, 1.0)
        theta = radius2theta(r1)
        self.assertLess(np.max(theta), 1.66)
        self.assertGreaterEqual(np.min(theta), 0.0)

        #- Test round trip
        r2 = theta2radius(theta)
        dr = 1000*(r2-r1)  #- difference in microns
        self.assertLess(np.sqrt(np.mean(dr**2)), 0.5)
        self.assertLess(np.max(np.abs(dr)), 1.0)

    def test_theta2radius(self):
        t1 = np.linspace(0, 1.63)
        radius = theta2radius(t1)
        self.assertLess(np.max(radius), 415.0)
        self.assertGreaterEqual(np.min(radius), 0.0)

        #- Test round trip
        t2 = radius2theta(radius)
        dt = 3600*(t2-t1)  #- different in arcsec
        self.assertLess(np.sqrt(np.mean(dt**2)), 0.05)
        self.assertLess(np.max(np.abs(dt)), 0.1)

if __name__ == '__main__':
    unittest.main()
