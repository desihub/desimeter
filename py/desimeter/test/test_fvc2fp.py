"""
Test FVC <-> FP transformation code
"""

import os, glob
import unittest
from pkg_resources import resource_filename
import tempfile

import numpy as np
from astropy.table import Table

from desimeter.transform.fvc2fp.poly2d import FVCFP_Polynomial
from desimeter.transform.fvc2fp.zb import FVCFP_ZhaoBurge

class _TestFVC2FP(object):
    
    @classmethod
    def setUpClass(cls):
        cls.TransformClass = None
        cls.spotfile = resource_filename('desimeter', 'test/data/test-spots.csv')
        cls.tempdir = tempfile.mkdtemp()

    def setUp(self):
        self.spotfile = resource_filename('desimeter', 'test/data/test-spots.csv')

    @classmethod
    def tearDownClass(cls):
        jsonfile = os.path.join(cls.tempdir, 'blat.json')
        if os.path.exists(jsonfile):
            os.remove(jsonfile)
        os.rmdir(cls.tempdir)

    def test_fit(self):
        tx = self.TransformClass()
        spots = Table.read(self.spotfile)
        spots_orig = spots.copy()
        
        #- basic fit isn't supposed to modify input
        tx.fit(spots)
        self.assertTrue(spots.colnames == spots_orig.colnames)
        for colname in spots.colnames:
            yep = np.all(spots[colname] == spots_orig[colname])
            self.assertTrue(yep, 'Column {} changed values'.format(colname))

        if 'X_FP_METRO' in spots.colnames:
            spots.remove_column('X_FP_METRO')

        if 'Y_FP_METRO' in spots.colnames:
            spots.remove_column('Y_FP_METRO')

        #- Updating the spots should update X_FP,Y_FP but not XPIX,YPIX
        tx.fit(spots, update_spots=True)
        self.assertTrue(np.all(spots['XPIX'] == spots_orig['XPIX']))
        self.assertTrue(np.all(spots['YPIX'] == spots_orig['YPIX']))
        self.assertFalse(np.all(spots['X_FP'] == spots_orig['X_FP']))
        self.assertFalse(np.all(spots['Y_FP'] == spots_orig['Y_FP']))

        self.assertIn('X_FP_METRO', spots.colnames)
        self.assertIn('Y_FP_METRO', spots.colnames)
        ii = spots['PINHOLE_ID'] > 0
        dx = (spots['X_FP'] - spots['X_FP_METRO'])[ii]
        dy = (spots['Y_FP'] - spots['Y_FP_METRO'])[ii]
        dr = np.sqrt(dx**2 + dy**2)
        self.assertLess(np.median(dr), 0.050)

    def test_transform(self):
        tx = self.TransformClass()
        spots = Table.read(self.spotfile)
        tx.fit(spots, update_spots=True)

        #- make a grid of different spots not used in the fit
        nspots = 20
        xpix = np.random.uniform(500, 5500, nspots)
        ypix = np.random.uniform(500, 5500, nspots)
        ii = (xpix-3000)**2 + (ypix-3000)**2 < 2500**2
        xpix = xpix[ii]
        ypix = ypix[ii]
        xfp, yfp = tx.fvc2fp(xpix, ypix)
        xpix_new, ypix_new = tx.fp2fvc(xfp, yfp)

        #- test roundtrip with close cuts in FVC pixel space
        dr = np.sqrt((xpix - xpix_new)**2 + (ypix - ypix_new)**2)
        self.assertLess(np.median(dr), 0.10)
        self.assertLess(np.max(dr), 0.20)

        #- Transform FVC -> FP should be close to the FP metrology
        ii = (spots['X_FP_METRO']>0) & (spots['Y_FP_METRO']>0)
        dx = spots['X_FP'][ii] - spots['X_FP_METRO'][ii]
        dy = spots['Y_FP'][ii] - spots['Y_FP_METRO'][ii]
        dr = np.sqrt(dx**2 + dy**2)
        self.assertLess(np.median(dr), 0.05)  #- median < 50 microns
        self.assertLess(np.max(dr), 0.10)     #- worst outlier < 100 microns


    def test_fit_transforms(self):
        '''
        Test fits stability with FVC offsets, rotations, and scales
        '''
        spots = Table.read(self.spotfile)
        tx = self.TransformClass()

        #- Original
        tx.fit(spots)
        xfp0, yfp0 = tx.fvc2fp(spots['XPIX'], spots['YPIX'])
        
        #- Translate by (1,2) pixels
        spots['XPIX'] += 1.0
        spots['YPIX'] += 2.0
        tx.fit(spots)
        xfp1, yfp1 = tx.fvc2fp(spots['XPIX'], spots['YPIX'])

        #- Scale by 1% from center of image (3000,3000)
        spots['XPIX'] = (spots['XPIX']-3000)*1.01 + 3000
        spots['YPIX'] = (spots['YPIX']-3000)*1.01 + 3000
        tx.fit(spots)
        xfp2, yfp2 = tx.fvc2fp(spots['XPIX'], spots['YPIX'])

        #- Rotate by 1 degree about (3000,3000)
        theta = np.radians(1)
        dx = spots['XPIX'] - 3000
        dy = spots['YPIX'] - 3000
        x = np.cos(theta)*dx - np.sin(theta)*dy + 3000
        y = np.sin(theta)*dx + np.cos(theta)*dy + 3000
        spots['XPIX'] = x
        spots['YPIX'] = y
        tx.fit(spots)
        xfp3, yfp3 = tx.fvc2fp(spots['XPIX'], spots['YPIX'])

        #- 0.01 micron roundtrip RMS consistency; 0.1 micron max error
        self.assertLess(np.std(xfp1-xfp0), 1e-5)
        self.assertLess(np.std(xfp2-xfp0), 1e-5)
        self.assertLess(np.std(xfp3-xfp0), 1e-5)
        self.assertLess(np.max(np.abs(xfp1-xfp0)), 1e-4)
        self.assertLess(np.max(np.abs(xfp2-xfp0)), 1e-4)
        self.assertLess(np.max(np.abs(xfp3-xfp0)), 1e-4)

    def test_json(self):
        spots = Table.read(self.spotfile)
        t1 = self.TransformClass()
        t1.fit(spots)
        filename = os.path.join(self.tempdir, 'blat.json')
        t1.write_jsonfile(filename)

        t2 = self.TransformClass.read_jsonfile(filename)

        #- Compare transforms
        nspots = 20
        xpix = np.random.uniform(500, 5500, nspots)
        ypix = np.random.uniform(500, 5500, nspots)
        x1, y1 = t1.fvc2fp(xpix, ypix)
        x2, y2 = t2.fvc2fp(xpix, ypix)

        self.assertTrue(np.allclose(x1, x2))
        self.assertTrue(np.allclose(y1, y2))

class TestPoly2d(_TestFVC2FP, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.spotfile = resource_filename('desimeter', 'test/data/test-spots.csv')
        cls.tempdir = tempfile.mkdtemp()
        cls.TransformClass = FVCFP_Polynomial

class TestZhaoBurge(_TestFVC2FP, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.spotfile = resource_filename('desimeter', 'test/data/test-spots.csv')
        cls.tempdir = tempfile.mkdtemp()
        cls.TransformClass = FVCFP_ZhaoBurge

    def test_reduce_expand(self):
        from desimeter.transform.fvc2fp.zb import _reduce_xyfvc, _expand_xyfvc
        x1, y1 = np.random.uniform(2000,4000, size=(2,100))
        rx, ry = _reduce_xyfvc(x1, y1)
        x2, y2 = _expand_xyfvc(rx, ry)
        self.assertTrue(np.all(np.abs(rx)<1))
        self.assertTrue(np.all(np.abs(ry)<1))
        self.assertTrue(np.allclose(x1, x2))
        self.assertTrue(np.allclose(y1, y2))

        from desimeter.transform.fvc2fp.zb import _reduce_xyfp, _expand_xyfp
        x1, y1 = np.random.uniform(-300, 300, size=(2,100))
        rx, ry = _reduce_xyfp(x1, y1)
        x2, y2 = _expand_xyfp(rx, ry)
        self.assertTrue(np.all(np.abs(rx)<1))
        self.assertTrue(np.all(np.abs(ry)<1))
        self.assertTrue(np.allclose(x1, x2))
        self.assertTrue(np.allclose(y1, y2))

if __name__ == '__main__':
    unittest.main()


