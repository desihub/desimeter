import unittest
from pkg_resources import resource_filename
from collections import Counter

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree as KDTree
from astropy.table import Table

from desimeter import detectspots

class TestDetectSpots(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        '''
        Create test image based upon input spots
        '''
        cls.spotfile = resource_filename('desimeter', 'test/data/test-spots.csv')
        cls.spots = spots = Table.read(cls.spotfile)
        np.random.seed(1)
        cls.img = np.random.normal(loc=1000.0, scale=35, size=(6000, 6000))

        spot = np.zeros((15,15))
        spot[7,7] = 1.0
        spot = detectspots.gaussian_convolve(spot)
        for x, y, counts in zip(spots['XPIX'], spots['YPIX'], spots['COUNTS']):
            dx = x % 1
            dy = y % 1
            x = int(x)
            y = int(y)
            cls.img[y-7:y+8, x-7:x+8] += counts*ndimage.shift(spot, (dy,dx))

    def setUp(self):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def test_detect(self):
        spots = detectspots.detectspots(self.img)
        self.assertTrue(len(spots) == len(self.spots))

        tree = KDTree(np.array((self.spots['XPIX'], self.spots['YPIX'])).T)
        distances, indices = tree.query(np.array((spots['XPIX'], spots['YPIX'])).T)

        #- all spots were matched
        self.assertEqual(len(set(indices)), len(indices))

        #- loose check on distances because original image wasn't constructed
        #- very exactly either; just check if we matched the right spot
        self.assertLess(np.max(distances), 0.02)

    def test_uint(self):
        #- Should also work with uint data (from raw FVC images)
        spots = detectspots.detectspots(self.img.astype(np.uint16))

if __name__ == '__main__':
    unittest.main()
