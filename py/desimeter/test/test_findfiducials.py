import unittest
from pkg_resources import resource_filename
from collections import Counter

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree as KDTree
from astropy.table import Table

from desimeter.findfiducials import findfiducials

class TestFindFiducials(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        '''
        Create test image based upon input spots
        '''
        cls.spotfile = resource_filename('desimeter', 'test/data/test-spots.csv')
        cls.spots = Table.read(cls.spotfile)
        cls.spots = cls.spots[:][cls.spots['PINHOLE_ID']>0] # keep only FIF or GIF
        cls.input_transform = resource_filename('desimeter',"data/canon-lens-fvc2fp.json")

    def setUp(self):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def test_findfiducials(self):
        spots = self.spots.copy()
        spots.remove_columns(['LOCATION', 'PINHOLE_ID', 'X_FP', 'Y_FP',
                              'X_FP_METRO', 'Y_FP_METRO'])
        spots = findfiducials(spots, input_transform=self.input_transform)
        #- All fiducial locations found
        self.assertTrue(np.all(np.in1d(self.spots['LOCATION'], spots['LOCATION'])))

        #- No extra locations found
        self.assertTrue(np.all(np.in1d(spots['LOCATION'], self.spots['LOCATION'])))

        #- FAILS: All individual pinholes found
        test_pinholes = list(zip(self.spots['LOCATION'], self.spots['PINHOLE_ID']))
        found_pinholes = list(zip(spots['LOCATION'], spots['PINHOLE_ID']))

        missing = set(test_pinholes) - set(found_pinholes)
        extra = set(found_pinholes) - set(test_pinholes)
        self.assertEqual(len(missing), 0, 'Missing pinholes: {}'.format(missing))

        #- No extra pinholes found
        self.assertEqual(len(extra), 0, 'Extra pinholes: {}'.format(extra))

    #- TODO: test that findfiducials is robust to translation, rotation, scale

if __name__ == '__main__':
    unittest.main()
