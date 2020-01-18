import unittest
from pkg_resources import resource_filename

import numpy as np

from desimeter.transform.gfa2fp import (
    apply_scale_rotation_offset, undo_scale_rotation_offset)
from desimeter.transform.gfa2fp import gfa2fp, fp2gfa
from desimeter import io

class TestTan2FP(unittest.TestCase):

    def test_gfa2fp(self):
        metrology = io.load_metrology()
        ii = (metrology['DEVICE_TYPE'] == 'GFA')
        metrology = metrology[ii]
        
        nx, ny = 2048, 1032
        xgfa_ref = np.array([0, nx-1, nx-1, 0])
        ygfa_ref = np.array([0, 0, ny-1, ny-1])
        
        for p in range(10):
            ii = (metrology['PETAL_LOC'] == p)
            if np.count_nonzero(ii) == 0:
                print('ERROR: missing GFA metrology for PETAL_LOC {}'.format(p))
                continue
            
            xfp = np.asarray(metrology['X_FP'][ii])
            yfp = np.asarray(metrology['Y_FP'][ii])
            
            xgfa, ygfa = fp2gfa(p, xfp, yfp)
            
            if p == 0:
                #--- DEBUG ---
                import IPython
                IPython.embed()
                #--- DEBUG ---

    def test_scale_rotation_offset(self):
        x1 = np.array([1.0, 2.0, 3.0, 4.0])
        y1 = np.array([2.0, 4.0, 3.0, 1.0])
        xoff = -1.0
        yoff = 2
        scale = 1.3
        rotation = 20

        #- Test round trip accuracy
        xx, yy = apply_scale_rotation_offset(x1, y1, scale, rotation, xoff, yoff)
        x2, y2 = undo_scale_rotation_offset(xx, yy, scale, rotation, xoff, yoff)
        
        self.assertTrue(np.allclose(x1, x2))
        self.assertTrue(np.allclose(y1, y2))
        
        #- Test scalar input
        xx, yy = apply_scale_rotation_offset(1.0, 2.0, scale, rotation, xoff, yoff)
        self.assertTrue(np.isscalar(xx))
        self.assertTrue(np.isscalar(yy))
        xx, yy = undo_scale_rotation_offset(1.0, 2.0, scale, rotation, xoff, yoff)
        self.assertTrue(np.isscalar(xx))
        self.assertTrue(np.isscalar(yy))

        #- Test non-array list
        xx, yy = apply_scale_rotation_offset([1.0, 2.0], [2, 3], scale, rotation, xoff, yoff)
        self.assertEqual(len(xx), 2)
        self.assertEqual(len(yy), 2)

        xx, yy = undo_scale_rotation_offset([1.0, 2.0], [2, 3], scale, rotation, xoff, yoff)
        self.assertEqual(len(xx), 2)
        self.assertEqual(len(yy), 2)

if __name__ == '__main__':
    unittest.main()
