import unittest
from pkg_resources import resource_filename
from collections import Counter

import numpy as np
from desimeter.simplecorr import SimpleCorr
        
class TestCorr(unittest.TestCase):
    
     def test_simple_corr(self):
        print("Testing match with arbitrary translation and dilatation")
        nn=12
        x1=np.random.uniform(size=nn)-0.5
        y1=np.random.uniform(size=nn)-0.5
        
        # arb. dilatation
        xscale = 12.
        yscale = 13.
        x3 = xscale*x1
        y3 = yscale*y1
        # arb. translation
        x3 += 33.
        y3 += 11.
        
        corr31=SimpleCorr()
        corr31.fit(x3,y3,x1,y1)
        x3b,y3b=corr31.apply(x3,y3)
        
        dist=np.sqrt((x3b-x1)**2+(y3b-y1)**2)
        assert(np.all(dist<1e-6))

        print("Testing inverse")
        x1b,y1b=corr31.apply_inverse(x1,y1)

        dist=np.sqrt((x1b-x3)**2+(y1b-y3)**2)
        assert(np.all(dist<1e-6))

        
if __name__ == '__main__':
    unittest.main()
    
