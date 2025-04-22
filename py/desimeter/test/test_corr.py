import unittest
from collections import Counter

import numpy as np
from desimeter.simplecorr import SimpleCorr
        
class TestCorr(unittest.TestCase):
    
     def test_simple_corr(self):
        print("Testing simple corr")
        nn=12
        x1=np.random.uniform(size=nn)-0.5
        y1=np.random.uniform(size=nn)-0.5
        
        # arb. angle
        a=10./180.*np.pi
        ca=np.cos(a)
        sa=np.sin(a)

        # arb. dilatation
        xscale = 12.
        yscale = 13.
        x3 =  xscale*ca*x1 + yscale*sa*y1
        y3 = -xscale*sa*x1 + yscale*ca*y1
        
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

     def test_fit_rotoff(self):
        print("Testing fit_rotoff")
        nn=12
        x1=np.random.uniform(size=nn)-0.5
        y1=np.random.uniform(size=nn)-0.5

        # arb. angle
        a=22./180.*np.pi
        ca=np.cos(a)
        sa=np.sin(a)
        
        # arb. translation
        x3 = 33.  + ca*x1 + sa*y1
        y3 = 11. - sa*x1 + ca*y1
        
        corr31=SimpleCorr()
        corr31.fit_rotoff(x3,y3,x1,y1)
        x3b,y3b=corr31.apply(x3,y3)
        dist=np.sqrt((x3b-x1)**2+(y3b-y1)**2)
        assert(np.all(dist<1e-8))

        corr13=SimpleCorr()
        corr13.fit_rotoff(x1,y1,x3,y3)
        x1b,y1b=corr13.apply(x1,y1)
        dist=np.sqrt((x3-x1b)**2+(y3-y1b)**2)
        assert(np.all(dist<1e-8))

        
if __name__ == '__main__':
    unittest.main()
    
