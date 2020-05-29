import unittest
from pkg_resources import resource_filename
from collections import Counter

import numpy as np
from desimeter.circles import fit_circle,robust_fit_circle

class TestCircles(unittest.TestCase):

     def test_fit_circle(self):
        print("Testing fit circle")
        xc=12
        yc=24
        r=3.

        nn=12
        a=np.linspace(0,2*np.pi,nn)
        x=xc+r*np.cos(a)
        y=yc+r*np.sin(a)

        xfit,yfit,rfit = fit_circle(x,y)
        print(xfit,yfit,rfit)
        assert(np.abs(xfit-xc)<1e-6)
        assert(np.abs(yfit-yc)<1e-6)
        assert(np.abs(rfit-r)<1e-6)

        xfit,yfit,rfit,ok = robust_fit_circle(x,y)

        nbad=np.sum(ok==False)
        print(xfit,yfit,rfit,nbad)
        assert(np.abs(xfit-xc)<1e-6)
        assert(np.abs(yfit-yc)<1e-6)
        assert(np.abs(rfit-r)<1e-6)
        assert(nbad==0)

        x[1] += 2.

        xfit,yfit,rfit,ok = robust_fit_circle(x,y)
        nbad=np.sum(ok==False)
        print(xfit,yfit,rfit,nbad)
        assert(np.abs(xfit-xc)<1e-6)
        assert(np.abs(yfit-yc)<1e-6)
        assert(np.abs(rfit-r)<1e-6)
        assert(nbad==1)



if __name__ == '__main__':
    unittest.main()
