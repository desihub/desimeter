import unittest
from pkg_resources import resource_filename
from collections import Counter

import numpy as np
from desimeter.match import match_same_system,match_arbitrary_translation_dilatation

class TestMatch(unittest.TestCase):

    def x1y1x2y2(self) :
        nn=12
        x1=np.random.uniform(size=nn)-0.5
        y1=np.random.uniform(size=nn)-0.5
        ii=np.arange(nn,dtype=int)
        np.random.shuffle(ii)
        x2=x1[ii].copy()
        y2=y1[ii].copy()
        return x1,y1,x2,y2
    
    def test_same_coordinate_system(self):
        print("Testing match in same coordinate system")
        x1,y1,x2,y2 = self.x1y1x2y2()
        indices_2,distances = match_same_system(x1,y1,x2,y2)
        ii2=indices_2[np.arange(x1.size)]
        dist=np.sqrt( (x1-x2[ii2])**2+(y1-y2[ii2])**2 )
        assert(np.all(dist==0.))

    def test_arbitrary_translation_dilatation(self):
        print("Testing match with arbitrary translation and dilatation")
        x1,y1,x2,y2 = self.x1y1x2y2()
        # arb. dilatation
        scale = 12.
        x3 = scale*x2
        y3 = scale*y2
        # arb. translation
        x3 += 33.
        y3 += 11.
        indices_2,distances = match_same_system(x1,y1,x2,y2)
        ii2=indices_2[np.arange(x1.size)]
        dist=np.sqrt( (x1-x2[ii2])**2+(y1-y2[ii2])**2 )
        assert(np.all(dist==0.))
        
        
if __name__ == '__main__':
    unittest.main()
    
