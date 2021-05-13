import unittest
from pkg_resources import resource_filename
from collections import Counter

import numpy as np

from desimeter.io import dm2pm_filename
from desimeter.transform.dm2pm import DM2PM

class TestDM2PM(unittest.TestCase):

     def test_dm2pm(self):
        print("Testing DM2PM")
        nn=200
        x1=200.*(np.random.uniform(size=nn)-0.5)
        y1=200.*(np.random.uniform(size=nn)-0.5)
        dm2pm = DM2PM.read(dm2pm_filename())
        x2,y2 = dm2pm.dm2pm(x1,y1)
        x3,y3 = dm2pm.pm2dm(x2,y2)
        dist=np.sqrt((x3-x1)**2+(y3-y1)**2)
        mdist=np.max(dist)
        print(mdist)
        assert(mdist<1e-3)

        # test fit
        x4=1.01*x1
        y4=1.01*y1
        dm2pm.fit(x1,y1,x4,y4)
        x5,y5=dm2pm.dm2pm(x1,y1)
        dist=np.sqrt((x5-x4)**2+(y5-y4)**2)
        mdist=np.max(dist)
        print(mdist)
        assert(mdist<1e-3)


if __name__ == '__main__':
    unittest.main()
