"""
Test RA Dec <-> Tangent plane transformation code
"""


import unittest

import numpy as np
import desimeter.transform.radec2tan as radec2tan


class TestRADEC2TAN(unittest.TestCase):
    
    
    def test_hadec2xy(self):

        cha=12
        cdec=24
        x1 = np.random.uniform(size=3)-0.5
        y1 = np.random.uniform(size=3)-0.5
        ha1,dec1 = radec2tan.xy2hadec(x1,y1,cha,cdec)
        x2,y2 = radec2tan.hadec2xy(ha1,dec1,cha,cdec)
        print("x1=",x1)
        print("y1=",y1)
        print("x2=",x2)
        print("y2=",y2)

        assert(np.all(np.abs(x1-x2)<1e-6))
        assert(np.all(np.abs(y1-y2)<1e-6))

        # check orientation
        x,y = radec2tan.hadec2xy(cha+1,cdec,cha,cdec)
        assert(x>0)
        x,y = radec2tan.hadec2xy(cha,cdec+1,cha,cdec)
        assert(y>0)

    def test_scalar_and_array(self):
        cra = 12.
        cdec = 24.
        lst = 13.
        ra = cra + 1.
        dec = cdec + 1.
        mjd = 58864.3
        
        x , y = radec2tan.radec2tan(ra,dec,tel_ra=cra,tel_dec=cdec,mjd=mjd,lst_deg=lst,hexrot_deg=0)
        # check we can convert to float
        x = float(x)
        y = float(y)

        npt = 12
        ra  = cra + np.linspace(1,1.,npt)
        dec = cdec + np.linspace(-1,1.,npt)
        x , y = radec2tan.radec2tan(ra,dec,tel_ra=cra,tel_dec=cdec,mjd=mjd,lst_deg=lst,hexrot_deg=0)
        # check dimension of output array
        assert(x.shape == ra.shape)
        assert(y.shape == ra.shape)


    def test_polar_misalignment(self):
        
        # circle of coordinates in FP
        phi = np.linspace(0,np.pi,4)
        theta = 1./180*np.pi
        x1 = np.sin(theta)*np.cos(phi)
        y1 = np.sin(theta)*np.sin(phi)

        #  pointings

        for cha in [-60,30,0,30,60] :
            for cdec in [80] :
                
                ha1,dec1 = radec2tan.xy2hadec(x1,y1,cha,cdec)

                # apply rotation
                M = radec2tan.compute_polar_misalignment_rotation_matrix()
                ha2,dec2 = radec2tan.getLONLAT(M.dot(radec2tan.getXYZ(ha1,dec1)))

                # back to FP coordinates
                cha2,cdec2 = radec2tan.getLONLAT(M.dot(radec2tan.getXYZ(cha,cdec)))
                x2,y2 = radec2tan.hadec2xy(ha2,dec2,cha2,cdec2)

                # measure field rotation angle
                x1 -= np.mean(x1)
                y1 -= np.mean(y1)
                x2 -= np.mean(x2)
                y2 -= np.mean(y2)
                angle =  np.mean(radec2tan.arcsind((x1*y2-x2*y1)/np.sqrt((x1**2+y1**2)*(x2**2+y2**2))))
                print("HA = {} , Dec = {}, mean rotation angle ={} deg ".format(cha,cdec,angle))
        
        
        
if __name__ == '__main__' :
    
    unittest.main()
    
