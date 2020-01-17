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

        for cha in [60] :
            for cdec in [80] :
                
                ha1,dec1 = radec2tan.xy2hadec(x1,y1,cha,cdec)

                # apply rotation
                M = radec2tan.compute_polar_misalignment_rotation_matrix(me_arcsec=radec2tan.ME_ARCSEC,ma_arcsec=radec2tan.MA_ARCSEC)
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
                print("at HA= {} deg and Dec = {} deg, the mean field rotation angle= {:4.3f} deg ".format(cha,cdec,angle))
        
    def test_precession(self):
        ra = 12.
        dec = 24.
        mjd = 58864.3
        ra2,dec2 = radec2tan.apply_precession_from_icrs(ra, dec, mjd, use_astropy = False)
        ra2b,dec2b = radec2tan.apply_precession_from_icrs(ra, dec, mjd, use_astropy = True)
        print("precession this code dRA= {:4.1f} arcsec , dDec= {:4.1f} arcsec".format((ra2-ra)*3600.,(dec2-dec)*3600.))
        print("precession astropy   dRA= {:4.1f} arcsec , dDec= {:4.1f} arcsec".format((ra2b-ra)*3600.,(dec2b-dec)*3600.))
        assert(np.abs(ra2-ra2b)<1.) # we care only about the relative variation of this in the focal plane
        assert(np.abs(dec2-dec2b)<1.) # we care only about the relative variation of this in the focal plane

        # test 1D
        ra = ra*np.ones(4)
        dec = dec*np.ones(4)
        ra2,dec2 = radec2tan.apply_precession_from_icrs(ra, dec, mjd, use_astropy = False)
        ra2b,dec2b = radec2tan.apply_precession_from_icrs(ra, dec, mjd, use_astropy = True)
        print("precession this code dRA= {} arcsec , dDec= {} arcsec".format((ra2-ra)*3600.,(dec2-dec)*3600.))
        print("precession astropy   dRA= {} arcsec , dDec= {} arcsec".format((ra2b-ra)*3600.,(dec2b-dec)*3600.))
        
    def test_aberration(self):
        ra = 12.
        dec = 24.
        mjd = 58864.3
        ra2,dec2 = radec2tan.apply_aberration(ra, dec, mjd, use_astropy = False)
        ra2b,dec2b = radec2tan.apply_aberration(ra, dec, mjd, use_astropy = True)
        print("aberration this code dRA= {:4.1f} arcsec , dDec= {:4.1f} arcsec".format((ra2-ra)*3600.,(dec2-dec)*3600.))
        print("aberration astropy   dRA= {:4.1f} arcsec , dDec= {:4.1f} arcsec".format((ra2b-ra)*3600.,(dec2b-dec)*3600.))
        assert(np.abs(ra2-ra2b)<1.) # we care only about the relative variation of this in the focal plane
        assert(np.abs(dec2-dec2b)<1.) # we care only about the relative variation of this in the focal plane

        # test 1D
        ra = ra*np.ones(4)
        dec = dec*np.ones(4)
        ra2,dec2 = radec2tan.apply_aberration(ra, dec, mjd, use_astropy = False)
        ra2b,dec2b = radec2tan.apply_aberration(ra, dec, mjd, use_astropy = True)
        print("aberration this code dRA= {} arcsec , dDec= {} arcsec".format((ra2-ra)*3600.,(dec2-dec)*3600.))
        print("aberration astropy   dRA= {} arcsec , dDec= {} arcsec".format((ra2b-ra)*3600.,(dec2b-dec)*3600.))
        
    def test_tan2radec(self):
        cra=12
        cdec=24
        lst = 13.
        ra = cra + 1.
        dec = cdec + 1.
        mjd = 58864.3
        x,y = radec2tan.radec2tan(ra,dec,tel_ra=cra,tel_dec=cdec,mjd=mjd,lst_deg=lst,hexrot_deg=0.1)
        ra2,dec2 = radec2tan.tan2radec(x,y,tel_ra=cra,tel_dec=cdec,mjd=mjd,lst_deg=lst,hexrot_deg=0.1)
        dra=ra2-ra
        ddec=dec2-dec
        print("dra={} arcsec , ddec={} arcsec".format(dra*3600,ddec*3600))
        assert(np.abs(dra*3600.)<0.001) # 1 milli-arcsec 
        assert(np.abs(ddec*3600.)<0.001) # 1 milli-arcsec 
        
    def test_refraction(self):
        alt = 30. # deg
        alt2 = radec2tan.apply_refraction(alt)
        alt3 = radec2tan.undo_refraction(alt2)
        print("At altitude={} deg, refraction={:4.1f} arcsec, and zero={:4.3f} arcsec".format(alt,(alt2-alt)*3600.,(alt3-alt)*3600.))
        assert(np.abs((alt3-alt)*3600.)<0.001) # 1 milli-arcsec 
        
if __name__ == '__main__' :
    
    unittest.main()
    
