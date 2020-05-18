"""
Test RA Dec <-> Tangent plane transformation code.
"""


import unittest

import numpy as np
from desimeter.fieldmodel import FieldModel,TanCorr
from desimeter.time import mjd2lst
from desimeter.transform.radec2tan import xy2hadec,hadec2xy
from desimeter.trig import cosd

class TestFieldModel(unittest.TestCase):


    def test_tancorr(self):
        fm = FieldModel()
        fm.sxx = 1.1 ; fm.syy = 1.2 ; fm.sxy = -0.05
        fm.fieldrot_deg = 10.
        x=np.linspace(-0.1,0.1,4)
        y=np.linspace(-0.2,0.2,4)
        x2,y2 = fm.tancorr_inst2sky(x,y)
        x3,y3 = fm.tancorr_sky2inst(x2,y2)
        mdiff=np.max((x3-x)**2+(y3-y)**2)
        assert(mdiff<1.e-12)

    def test_tancorr_fit(self):

        tcorr = TanCorr()

        tel_ha1 = 12.
        tel_dec1 = 80.

        npts = 24
        x1 = 0.1*np.random.uniform(size=npts)
        y1 = 0.1*np.random.uniform(size=npts)
        ha,dec = xy2hadec(x1,y1,tel_ha1,tel_dec1)

        true_dha_arcsec=3.
        true_ddec_arcsec=1.

        tel_ha2  = tel_ha1 + true_dha_arcsec/3600.
        tel_dec2 = tel_dec1 + true_ddec_arcsec/3600.
        x2,y2 = hadec2xy(ha,dec,tel_ha2,tel_dec2)

        tcorr.fit(x1, y1, x2, y2)
        meas_dha_arcsec  = tcorr.dha/cosd(tel_dec1)*3600.
        meas_ddec_arcsec = tcorr.ddec*3600.

        assert(np.abs(true_dha_arcsec-meas_dha_arcsec)<0.01)
        assert(np.abs(true_ddec_arcsec-meas_ddec_arcsec)<0.01)







    def test_fm(self):
        fm = FieldModel()
        fm.sxx = 1.1 ; fm.syy = 1.2 ; fm.sxy = -0.05
        fm.ra = 12.
        fm.dec = 24.
        fm.fieldrot_deg = 10.
        fm.hexrot_deg = 0.
        fm.mjd = 58800.
        fm.lst = mjd2lst(fm.mjd)
        fm.adc1 = 0.
        fm.adc2 = 0.

        x=np.linspace(0,200.,4) # mm
        y=np.linspace(0,100.,4) # mm
        ra,dec = fm.fp2radec(x,y)
        x2,y2 = fm.radec2fp(ra,dec)
        mdiff=np.max((x2-x)**2+(y2-y)**2)
        assert(mdiff<1e-4) # 0.1 um

    def test_fm_io(self):
        fm = FieldModel()
        fm.sxx = 1.1 ; fm.syy = 1.2 ; fm.sxy = -0.05
        fm.ra = 12.
        fm.dec = 24.
        fm.fieldrot_deg = 10.
        fm.hexrot_deg = 0.
        fm.mjd = 58800.
        fm.lst = mjd2lst(fm.mjd)
        fm.adc1 = 0.
        fm.adc2 = 0.

        string = fm.tojson()
        print(string)

        fm2 = FieldModel.fromjson(string)
        print(fm2.sxx)

        assert(np.abs(fm.sxx-fm2.sxx)<1e-8)






if __name__ == '__main__' :

    unittest.main()
