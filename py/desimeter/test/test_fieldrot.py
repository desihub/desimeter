import unittest
import numpy as np
from desimeter.fieldmodel import fieldrot, dfieldrotdt_physical_model, dfieldrotdt_empirical_model

class TestFieldRot(unittest.TestCase):

    def test_fieldrot(self):

        # test scalars and compare with physical model
        mjd = 58000.
        hexrot_deg = 0.
        for dec in [-30,0,30,60] :
            for ha in [-40,0,40] :
                ra=0
                lst_deg=ha
                angle = fieldrot(ra=ra,dec=dec,mjd=mjd,lst_deg=lst_deg,hexrot_deg=hexrot_deg)
                rate1 = dfieldrotdt_physical_model(ra=ra,dec=dec,mjd=mjd,lst_deg=lst_deg)
                rate2 = dfieldrotdt_empirical_model(ra=ra,dec=dec,lst_deg=lst_deg)
                print("HA={} Dec={} ANGLE={:4.3f} deg RATE1={:3.2f} arcsec/min RATE2={:3.2f} arcsec/min  DIFF={:4.3f} arcsec/min".format(ha,dec,angle,rate1,rate2,rate1-rate2))

        # test 1d arrays
        ha=np.linspace(-60,60,3*4)
        dec=np.repeat(30.,ha.shape)
        rate=dfieldrotdt_empirical_model(ra=ra,dec=dec,lst_deg=lst_deg)
        # test 2d arrays
        ha=ha.reshape((3,4))
        dec=dec.reshape((3,4))
        rate=dfieldrotdt_empirical_model(ra=ra,dec=dec,lst_deg=lst_deg)


if __name__ == '__main__' :

    unittest.main()
