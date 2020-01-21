
"""
Utility functions to fit guide stars
"""

import json
import numpy as np

from desimeter.transform.radec2tan import xy2hadec,hadec2xy

class GuideStarsCorr(object):

    def __init__(self) :
        self.dha  = 0.
        self.ddec = 0.
        self.sxx = 1.
        self.syy = 1.
        self.sxy = 0.
        self.fieldrot_deg = 0.
        self.expid  = 0
        self.nstars = 0
        self.rms_arcsec = 0.
    
    def tojson(self):
        params = dict()
        params['name'] = "Tangent Plane Adjustment"
        params['version'] = '1'
        params['dha'] = self.dha
        params['ddec'] = self.ddec
        params['sxx'] = self.sxx
        params['syy'] = self.syy
        params['sxy'] = self.sxy
        params['fieldrot_deg'] = self.fieldrot_deg        
        params['expid']=self.expid
        params['nstars'] = self.nstars
        params['rms_arcsec'] = self.rms_arcsec        
        return json.dumps(params)

    @classmethod
    def fromjson(cls, jsonstring):
        tx = cls()
        params = json.loads(jsonstring)
        assert params['name'] == "Tangent Plane Adjustment"
        assert params['version'] == '1'
        tx.dha = params['dha']
        tx.ddec = params['ddec']
        tx.sxx = params['sxx']
        tx.syy = params['syy']
        tx.sxy = params['sxy']
        tx.fieldrot_deg = params['fieldrot_deg']
        tx.expid = params['expid']
        tx.nstars = params['nstars']
        tx.rms_arcsec = params['rms_arcsec']        
        return tx

    def fit(self, x1, y1, x2, y2) :
        """
        Adjust tranformation from focal plane x1,y1 to x2,y2
        
        Args:
           x1,y1,x2,y2 : 1D np.arrays of coordinates in tangent plane

        Returns:
           None

        """

        assert((x1.shape == y1.shape)&(x2.shape == y2.shape)&(x1.shape == x2.shape))
        
        # first ajust an offset using spherical coordinates
        # assume fiducial pointing of telescope to convert
        # tangent plane coords to angles
        self.dha=0.
        self.ddec=0.
        x1t=x1+0.
        y1t=y1+0.
        ha,dec = xy2hadec(x1,y1,0,0)
        for i in range(4) :
            x1t,y1t  = hadec2xy(ha,dec,self.dha,self.ddec)
            dx = np.mean(x2-x1t)
            dy = np.mean(y2-y1t)
            #print(i,dx,dy)
            self.dha -= dx/np.cos(self.ddec)*180./np.pi
            self.ddec -= dy*180./np.pi
        x1t,y1t  = hadec2xy(ha,dec,self.dha,self.ddec)

        # now fit simultaneously extra offset, rotation, scale
        self.nstars=x1t.size
        H=np.zeros((3,self.nstars))
        H[0] = 1.
        H[1] = x1t
        H[2] = y1t
        A = H.dot(H.T)
        Ai = np.linalg.inv(A)
        ax = Ai.dot(np.sum(x2*H,axis=1))
        x2p = ax[0] + ax[1]*x1t + ax[2]*y1t # x2p = predicted x2 from x1t (=x1 after telescope pointing offset)
        ay = Ai.dot(np.sum(y2*H,axis=1))
        y2p = ay[0] + ay[1]*x1t + ay[2]*y1t # y2p = predicted y2 from y1t

         # tangent plane coordinates are in radians
        self.rms_arcsec = np.sqrt( np.mean( (x2-x2p)**2 + (y2-y2p)**2 ) )*(180*3600)/np.pi
        
        # interpret this back into telescope pointing offset, field rotation, dilatation

        # pointing offset
        # increasing gaia stars x means telescope is more to the left so tel_ha should be decreased
        # increasing gaia stars y means telescope is more to the bottom so tel_dec should be decreased
        # tangent plane coordinates are in rad
        ddha  = -ax[0]*180./np.pi/np.cos(self.ddec/180*np.pi)
        dddec = -ay[0]*180./np.pi
        self.dha  += ddha
        self.ddec += dddec
        
        # dilatation and rotation
        # |ax1 ax2| |sxx sxy| |ca  -sa|
        # |ay1 ay2|=|syx syy|*|sa   ca|
        # ax1=sxx*ca+sxy*sa ; ax2=-sxx*sa+sxy*ca
        # ay1=syx*ca+syy*sa ; ay2=-syx*sa+syy*ca
        # ax1+ay2 = (sxx+syy)*ca
        # ay1-ax2 = (sxx+syy)*sa

        sxx_p_syy = np.sqrt( (ax[1]+ay[2])**2+(ay[1]-ax[2])**2 )
        sa=(ay[1]-ax[2])/sxx_p_syy
        ca=(ax[1]+ay[2])/sxx_p_syy

        self.fieldrot_deg = np.arctan2(sa,ca)*180/np.pi

        sxy = sa*ax[1]+ca*ay[1] - sxx_p_syy*ca*sa
        sxx =(ax[1]-sxy*sa)/ca
        syy = (ay[1]-sxy*ca)/sa

        self.sxx = sxx
        self.syy = syy
        self.sxy = sxy

    def apply(self,x,y) :
        scale_matrix = np.array([[self.sxx,self.sxy],[self.sxy,self.syy]])
        ca=np.cos(self.fieldrot_deg/180*np.pi)
        sa=np.sin(self.fieldrot_deg/180*np.pi)
        rot_matrix = np.array([[ca,-sa],[sa,ca]])
        ha,dec  = xy2hadec(x,y,0,0)
        x1t,y1t = hadec2xy(ha,dec,self.dha,self.ddec)
        xy=scale_matrix.dot(rot_matrix.dot(np.array([x1t,y1t])))
        return xy[0],xy[1]
        
        
        
