
"""
Utility class to fit a translation , dilatation , rotation in a cartesian plane 
"""

import json
import numpy as np

class SimpleCorr(object):

    def __init__(self) :
        self.dx = 0.
        self.dy = 0.
        self.sxx = 1.
        self.syy = 1.
        self.sxy = 0.
        self.rot_deg = 0.
        self.nmatch = 0
        self.rms = 0.
    
    def tojson(self):
        params = dict()
        params['name'] = "Simple Correction"
        params['version'] = '1'
        params['dx'] = self.dx
        params['dy'] = self.dy
        params['sxx'] = self.sxx
        params['syy'] = self.syy
        params['sxy'] = self.sxy
        params['rot_deg'] = self.rot_deg        
        params['nmatch'] = self.nmatch
        params['rms'] = self.rms      
        return json.dumps(params)

    @classmethod
    def fromjson(cls, jsonstring):
        tx = cls()
        params = json.loads(jsonstring)
        assert params['name'] == "Tangent Plane Adjustment"
        assert params['version'] == '1'
        tx.dx = params['dx']
        tx.dy = params['dy']
        tx.sxx = params['sxx']
        tx.syy = params['syy']
        tx.sxy = params['sxy']
        tx.rot_deg = params['rot_deg']
        tx.nmatch = params['nmatch']
        tx.rms = params['rms']        
        return tx

    def __str__(self):
        return "desimeter.simplecorr.SimpleCorr\n dx= {:7.6f}\n dy= {:7.6f}\n rot_deg= {:5.4f}\n sxx= {:5.4f}\n syy= {:5.4f}\n sxy= {:5.4f}\n nmatch= {:d}\n rms= {:7.6f} \n".format(self.dx,self.dy,self.rot_deg,self.sxx,self.syy,self.sxy,self.nmatch,self.rms)

    def fit(self, x1, y1, x2, y2) :
        """
        Adjust tranformation from x1,y1 to x2,y2
        
        Args:
           x1,y1,x2,y2 : 1D np.arrays of coordinates in tangent plane

        Returns:
           None

        """

        assert((x1.shape == y1.shape)&(x2.shape == y2.shape)&(x1.shape == x2.shape))
                
        # now fit simultaneously extra offset, rotation, scale
        self.nmatch=x1.size
        H=np.zeros((3,self.nmatch))
        H[0] = 1.
        H[1] = x1
        H[2] = y1
        A = H.dot(H.T)
        Ai = np.linalg.inv(A)
        ax = Ai.dot(np.sum(x2*H,axis=1))
        x2p = ax[0] + ax[1]*x1 + ax[2]*y1 # x2p = predicted x2 from x1 (=x1 after telescope pointing offset)
        ay = Ai.dot(np.sum(y2*H,axis=1))
        y2p = ay[0] + ay[1]*x1 + ay[2]*y1 # y2p = predicted y2 from y1

         # tangent plane coordinates are in radians
        self.rms = np.sqrt( np.mean( (x2-x2p)**2 + (y2-y2p)**2 ) )
        
        # pointing offset
        self.dx = ax[0]
        self.dy = ay[0]
        
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

        self.rot_deg = np.arctan2(sa,ca)*180/np.pi

        sxy = sa*ax[1]+ca*ay[1] - sxx_p_syy*ca*sa
        sxx =(ax[1]-sxy*sa)/ca
        syy = (ay[1]-sxy*ca)/sa

        self.sxx = sxx
        self.syy = syy
        self.sxy = sxy

    def apply(self,x,y) :
        scale_matrix = np.array([[self.sxx,self.sxy],[self.sxy,self.syy]])
        ca=np.cos(self.rot_deg/180*np.pi)
        sa=np.sin(self.rot_deg/180*np.pi)
        rot_matrix = np.array([[ca,-sa],[sa,ca]])
        xy=scale_matrix.dot(rot_matrix.dot(np.array([x,y])))
        return xy[0]+self.dx,xy[1]+self.dy
  
