"""
Utility functions to fit and apply coordinates transformation from desimeter to PlateMaker focal plane coordinates XY.
They are supposed to be the same in theory but have in practice some systematic differences due to different interpolation functions
and different coordinates for fiducials, GFAs.
"""

import json
import numpy as np
from desimeter.transform.zhaoburge import getZhaoBurgeXY,fitZhaoBurge

#-------------------------------------------------------------------------

class DM2PM(object):

    def __init__(self):
        """
        init
        """
        self.zbpolids = np.array([0,1,2,3,4,5,6,9,20,27,28,29,30],dtype=int)
        self.dm2pm_zbcoeffs = np.zeros(self.zbpolids.shape,dtype=float)
        self.pm2dm_zbcoeffs = np.zeros(self.zbpolids.shape,dtype=float)

    @classmethod
    def read_jsonfile(cls, filename):
        with open(filename) as fx:
            s = fx.read()
        return cls.fromjson(s)

    @classmethod
    def read(cls,filename):
        if filename.find(".json")>=0 :
            return cls.read_jsonfile(filename)
        else :
            raise RuntimeError("don't know how to read {}".format(filename))

    def write_jsonfile(self, filename):
        with open(filename, 'w') as fx:
            fx.write(self.tojson())

    def write(self,filename):
        if filename.find(".json")>=0 :
            return self.write_jsonfile(filename)
        else :
            raise RuntimeError("don't know how to write {}".format(filename))

    #- Utility transforms to/from reduced [-1,1] coordinates
    def _reduce_xy(self, x, y):
        """
        Rescale FP xy coordinates [-420,420] -> [-1,1] and flip x axis
        """
        a = 420.0
        return x/a, y/a

    def _expand_xy(self, x, y):
        """
        Undo _redux_xyfp() transform
        """
        a = 420.0
        return x*a, y*a

    def tojson(self):
        params = dict()
        params['method'] = 'Zhao-Burge'
        params['version'] = '1'
        params['zbpolids'] = [int(polid) for polid in self.zbpolids]
        params['dm2pm_zbcoeffs'] = list(self.dm2pm_zbcoeffs)
        params['pm2dm_zbcoeffs'] = list(self.pm2dm_zbcoeffs)
        return json.dumps(params)

    def __str__(self) :
        return self.tojson()

    @classmethod
    def fromjson(cls, jsonstring):
        tx = cls()
        params = json.loads(jsonstring)
        assert params['method'] == 'Zhao-Burge'
        if params['version'] == '1' :
            tx.zbpolids = np.asarray(params['zbpolids'])
            tx.dm2pm_zbcoeffs = np.asarray(params['dm2pm_zbcoeffs']).astype(float)
            tx.pm2dm_zbcoeffs = np.asarray(params['pm2dm_zbcoeffs']).astype(float)
        else :
            raise RuntimeError("don't know version {}".format(params['version']))
        return tx

    def fit(self,xdm,ydm,xpm,ypm):
        """ Fit transformation.
        """
        rxdm, rydm = self._reduce_xy(xdm,ydm)
        rxpm, rypm = self._reduce_xy(xpm,ypm)
        _ , self.dm2pm_zbcoeffs = fitZhaoBurge(rxdm, rydm, rxpm, rypm, polids=self.zbpolids)
        _ , self.pm2dm_zbcoeffs = fitZhaoBurge(rxpm, rypm, rxdm, rydm, polids=self.zbpolids)

    def dm2pm(self, x, y):
        """
        Converts desimeter to platemaker (both CS5 coordinates, in mm)
        """
        rx,ry   = self._reduce_xy(x,y)
        dx,dy = getZhaoBurgeXY(self.zbpolids,self.dm2pm_zbcoeffs,rx,ry)
        x2,y2   = self._expand_xy(rx+dx,ry+dy)
        return x2,y2

    def pm2dm(self, x, y):
        """
        Converts platemaker to desimeter (both CS5 coordinates, in mm)
        """
        rx,ry   = self._reduce_xy(x,y)
        dx,dy = getZhaoBurgeXY(self.zbpolids,self.pm2dm_zbcoeffs,rx,ry)
        x2,y2   = self._expand_xy(rx+dx,ry+dy)
        return x2,y2
