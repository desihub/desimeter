"""
Utility functions to fit and apply coordinates transformation from FVC to FP
"""
import json
from pkg_resources import resource_filename

import numpy as np
from scipy.interpolate import interp1d

from desimeter.transform.zhaoburge import getZhaoBurgeXY, transform, fitZhaoBurge
from desimeter.trig import average_angles_deg, put360

#- Utility transforms to/from reduced [-1,1] coordinates
def _reduce_xyfp(x, y):
    """
    Rescale FP xy coordinates
    """
    a =  420.0 #mm
    return x/a, y/a

def _expand_xyfp(x, y):
    """
    Undo _redux_xyfp() transform
    """
    a = 420.0 #mm
    return x*a, y*a

def _reduce_xytan(x, y):
    """
    Rescale tangent plane xy pix coords
    """
    a = 0.027
    return x/a, y/a

def _expand_xytan(x, y):
    """
    Undo _redux_xytan() transform
    """
    a = 0.027
    return x*a, y*a


#-------------------------------------------------------------------------

class TAN2FP_RayTraceFit(object) :
    def tojson(self):
        params = dict()
        params['method'] = 'Ray Trace fit with Zhao-Burge'
        params['version'] = '1'
        params['adc1'] = list(self.adc1)
        params['adc2'] = list(self.adc2)
        params['scale'] = list(self.scale)
        params['rotation'] = list(self.rotation)
        params['offset_x'] = list(self.offset_x)
        params['offset_y'] = list(self.offset_y)
        params['zbpolids'] = [int(polid) for polid in self.zbpolids]
        params['zbcoeffs'] = list(self.zbcoeffs.ravel())
        return json.dumps(params)

    @classmethod
    def fromjson(cls, jsonstring):
        tx = cls()
        params = json.loads(jsonstring)
        assert params['method'] == 'Ray Trace fit with Zhao-Burge'
        assert params['version'] == '1'
        tx.adc1 = np.asarray(params['adc1'])
        tx.adc2 = np.asarray(params['adc2'])
        tx.scale = np.asarray(params['scale'])
        tx.rotation = np.asarray(params['rotation'])
        tx.offset_x = np.asarray(params['offset_x'])
        tx.offset_y = np.asarray(params['offset_y'])
        tx.zbpolids = np.asarray(params['zbpolids'])
        tx.zbcoeffs = np.asarray(params['zbcoeffs']).reshape((tx.adc1.size,tx.zbpolids.size))
        return tx

    @classmethod
    def read_jsonfile(cls, filename):
        with open(filename) as fx:
            s = fx.read()
        return cls.fromjson(s)

    def write_jsonfile(self, filename):
        with open(filename, 'w') as fx:
            fx.write(self.tojson())

    def fit(self, table ) :
        """TODO: document"""
        #- identify ADC setups
        adc12=(np.array(table['ADC1'])+1000*np.array(table['ADC2']))
        uadc12 = np.unique(adc12)
        nconfig = len(uadc12)

        self.adc1     = np.zeros(nconfig,dtype=float)
        self.adc2     = np.zeros(nconfig,dtype=float)
        self.scale    = np.zeros(nconfig,dtype=float)
        self.rotation = np.zeros(nconfig,dtype=float)
        self.offset_x = np.zeros(nconfig,dtype=float)
        self.offset_y = np.zeros(nconfig,dtype=float)
        self.zbpolids = None
        self.zbcoeffs = list()

        for config, adc12v in enumerate(uadc12) :
            selection = (adc12==adc12v)
            self.adc1[config]=table['ADC1'][selection][0]
            self.adc2[config]=table['ADC2'][selection][0]
            print("Fitting ADC1={} ADC2={}".format(self.adc1[config],self.adc2[config]))

            #- Get reduced coordinates
            rxtan, rytan = _reduce_xytan(table['X_TAN'][selection], table['Y_TAN'][selection])
            rxfp, ryfp = _reduce_xyfp(table['X_FP'][selection], table['Y_FP'][selection])



            #################################################################
            ## CHOICE OF POLYNOMIALS IS HERE
            ##
            polids = np.array([0,1,2,3,4,5,6,9,14,15,20,27,28,29,30],dtype=int) #
            #################################################################
            #- Perform fit
            zbpolids, zbcoeffs =  fitZhaoBurge(rxtan, rytan, rxfp, ryfp, polids=polids)
            self.scale[config] = 1
            self.rotation[config] = 0.
            self.offset_x[config] = 0.
            self.offset_y[config] = 0.

            if self.zbpolids is None :
                self.zbpolids = zbpolids
            else :
                assert(np.all(self.zbpolids==zbpolids))
            self.zbcoeffs.append(zbcoeffs)

            #- Goodness of fit
            #xfp_fit, yfp_fit = self.tan2fp(table['X_TAN'][selection], table['Y_TAN'][selection])
            #dx = (table['X_FP'][selection] - xfp_fit)
            #dy = (table['Y_FP'][selection] - yfp_fit)
            #dr = np.sqrt(dx**2 + dy**2)
            #log.info('Mean, median, RMS distance = {:.1f}, {:.1f}, {:.1f} um'.format(
            #    1000*np.mean(dr), 1000*np.median(dr), 1000*np.sqrt(np.mean(dr**2))))

        self.zbcoeffs = np.vstack(self.zbcoeffs)

    def interpolate_coeffs(self,adc1, adc2):
        # interpolate transform parameters
        dadc_array = self.adc2-self.adc1
        dadc_array[dadc_array<0] += 360.
        sorted_indices=np.argsort(dadc_array)
        dadc_array = dadc_array[sorted_indices]

        dadc_arg   = put360(adc2 - adc1)
        # In principle, the two ADCs could be *set* to the same angle
        # but *measured* to be ADC2 - ADC1 = -1 degrees -> 359
        # degrees.  Convert that to zero, with an arbitrary limit of 1
        # degree error.
        if dadc_arg > 359:
            dadc_arg = 0.
        dadc_arg = np.clip(dadc_arg, 0., 180.)

        scale    = interp1d(dadc_array,self.scale[sorted_indices],'cubic')(dadc_arg)
        rotation = interp1d(dadc_array,self.rotation[sorted_indices],'cubic')(dadc_arg)
        offset_x = interp1d(dadc_array,self.offset_x[sorted_indices],'cubic')(dadc_arg)
        offset_y = interp1d(dadc_array,self.offset_y[sorted_indices],'cubic')(dadc_arg)
        zbcoeffs = np.zeros(self.zbcoeffs.shape[1])
        for i in range(self.zbcoeffs.shape[1]) :
            zbcoeffs[i] = interp1d(dadc_array,self.zbcoeffs[sorted_indices,i],'cubic')(dadc_arg)

        return scale,rotation,offset_x,offset_y,zbcoeffs

    def tan2fp(self, xtan, ytan, adc1, adc2):
        """
        Converts tangent plane coordinates xtan,ytan -> focal plane xfp,yfp
        """
        scale,rotation,offset_x,offset_y,zbcoeffs = self.interpolate_coeffs(adc1,adc2)

        mean_adc_rad = np.deg2rad(average_angles_deg(adc1, adc2))
        ca = np.cos(mean_adc_rad)
        sa = np.sin(mean_adc_rad)

        rxtan, rytan = _reduce_xytan(xtan,ytan)

        # rotate then derotate to account
        # for average ADC angle
        # this leave 20 microns residuals, I don't understand why ...
        rrxtan = ca*rxtan + sa*rytan
        rrytan = -sa*rxtan + ca*rytan
        rrxfp, rryfp   = transform(rrxtan, rrytan, scale, rotation, offset_x, offset_y, self.zbpolids, zbcoeffs)
        rxfp = ca*rrxfp - sa*rryfp
        ryfp = sa*rrxfp + ca*rryfp

        xfp, yfp     = _expand_xyfp(rxfp, ryfp)

        return xfp, yfp

    def fp2tan(self, xfp, yfp, adc1, adc2):
        """
        Converts focal plane xfp,yfp -> tangent plane xtan,ytan
        """

        scale,rotation,offset_x,offset_y,zbcoeffs = self.interpolate_coeffs(adc1, adc2)

        mean_adc_rad = np.deg2rad(average_angles_deg(adc1, adc2))
        ca = np.cos(mean_adc_rad)
        sa = np.sin(mean_adc_rad)

        rxfp, ryfp = _reduce_xyfp(xfp, yfp)

        # average ADC rotation
        rrxfp = ca*rxfp + sa*ryfp
        rryfp = -sa*rxfp + ca*ryfp

        #- first undo Zhao-Burge terms
        #- Iteratively find the correction, since we aren't interested
        #- in the correction at rxfp,ryfp but rather the correction at
        #- a different rx,ry that when applies becomes rxfp, ryfp
        dx = dy = 0.0
        for _ in range(20):
            dx2, dy2 = getZhaoBurgeXY(self.zbpolids, zbcoeffs, rrxfp-dx, rryfp-dy)
            dmax = max(np.max(np.abs(dx2-dx)), np.max(np.abs(dy2-dy)))
            dx, dy = dx2, dy2
            if dmax < 1e-12:
                break

        rrxfp -= dx
        rryfp -= dy

        #- Then apply inverse scale, roation, offset
        rrxfp /= scale
        rryfp /= scale

        rxx = (rrxfp*np.cos(-rotation) - rryfp*np.sin(-rotation))
        ryy = (rrxfp*np.sin(-rotation) + rryfp*np.cos(-rotation))

        rxx -= offset_x
        ryy -= offset_y

        # undo average ADC rotation
        xx = ca*rxx - sa*ryy
        yy = +sa*rxx + ca*ryy

        xtan, ytan = _expand_xytan(xx, yy)

        return xtan, ytan

# generic function
raytracefit_instance = None

def get_raytracefit() :
    global raytracefit_instance
    if raytracefit_instance is None :
        filename = resource_filename('desimeter', 'data/raytrace-tan2fp.json')
        raytracefit_instance = TAN2FP_RayTraceFit.read_jsonfile(filename=filename)
    return raytracefit_instance

def tan2fp(xtan, ytan, adc1, adc2):
    return get_raytracefit().tan2fp(xtan, ytan, adc1, adc2)


def fp2tan(xtan, ytan, adc1, adc2):
    return get_raytracefit().fp2tan(xtan, ytan, adc1, adc2)
