
"""
Utility functions to fit guide stars
"""

import json
import numpy as np

from desimeter.transform.radec2tan import xy2hadec,hadec2xy,tan2radec,radec2tan
from desimeter.transform.tan2fp import tan2fp,fp2tan
from desimeter.transform.gfa2fp import gfa2fp

from desimeter.log import get_logger
from astropy.table import Table

class FieldModel(object):

    def __init__(self) :

        self.ra  = None
        self.dec = None
        self.mjd = None
        self.lst = None
        self.hexrot_deg = None
        
        self.sxx = 1.
        self.syy = 1.
        self.sxy = 0.
        self.fieldrot_zp_deg = 0. # zero point of field rotation (should be constant across exposures)
        self.fieldrot_deg = 0. # rotation of field as defined in DESI-5190, to compare to fiberassign FIELDROT 
        self.expid  = 0
        self.nstars = 0
        self.rms_arcsec = 0.
    
    def tojson(self):
        params = dict()
        params['name'] = "Field model"
        params['version'] = '1'
        params['ra'] = self.ra
        params['dec'] = self.dec
        params['mjd'] = self.mjd
        params['lst'] = self.lst
        params['hexrot_deg'] = self.hexrot_deg
        
        params['sxx'] = self.sxx
        params['syy'] = self.syy
        params['sxy'] = self.sxy
        params['fieldrot_zp_deg'] = self.fieldrot_zp_deg        
        params['fieldrot_deg'] = self.fieldrot_deg        
        params['expid']=self.expid
        params['nstars'] = self.nstars
        params['rms_arcsec'] = self.rms_arcsec        
        return json.dumps(params)

    @classmethod
    def fromjson(cls, jsonstring):
        tx = cls()
        params = json.loads(jsonstring)
        assert params['name'] == "Field model"
        assert params['version'] == '1'
        tx.ra = params['ra']
        tx.dec = params['dec']
        tx.mjd = params['mjd']
        tx.lst = params['lst']
        tx.hexrot_deg = params['hexrot_deg']
        
        tx.sxx = params['sxx']
        tx.syy = params['syy']
        tx.sxy = params['sxy']
        tx.fieldrot_zp_deg = params['fieldrot_zp_deg']
        tx.fieldrot_deg = params['fieldrot_deg']
        tx.expid = params['expid']
        tx.nstars = params['nstars']
        tx.rms_arcsec = params['rms_arcsec']        
        return tx

    
    def read_guide_stars_catalog(self,filename,max_sep_arcsec = 2.) :

        log = get_logger()
        
        log.info("reading guide stars in {}".format(filename))

        # here we could do some conversion of column names
        catalog = Table.read(filename)

        if ( not  "xcentroid" in catalog.dtype.names ) or ( not "ra_gaia" in catalog.dtype.names ) :
            log.error("I can only deal with Aaron's catalogs with columns xcentroid,ycentroid,ra_gaia,dec_gaia, sorry")
            raise RuntimeError("I can only deal with Aaron's catalogs with columns xcentroid,ycentroid,ra_gaia,dec_gaia, sorry")

        max_sep_arcsec = 2.
        log.info("selection stars for which we have a good match (< {} arcsec)".format(max_sep_arcsec))
        
        dra  = (catalog["ra"]-catalog["ra_gaia"])*np.cos(catalog["dec_gaia"]/180*np.pi)*3600. # arcsec
        ddec = (catalog["dec"]-catalog["dec_gaia"])*3600. # arcsec
        dr = np.sqrt(dra**2+ddec**2)
        selection = (dr<2) # arcsec
        if np.sum(selection)==0 :
            log.error("no star is matched with sufficient precision!")
            raise RuntimeError("no star is matched with sufficient precision!")

        return catalog[:][selection]


    
    def fit_tancorr(self,catalog,mjd=None,hexrot_deg=None,lst=None) :

        log = get_logger()
        
        x_gfa  = catalog["xcentroid"]
        y_gfa  = catalog["ycentroid"]
        ra_gaia  = catalog["ra_gaia"]
        dec_gaia = catalog["dec_gaia"]
        
        # mjd,hexprot_deg,lst could have been set before
        if mjd is not None :
            self.mjd = mjd
            log.info("Use argument MJD={}".format(self.mjd))
        elif "mjd_obs" in catalog.keys():
            self.mjd    = np.mean(catalog["mjd_obs"])
            log.info("Use 'mjd_obs' in catalog, MJD={}".format(self.mjd))
        elif self.mjd is None :
            log.error("mjd is None")
            raise RuntimeError("mjd is None")
        else :
            log.info("Use MJD={}".format(self.mjd))
        
        if hexrot_deg is not None :
            self.hexrot_deg = hexrot_deg
        elif self.hexrot_deg is None :
            log.error("hexrot_deg is None")
            raise RuntimeError("hexrot_deg is None")

        if lst is not None :
            self.lst = lst
        elif self.lst is None :
            log.warning("Compute LST from MJD={}".format(self.mjd))
            self.lst = mjd2lst(self.mjd)
        log.info("Use LST={}".format(self.lst))
            
        # first transfo: gfa2fp
        x_fp,y_fp = self.all_gfa2fp(x_gfa,y_gfa,petal_loc=catalog["petal_loc"])
        
        # keep only petal data for which we have the metrology        
        selection = (x_fp!=0)
        x_gfa    = x_gfa[selection]
        y_gfa    = y_gfa[selection]
        x_fp     = x_fp[selection]
        y_fp     = y_fp[selection]
        ra_gaia  = ra_gaia[selection]
        dec_gaia = dec_gaia[selection]

        # transform focal plane to tangent plane
        x_tan_meas,y_tan_meas = fp2tan(x_fp,y_fp)

        correction = TanCorr()

        for loop in range(3) : # loop because change of pointing induces a rotation of the field
            
            # we transform GAIA coordinates to the tangent plane
            x_tan_gaia,y_tan_gaia = radec2tan(ra_gaia,dec_gaia,self.ra,self.dec,mjd=self.mjd,lst_deg=self.lst,hexrot_deg = self.hexrot_deg)
        
            # now that we have both sets of coordinates, we fit a transformation from one to the other
            correction.fit(x_tan_meas,y_tan_meas,x_tan_gaia,y_tan_gaia)

            # opposite sign for the telescope offset because I have converted GAIA RA Dec to tangent plane ...
            self.dec -= correction.ddec
            self.ra  += correction.dha/np.cos(self.dec*np.pi/180.) # HA = LST-RA

        
        # save params to this
        log.info("RMS coord. residual = {:3.2f} arcsec".format(correction.rms_arcsec))
        log.info("Rotation angle (field rot ZP) ={:4.3f} deg".format(correction.rot_deg))
        log.info("Pointing correction dHA={:3.2f} arcsec, dDec={:3.2f} arcsec".format(correction.dha*3600.,correction.ddec*3600.))
        log.info("Scales sxx={:5.4f} syy={:5.4f} sxy={:5.4f}".format(correction.sxx,correction.syy,correction.sxy))

        # now I just copy the correction parameters in this class
        self.sxx = correction.sxx
        self.syy = correction.syy
        self.sxy = correction.sxy
        self.fieldrot_zp_deg = correction.rot_deg
        self.nstars = correction.nstars
        self.rms_arcsec = correction.rms_arcsec

        # I now derive the field rotation
        self.fieldrot_deg = self.compute_fieldrot()


    def compute_fieldrot(self) :

        # cross of side length 1 degree in tangent plane
        phi = np.arange(4)*np.pi/2.
        x1 = np.append(0.,np.pi/180.*np.cos(phi))
        y1 = np.append(0.,np.pi/180.*np.sin(phi))

        # convert to sky 
        xfp,yfp = tan2fp(x1,y1)
        ra,dec  = self.fp2radec(xfp,yfp)
        
        # vanilla transformation from ha,dec to tangent plane
        ha      = self.lst - ra
        x2,y2   = hadec2xy(ha,dec,ha[0],dec[0])

        return  180./np.pi*np.mean((y1[1:]*x2[1:]-x1[1:]*y2[1:])/np.sqrt((x1[1:]**2+y1[1:]**2)*(x2[1:]**2+y2[1:]**2)))
            
    # tangent plane correction from the instrument to the sky
    def tancorr_inst2sky(self,x_tan,y_tan) :
        t2t = TanCorr()
        t2t.dha  = 0. # because already included in self.ra
        t2t.ddec = 0. # because already included in self.dec
        # dilatation params
        t2t.sxx = self.sxx
        t2t.syy = self.syy
        t2t.sxy = self.sxy
        # rotation params
        t2t.rot_deg = self.fieldrot_zp_deg
        return t2t.apply(x_tan,y_tan)

    # tangent plane correction from the sky to the instrument
    # this is the inverse of the recorded values
    def tancorr_sky2inst(self,x_tan,y_tan) :
        t2t = TanCorr()
        t2t.dha  = 0. # because already included in self.ra
        t2t.ddec = 0. # because already included in self.dec
        # dilatation params
        t2t.sxx = self.sxx
        t2t.syy = self.syy
        t2t.sxy = self.sxy
        # rotation params
        t2t.rot_deg = self.fieldrot_zp_deg
        return t2t.apply_inverse(x_tan,y_tan)
        
    def all_gfa2fp(self,x_gfa,y_gfa,petal_loc) :

        log = get_logger()
        
        # simple loop on petals
        petals = np.unique(petal_loc)
        x_fp = np.zeros(x_gfa.shape)
        y_fp = np.zeros(y_gfa.shape)
        for petal in petals :
            ii = (petal_loc==petal)
            try :
                x,y = gfa2fp(petal,x_gfa[ii],y_gfa[ii])
                x_fp[ii]  = x
                y_fp[ii]  = y
            except KeyError as e :
                log.warning("missing metrology")
                pass
        return x_fp,y_fp
    
    def fp2radec(self,x_fp,y_fp) :
        x_tan,y_tan = fp2tan(x_fp,y_fp)
        x_tan,y_tan = self.tancorr_inst2sky(x_tan,y_tan) # correction
        ra,dec = tan2radec(x_tan,y_tan,self.ra,self.dec,mjd=self.mjd,lst_deg=self.lst,hexrot_deg = self.hexrot_deg)
        return ra,dec

    def radec2fp(self,ra,dec) :
        x_tan,y_tan = radec2tan(ra,dec,self.ra,self.dec,mjd=self.mjd,lst_deg=self.lst,hexrot_deg = self.hexrot_deg)
        x_tan,y_tan = self.tancorr_sky2inst(x_tan,y_tan) # correction
        x_fp,y_fp   = tan2fp(x_tan,y_tan)
        return x_fp,y_fp

    
    
class TanCorr(object):

    def __init__(self) :
        self.dha  = 0.
        self.ddec = 0.
        self.sxx = 1.
        self.syy = 1.
        self.sxy = 0.
        self.rot_deg = 0.
        self.nstars = 0
        self.rms_arcsec = 0.
    
    def tojson(self):
        params = dict()
        params['name'] = "Tangent Plane Correction"
        params['version'] = '1'
        params['dha'] = self.dha
        params['ddec'] = self.ddec
        params['sxx'] = self.sxx
        params['syy'] = self.syy
        params['sxy'] = self.sxy
        params['rot_deg'] = self.rot_deg        
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
        tx.rot_deg = params['rot_deg']
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
        ddha  = -ax[0]*180./np.pi
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
        ha,dec  = xy2hadec(x,y,0,0)
        x1t,y1t = hadec2xy(ha,dec,self.dha,self.ddec)
        xy=scale_matrix.dot(rot_matrix.dot(np.array([x1t,y1t])))
        return xy[0],xy[1]

    def apply_inverse(self,x,y) :

        det = self.sxx*self.syy - self.sxy**2
        scale_matrix = np.array([[self.syy,-self.sxy],[-self.sxy,self.sxx]])/det
        
        ca=np.cos(self.rot_deg/180*np.pi)
        sa=np.sin(self.rot_deg/180*np.pi)
        rot_matrix = np.array([[ca,sa],[-sa,ca]])
        ha,dec  = xy2hadec(x,y,0,0)
        x1t,y1t = hadec2xy(ha,dec,self.dha,self.ddec)
        xy=rot_matrix.dot(scale_matrix.dot(np.array([x1t,y1t])))
        return xy[0],xy[1]
        
        
        
