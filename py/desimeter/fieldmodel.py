"""
Utility functions to fit guide stars
"""

import json
import numpy as np

from desimeter.transform.radec2tan import xy2hadec,hadec2xy,tan2radec,radec2tan
from desimeter.transform.tan2fp import tan2fp,fp2tan
from desimeter.transform.gfa2fp import gfa2fp
from desimeter.trig import cosd, sind, arctan2d, rot2deg

from desimeter.log import get_logger
from astropy.table import Table
from astropy.time import Time

class FieldModel(object):

    def __init__(self) :

        self.ra  = None
        self.dec = None
        self.mjd = None
        self.lst = None
        self.hexrot_deg = None
        self.adc1 = None
        self.adc2 = None

        self.precession = True
        self.aberration = True
        self.polar_misalignment = True

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
        params['adc1'] = self.adc1
        params['adc2'] = self.adc2

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
        tx.adc1 = params['adc1']
        tx.adc2 = params['adc2']

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

        if "mjd_obs" in catalog.dtype.names :
            self.mjd = np.mean(catalog["mjd_obs"])
            log.info("use mjd={} from catalog['mjd_obs']".format(self.mjd))

        if ( not  "xcentroid" in catalog.dtype.names ) or ( not "ra_gaia" in catalog.dtype.names ) :
            log.error("I can only deal with Aaron's catalogs with columns xcentroid,ycentroid,ra_gaia,dec_gaia, sorry")
            raise RuntimeError("I can only deal with Aaron's catalogs with columns xcentroid,ycentroid,ra_gaia,dec_gaia, sorry")

        log.info("selection stars for which we have a good match (< {} arcsec)".format(max_sep_arcsec))

        if all([_ in catalog.columns for _ in ['pmra', 'pmdec', 'ref_epoch']]):

            if self.mjd is None :
                log.error("Cannot compute proper motion correction because mjd=None")
                raise RuntimeError("Cannot compute proper motion correction because mjd=None")

            # if proper motions and reference epochs are there
            pmra = catalog['pmra']
            pmdec = catalog['pmdec']
            # if unknown, zero out the pms
            pmra[~np.isfinite(pmra)] = 0
            pmdec[~np.isfinite(pmdec)] = 0

            cur_year = Time(self.mjd, format='mjd').to_value(format='jyear')
            # observation time in decimal years (like 2020.3)
            ref_epoch = catalog['ref_epoch']
            dra  = (cur_year - ref_epoch) * pmra / 3600e3 / cosd(catalog['dec_gaia'])
            ddec = (cur_year - ref_epoch) * pmdec / 3600e3

            # add pm and rename columns
            catalog['ra_gaia'] += dra
            catalog['dec_gaia'] += ddec
            catalog.rename_column('ra_gaia','ra_gaia_with_pm')
            catalog.rename_column('dec_gaia','dec_gaia_with_pm')
            ra_column = 'ra_gaia_with_pm'
            dec_column = 'dec_gaia_with_pm'

        else :
            ra_column = 'ra_gaia'
            dec_column = 'dec_gaia'
            log.warning("No proper motion info in catalog")

        match_dra  = (catalog["ra"]-catalog[ra_column])*cosd(catalog[dec_column]) * 3600. # arcsec
        match_ddec = (catalog["dec"]-catalog[dec_column]) * 3600. # arcsec

        dr = np.hypot(match_dra, match_ddec)
        selection = (dr<max_sep_arcsec) # arcsec
        if np.sum(selection)==0 :
            log.error("no star is matched with sufficient precision!")
            raise RuntimeError("no star is matched with sufficient precision!")

        return catalog[:][selection]

    def fit_tancorr(self,catalog,mjd=None,hexrot_deg=None,lst=None) :
        log = get_logger()

        x_gfa  = catalog["xcentroid"]
        y_gfa  = catalog["ycentroid"]
        if 'ra_gaia_with_pm' in catalog.dtype.names :
            ra_gaia  = catalog["ra_gaia_with_pm"]
            dec_gaia = catalog["dec_gaia_with_pm"]
        else :
            ra_gaia  = catalog["ra_gaia"]
            dec_gaia = catalog["dec_gaia"]

        # mjd,hexprot_deg,lst could have been set before
        if mjd is not None :
            self.mjd = mjd
            log.info("Use argument MJD={}".format(self.mjd))
        if self.mjd is None :
            if "mjd_obs" in catalog.keys():
                self.mjd    = np.mean(catalog["mjd_obs"])
                log.info("Use 'mjd_obs' in catalog, MJD={}".format(self.mjd))
            else :
                log.error("mjd is None")
                raise RuntimeError("mjd is None")

        log.info("Use MJD={}".format(self.mjd))

        if hexrot_deg is not None :
            self.hexrot_deg = hexrot_deg
        elif self.hexrot_deg is None :
            log.error("hexrot_deg is None")
            raise RuntimeError("hexrot_deg is None")

        if lst is not None :
            self.lst = lst
        elif self.lst is None :
            from desimeter.time import mjd2lst
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
        x_tan_meas,y_tan_meas = fp2tan(x_fp,y_fp,self.adc1,self.adc2)

        correction = TanCorr()

        for _ in range(3) : # loop because change of pointing induces a rotation of the field

            # we transform GAIA coordinates to the tangent plane
            x_tan_gaia,y_tan_gaia = radec2tan(ra_gaia,dec_gaia,self.ra,self.dec,mjd=self.mjd,lst_deg=self.lst,hexrot_deg = self.hexrot_deg, precession = self.precession, aberration = self.aberration, polar_misalignment = self.polar_misalignment)

            # now that we have both sets of coordinates, we fit a transformation from one to the other
            correction.fit(x_tan_meas,y_tan_meas,x_tan_gaia,y_tan_gaia)

            # opposite sign for the telescope offset because I have converted GAIA RA Dec to tangent plane ...
            self.dec -= correction.ddec
            self.ra  += correction.dha/cosd(self.dec) # HA = LST-RA


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
        x1 = np.append(0., np.deg2rad(np.cos(phi)))
        y1 = np.append(0., np.deg2rad(np.sin(phi)))

        # convert to sky
        xfp,yfp = tan2fp(x1,y1,self.adc1,self.adc2)
        ra,dec  = self.fp2radec(xfp,yfp)

        # vanilla transformation from ha,dec to tangent plane
        ha      = self.lst - ra
        x2,y2   = hadec2xy(ha,dec,ha[0],dec[0])

        return  np.rad2deg(np.mean((y1[1:]*x2[1:]-x1[1:]*y2[1:])/np.sqrt((x1[1:]**2+y1[1:]**2)*(x2[1:]**2+y2[1:]**2))))

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
            except KeyError:
                log.warning("missing metrology")

        return x_fp,y_fp

    def fp2radec(self,x_fp,y_fp) :
        x_tan,y_tan = fp2tan(x_fp,y_fp,self.adc1,self.adc2)
        x_tan,y_tan = self.tancorr_inst2sky(x_tan,y_tan) # correction
        ra,dec = tan2radec(x_tan,y_tan,self.ra,self.dec,mjd=self.mjd,lst_deg=self.lst,hexrot_deg = self.hexrot_deg, precession = self.precession, aberration = self.aberration, polar_misalignment = self.polar_misalignment)
        return ra,dec

    def radec2fp(self,ra,dec) :
        x_tan,y_tan = radec2tan(ra,dec,self.ra,self.dec,mjd=self.mjd,lst_deg=self.lst,hexrot_deg = self.hexrot_deg, precession = self.precession, aberration = self.aberration, polar_misalignment = self.polar_misalignment)
        x_tan,y_tan = self.tancorr_sky2inst(x_tan,y_tan) # correction
        x_fp,y_fp   = tan2fp(x_tan,y_tan,self.adc1,self.adc2)
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
        for _ in range(4):
            x1t,y1t  = hadec2xy(ha,dec,self.dha,self.ddec)
            dx = np.mean(x2-x1t)
            dy = np.mean(y2-y1t)
            self.dha  -= np.rad2deg(dx)
            self.ddec -= np.rad2deg(dy)
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
        rms = np.sqrt( np.mean( (x2-x2p)**2 + (y2-y2p)**2 ) )
        self.rms_arcsec = np.rad2deg(rms) * 3600.

        # interpret this back into telescope pointing offset, field rotation, dilatation

        # pointing offset
        # increasing gaia stars x means telescope is more to the left so tel_ha should be decreased
        # increasing gaia stars y means telescope is more to the bottom so tel_dec should be decreased
        # tangent plane coordinates are in rad
        ddha  = -np.rad2deg(ax[0])
        dddec = -np.rad2deg(ay[0])

        self.dha  += ddha
        self.ddec += dddec

        # dilatation and rotation
        # |ax1 ax2| |sxx sxy| |ca  -sa|
        # |ay1 ay2|=|syx syy|*|sa   ca|
        # ax1=sxx*ca+sxy*sa ; ax2=-sxx*sa+sxy*ca
        # ay1=syx*ca+syy*sa ; ay2=-syx*sa+syy*ca
        # ax1+ay2 = (sxx+syy)*ca
        # ay1-ax2 = (sxx+syy)*sa

        sxx_p_syy = np.hypot(ax[1] + ay[2], ay[1] - ax[2])
        sa=(ay[1]-ax[2])/sxx_p_syy
        ca=(ax[1]+ay[2])/sxx_p_syy

        self.rot_deg = arctan2d(sa,ca)

        sxy = sa*ax[1]+ca*ay[1] - sxx_p_syy*ca*sa
        sxx =(ax[1]-sxy*sa)/ca
        syy = (ay[1]-sxy*ca)/sa

        self.sxx = sxx
        self.syy = syy
        self.sxy = sxy

    def apply(self,x,y) :
        scale_matrix = np.array([[self.sxx,self.sxy],[self.sxy,self.syy]])
        rot_matrix = rot2deg(self.rot_deg)
        ha,dec  = xy2hadec(x,y,0,0)
        x1t,y1t = hadec2xy(ha,dec,self.dha,self.ddec)
        xy=scale_matrix.dot(rot_matrix.dot(np.array([x1t,y1t])))
        return xy[0],xy[1]

    def apply_inverse(self,x,y) :
        det = self.sxx*self.syy - self.sxy**2
        scale_matrix = np.array([[self.syy,-self.sxy],[-self.sxy,self.sxx]])/det
        rot_matrix = rot2deg(-self.rot_deg)
        ha,dec  = xy2hadec(x,y,0,0)
        x1t,y1t = hadec2xy(ha,dec,self.dha,self.ddec)
        xy=rot_matrix.dot(scale_matrix.dot(np.array([x1t,y1t])))
        return xy[0],xy[1]

def fieldrot(ra,dec,mjd,lst_deg,hexrot_deg=0) :
    """
    Computes the field rotation in degrees.

    Args:
      ra: scalar, float, RA in degrees
      dec: scalar, float, Dec in degrees
      mjd: scalar, float, MJD in days
      lst_deg: scalar, float, LST in degrees (ha = lst_deg - ra)

    Optional:
      hexrot_deg: scalar, float, Hexapod rotation in degrees (default=0)

    Returns:
      field rotation angle in degrees

    """
    fm = FieldModel()
    fm.ra=ra
    fm.dec=dec
    fm.mjd=mjd
    fm.lst=lst_deg
    fm.adc1=0
    fm.adc2=0
    fm.hexrot_deg=hexrot_deg
    return fm.compute_fieldrot()

def dfieldrotdt_physical_model(ra,dec,mjd,lst_deg) :
    """
    Computes the derivative with time of the field rotation in arcsec per minute.

    Args:
      ra: scalar, float, RA in degrees
      dec: scalar, float, Dec in degrees
      mjd: scalar, float, MJD in days
      lst_deg: scalar, float, LST in degrees (ha = lst_deg - ra)

    Returns:
      field rotation angle derivative with time in arcsec per minute

    """
    one_minute_in_degrees = 1./60./24.*360 #
    return 3600.*(fieldrot(ra,dec,mjd,lst_deg+one_minute_in_degrees,hexrot_deg=0) - fieldrot(ra,dec,mjd,lst_deg,hexrot_deg=0))


def dfieldrotdt_empirical_model(ra,dec,lst_deg) :
    """
    Computes the derivative with time of the field rotation in arcsec per minute.

    Args:
      ra: scalar, float, RA in degrees
      dec: scalar, float, Dec in degrees
      lst_deg: scalar, float, LST in degrees (ha = lst_deg - ra)

    Returns:
      field rotation angle derivative with time in arcsec per minute

    """



    # Empirical polynomial model
    # based on the output of the script desi_fit_field_rotation_rate.
    # Result of fit on data/guide_data_20200415.csv + guide_20201215_with_hexrot.csv

    scalar = np.isscalar(ra)
    if scalar :
        ra  = np.atleast_1d(ra)
        dec = np.atleast_1d(dec)

    if np.any(np.abs(dec)>90.) :
        raise ValueError("Unphysical Dec in degrees: {}".format(dec))

    ha = (lst_deg - ra)
    ha = ha%360.
    ha[ha>180.] -= 360.

    x=ha/60.
    y=(dec-30)/30.

    rotation_rate_arcsec_per_min = -0.396 + 0.012*x + 0.012*y + 0.124*x*y + 0.439*x**2 - 0.144*y**2 + 0.012*x**3 - 0.021*x**2*y + 0.095*x*y**2 -0.062*y**3

    # saturation to avoid crazy values if we run outside of the validity range of the polynomial
    min_val = -1.5
    max_val = 1.5
    rotation_rate_arcsec_per_min[rotation_rate_arcsec_per_min<min_val] = min_val
    rotation_rate_arcsec_per_min[rotation_rate_arcsec_per_min>max_val] = max_val

    if scalar :
        return rotation_rate_arcsec_per_min[0]
    else :
        return rotation_rate_arcsec_per_min


def dfieldrotdt(ra,dec,mjd,lst_deg) :
    """
    Computes the derivative with time of the field rotation in arcsec per minute.

    Args:
      ra: scalar, float, RA in degrees
      dec: scalar, float, Dec in degrees
      mjd: scalar, float, MJD in days
      lst_deg: scalar, float, LST in degrees (ha = lst_deg - ra)

    Returns:
      field rotation angle derivative with time in arcsec per minute

    """

    # choose to use the totally empirical model, ignore mjd
    return dfieldrotdt_empirical_model(ra,dec,lst_deg)
