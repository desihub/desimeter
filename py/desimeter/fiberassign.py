import numpy as np

from desimeter.transform.radec2tan import radec2tan,hadec2xy,hadec2altaz,tan2radec,LATITUDE
from desimeter.transform.tan2fp.raytracefit import tan2fp,fp2tan
from desimeter.trig import sincosd

def _measure_fieldrot_deg(ha,dec,tel_ha,tel_dec,xfp_mm,yfp_mm) :

    ok = (~(np.isnan(ha*dec*xfp_mm*yfp_mm)))&(xfp_mm**2+yfp_mm**2>10**2)
    x2,y2=hadec2xy(ha[ok],dec[ok],tel_ha,tel_dec) # rad
    return np.rad2deg(np.mean((yfp_mm[ok]*x2-xfp_mm[ok]*y2)/np.sqrt((xfp_mm[ok]**2+yfp_mm[ok]**2)*(x2**2+y2**2))))


# This comes straight from PlateMaker: python/PlateMaker/astron.py
# And it gives a value that = -1 * the DESI 'PARALLAC' header card.
def parallactic_angle(ha, decl):
    """Calculate the parallactic angle.
        Args:
        ha: Hour Angle (decimal degrees)
        decl: declination (decimal degrees)
        Returns:
        the parallactic angle in decimal degrees
    """
    rha = np.deg2rad(ha)
    rdecl = np.deg2rad(decl)
    rphi = np.deg2rad(LATITUDE)
    rpsi = -1 * np.arctan2(np.sin(rha) * np.cos(rphi),
                           np.cos(rdecl) * np.sin(rphi) - np.sin(rdecl) * np.cos(rphi) * np.cos(rha))
    return np.rad2deg(rpsi)

def zd2deltaadc(zd):
        # This comes from the PlateMaker code --
        #  https://desi.lbl.gov/trac/browser/code/online/DervishTools/trunk/desi/etc/desi/PRISM.par
        #  https://desi.lbl.gov/trac/browser/code/online/DervishTools/trunk/desi/etc/common.tcl#L839
        t = np.tan(np.deg2rad(zd))
        A = -0.0183 + -0.3795*t + -0.1939*t**2
        A = np.rad2deg(A)
        return 2*A

# Direct copy from https://github.com/desihub/fiberassign/blob/radec2xy/py/fiberassign/targets.py
# These are reproductions of PlateMaker Tcl functions in
# https://desi.lbl.gov/trac/browser/code/online/DervishTools/trunk/desi/etc/nfs.tcl#L43
# def pm_sidtim(mjd):
#     return fmod(mjd + (mjd-52903.54875)*(366.24125/365.24125-1.),1.)*360.
def pm_zd(ha, dec, lat):
    rha = np.deg2rad(ha)
    rdec = np.deg2rad(dec)
    rphi = np.deg2rad(lat)
    rzen = np.arccos(np.cos(rha) * np.cos(rdec) * np.cos(rphi) +
            np.sin(rdec) * np.sin(rphi))
    return np.rad2deg(rzen)

def pm_psi(ha, dec, latitude):
    rha = np.deg2rad(ha)
    rdec = np.deg2rad(dec)
    rphi = np.deg2rad(latitude)
    rpsi = np.arctan2(np.sin(rha) * np.cos(rphi), np.cos(rdec) * np.sin(rphi) - np.sin(rdec) * np.cos(rphi) * np.cos(rha))
    return np.rad2deg(rpsi)

def pm_zd2deltaadc(zd):
    t = np.tan(np.deg2rad(zd))
    A = -0.0183 + -0.3795*t + -0.1939*t**2
    A = np.rad2deg(A)
    return 2*A

def pm_get_adc_angles(ha, dec):
    # Here we're reproducing PlateMaker's astrometric
    # transformations to get to the ADC angles it's going to set.
    # These have degree-level disagreements with other methods.
    pm_longitude = 111.6003
    pm_latitude  =  31.9634
    #lst0 = pm_sidtim(mjd)
    #st = lst0 - longitude
    #ha = st - tile_ra
    zd  = pm_zd (ha, dec, pm_latitude)
    psi = pm_psi(ha, dec, pm_latitude)
    dadc = pm_zd2deltaadc(zd)
    adc1 = psi + dadc/2.
    adc2 = psi - dadc/2.
    return adc1, adc2

def fiberassign_radec2xy_cs5(ra,dec,tile_ra,tile_dec,tile_mjd,tile_ha,tile_fieldrot,adc1=None,adc2=None) :

    # LST from HA
    lst=tile_ha+tile_ra

    if ( (adc1 is None) and (adc2 is not None) ) or ( (adc1 is not None) and (adc2 is None) ) :
        raise ValueError("either set both adc angles or none")

    if adc1 is None :
        adc1,adc2 =  pm_get_adc_angles(tile_ha,tile_dec)

    # start with pointing = tile center
    tel_ra=tile_ra+0.
    tel_dec=tile_dec+0.

    # tune telescope pointing given ADC angle
    # in order to fp coordinates of tile RA and DEC, it's not zero because of the ADC angles
    for iteration in range(2) :

        xtan,ytan = radec2tan(np.array([tile_ra]),np.array([tile_dec]),tel_ra,tel_dec,tile_mjd,lst,hexrot_deg=0)
        xfp_0,yfp_0   = tan2fp(xtan,ytan,adc1,adc2) #mm
        #print("Temp tile center in FP coordinates = {},{} mm".format(xfp_0[0],yfp_0[0]))

        # numeric derivative
        eps = 1./3600. #
        xtan,ytan = radec2tan(np.array([tile_ra+eps]),np.array([tile_dec]),tel_ra,tel_dec,tile_mjd,lst,hexrot_deg=0)
        xfp_dra,yfp_dra   = tan2fp(xtan,ytan,adc1,adc2) #mm
        xtan,ytan = radec2tan(np.array([tile_ra]),np.array([tile_dec+eps]),tel_ra,tel_dec,tile_mjd,lst,hexrot_deg=0)
        xfp_ddec,yfp_ddec   = tan2fp(xtan,ytan,adc1,adc2) #mm
        dxdra=(xfp_dra[0]-xfp_0[0])/eps
        dydra=(yfp_dra[0]-yfp_0[0])/eps
        dxddec=(xfp_ddec[0]-xfp_0[0])/eps
        dyddec=(yfp_ddec[0]-yfp_0[0])/eps
        J=[[dxdra,dxddec],[dydra,dyddec]]

        # solve linear system to get tile RA Dec at center of fov
        Jinv=np.linalg.inv(J)
        X=Jinv.dot([xfp_0[0],yfp_0[0]])
        dra=X[0]
        ddec=X[1]

        # apply offset to telescope pointing
        tel_ra += dra
        tel_dec += ddec

    # verify
    xtan,ytan = radec2tan(np.array([tile_ra]),np.array([tile_dec]),tel_ra,tel_dec,tile_mjd,lst,hexrot_deg=0)
    xfp_0,yfp_0   = tan2fp(xtan,ytan,adc1,adc2) #mm
    print("Tile center in FP coordinates = {},{} mm".format(xfp_0[0],yfp_0[0]))

    # now compute coordinates of all targets
    xtan,ytan = radec2tan(ra,dec,tel_ra,tel_dec,tile_mjd,lst,hexrot_deg=0)
    tmp_xfp,tmp_yfp = tan2fp(xtan,ytan,adc1,adc2)

    # measure field rotation
    tmp_fieldrot = _measure_fieldrot_deg(-ra,dec,-tile_ra,tile_dec,tmp_xfp,tmp_yfp)

    # apply field rotation to match request
    drot = tile_fieldrot-tmp_fieldrot
    s,c = sincosd(drot)
    xfp = c * tmp_xfp - s * tmp_yfp
    yfp = s * tmp_xfp + c * tmp_yfp

    # verify
    realised_fieldrot = _measure_fieldrot_deg(-ra,dec,-tile_ra,tile_dec,xfp,yfp)
    print("Requested fieldrot={:3.1f} arcsec delta={:3.1f} arcsec".format(tile_fieldrot*3600.,(tile_fieldrot-realised_fieldrot)*3600.))

    return xfp,yfp

def fiberassign_radec2xy_flat(ra,dec,tile_ra,tile_dec,tile_mjd,tile_ha,tile_fieldrot,adc1,adc2) :

    xfp,yfp = fiberassign_radec2xy_cs5(ra,dec,tile_ra,tile_dec,tile_mjd,tile_ha,tile_fieldrot,adc1,adc2)

    # fiber assign coordinates are on the curved coordinates that follow the curved focal surface
    # the curved coordinates are called 'flat' in the focalplane parlance.
    xflat,yflat = ptl2flat(xfp,yfp)

    return xflat,yflat
