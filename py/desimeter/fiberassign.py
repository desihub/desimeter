"""
desimeter.fiberassign provides ra,dec <-> x,y transforms for fiberassign
which must predict these without using as-observed GFA data for final scale
and rotations.  See bin/plot_pm_coords for example using these functions.

The primary user-facing functions are:

  * Transforms to/from CS5 coordinate system:
    - fiberassign_radec2xy_cs5
    - fiberassign_cs5_xy2radec
  * Transforms to/from "flattened" coordinate system, treating the curved
    focal surface as if it was flat:
    - fiberassign_radec2xy_flat
    - fiberassign_flat_xy2radec
"""

import numpy as np

from desimeter.io import dm2pm_filename
from desimeter.transform.radec2tan import radec2tan,tan2radec,hadec2xy
from desimeter.transform.tan2fp.raytracefit import tan2fp,fp2tan
from desimeter.transform.pos2ptl import ptl2flat,flat2ptl
from desimeter.transform.dm2pm import DM2PM
from desimeter.trig import sincosd
import desimeter.log


log = desimeter.log.get_logger()


def _measure_fieldrot_deg(ha,dec,tel_ha,tel_dec,xfp_mm,yfp_mm) :

    ok = (~(np.isnan(ha*dec*xfp_mm*yfp_mm)))&(xfp_mm**2+yfp_mm**2>10**2)
    x2,y2=hadec2xy(ha[ok],dec[ok],tel_ha,tel_dec) # rad
    return np.rad2deg(np.mean((yfp_mm[ok]*x2-xfp_mm[ok]*y2)/np.sqrt((xfp_mm[ok]**2+yfp_mm[ok]**2)*(x2**2+y2**2))))


# This comes straight from PlateMaker: python/PlateMaker/astron.py
# And it gives a value that = -1 * the DESI 'PARALLAC' header card.

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
    zd = np.minimum(zd, 60.) # saturate zenith angle
    t = np.tan(np.deg2rad(zd))
    A = -0.0183 + -0.3795*t + -0.1939*t**2
    A = np.rad2deg(A)
    return 2*A

def pm_get_adc_angles(ha, dec):
    # Here we're reproducing PlateMaker's astrometric
    # transformations to get to the ADC angles it's going to set.
    # These have degree-level disagreements with other methods.
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

def fiberassign_radec2xy_cs5(ra,dec,tile_ra,tile_dec,tile_mjd,tile_ha,tile_fieldrot,adc1=None,adc2=None,to_platemaker=True, use_hardcoded_polmis_rotmat=False) :
    """Computes X Y focal plane coordinates of targets in CS5 coordinate system.
    Args:
      ra  : 1D numpy.array with target RA in degrees (need an array of points to compute the field rotation)
      dec : 1D numpy.array with target Dec in degrees, same size as ra
      tile_ra : float, RA of tile center in degrees, such that after pointing and ADC shift corrections, this corresponds to X=0,Y=0
      tile_dec : float, Dec of tile center in degrees, such that after pointing and ADC shift corrections, this corresponds to X=0,Y=0
      tile_mjd : float, days, used for precession (does not need to be super accurate)
      tile_ha :  float, hour angle of observation, in degrees
      tile_fieldrot :  float, design/requested field rotation, in degrees
      adc1 : optional, float, ADC1 angle in degrees: by default it is computed based on the HA and Dec.
      adc2 : optional, float, ADC2 angle in degrees: by default it is computed based on the HA and Dec.
      to_platemaker : assume output coordinates are in platemaker system and use dm2pm transform
      use_hardcoded_polmis_rotmat (optional, defaults to False): use pre-computed polar misalignment matrix, in desimeter.transform.radec2tan.get_hardcoded_polmis_rotmat() (bool)
    Returns xfp,yfp, focal plane coordinates in mm in CS5, 1D numpy arrays of same size as ra,dec
    """

    assert(ra.size == dec.size)

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
    for _ in range(2) :

        xtan,ytan = radec2tan(np.array([tile_ra]),np.array([tile_dec]),tel_ra,tel_dec,tile_mjd,lst,hexrot_deg=0, use_hardcoded_polmis_rotmat=use_hardcoded_polmis_rotmat)
        xfp_0,yfp_0   = tan2fp(xtan,ytan,adc1,adc2) #mm
        #log.info("Temp tile center in FP coordinates = {},{} mm".format(xfp_0[0],yfp_0[0]))

        # numeric derivative
        eps = 1./3600. #
        xtan,ytan = radec2tan(np.array([tile_ra+eps]),np.array([tile_dec]),tel_ra,tel_dec,tile_mjd,lst,hexrot_deg=0, use_hardcoded_polmis_rotmat=use_hardcoded_polmis_rotmat)
        xfp_dra,yfp_dra   = tan2fp(xtan,ytan,adc1,adc2) #mm
        xtan,ytan = radec2tan(np.array([tile_ra]),np.array([tile_dec+eps]),tel_ra,tel_dec,tile_mjd,lst,hexrot_deg=0, use_hardcoded_polmis_rotmat=use_hardcoded_polmis_rotmat)
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
    xtan,ytan = radec2tan(np.array([tile_ra]),np.array([tile_dec]),tel_ra,tel_dec,tile_mjd,lst,hexrot_deg=0, use_hardcoded_polmis_rotmat=use_hardcoded_polmis_rotmat)
    xfp_0,yfp_0   = tan2fp(xtan,ytan,adc1,adc2) #mm
    log.info("Tile center in FP coordinates = {},{} mm".format(xfp_0[0],yfp_0[0]))

    # now compute coordinates of all targets
    xtan,ytan = radec2tan(ra,dec,tel_ra,tel_dec,tile_mjd,lst,hexrot_deg=0, use_hardcoded_polmis_rotmat=use_hardcoded_polmis_rotmat)
    tmp_xfp,tmp_yfp = tan2fp(xtan,ytan,adc1,adc2)

    if to_platemaker :
        # apply tranformation from desimeter to platemater
        dm2pm = DM2PM.read(dm2pm_filename())
        tmp_xfp,tmp_yfp = dm2pm.dm2pm(tmp_xfp,tmp_yfp)

    # measure field rotation
    tmp_fieldrot = _measure_fieldrot_deg(-ra,dec,-tile_ra,tile_dec,tmp_xfp,tmp_yfp)

    # apply field rotation to match request
    drot = tile_fieldrot-tmp_fieldrot
    s,c = sincosd(drot)
    xfp = c * tmp_xfp - s * tmp_yfp
    yfp = s * tmp_xfp + c * tmp_yfp

    # verify
    realised_fieldrot = _measure_fieldrot_deg(-ra,dec,-tile_ra,tile_dec,xfp,yfp)
    log.info("Requested fieldrot={:3.1f} arcsec delta={:3.1f} arcsec".format(tile_fieldrot*3600.,(tile_fieldrot-realised_fieldrot)*3600.))

    return xfp,yfp

def fiberassign_cs5_xy2radec(xfp,yfp,tile_ra,tile_dec,tile_mjd,tile_ha,tile_fieldrot,adc1=None,adc2=None,from_platemaker=True, use_hardcoded_polmis_rotmat=False) :
    """Computes RA Dec from focal plane coordinates of targets in CS5 coordinate system (inverse of fiberassign_radec2xy_cs5)
    Args:
      xfp  : 1D numpy.array with target X in mm (need an array of points to compute the field rotation)
      yfp : 1D numpy.array with target Y in mm
      tile_ra : float, RA of tile center in degrees, such that after pointing and ADC shift corrections, this corresponds to X=0,Y=0
      tile_dec : float, Dec of tile center in degrees, such that after pointing and ADC shift corrections, this corresponds to X=0,Y=0
      tile_mjd : float, days, used for precession (does not need to be super accurate)
      tile_ha :  float, hour angle of observation, in degrees
      tile_fieldrot :  float, design/requested field rotation, in degrees
      adc1 : optional, float, ADC1 angle in degrees: by default it is computed based on the HA and Dec.
      adc2 : optional, float, ADC2 angle in degrees: by default it is computed based on the HA and Dec.
      from_platemaker : assume input coordinates are from platemaker and use pm2dm transform
      use_hardcoded_polmis_rotmat (optional, defaults to False): use pre-computed polar misalignment matrix, in desimeter.transform.radec2tan.get_hardcoded_polmis_rotmat() (bool)
    Returns RA,Dec, 1D numpy arrays of same size as X,Y
    """

    assert(xfp.size == yfp.size)

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
    for _ in range(2) :

        xtan,ytan = radec2tan(np.array([tile_ra]),np.array([tile_dec]),tel_ra,tel_dec,tile_mjd,lst,hexrot_deg=0, use_hardcoded_polmis_rotmat=use_hardcoded_polmis_rotmat)
        xfp_0,yfp_0   = tan2fp(xtan,ytan,adc1,adc2) #mm
        #log.info("Temp tile center in FP coordinates = {},{} mm".format(xfp_0[0],yfp_0[0]))

        # numeric derivative
        eps = 1./3600. #
        xtan,ytan = radec2tan(np.array([tile_ra+eps]),np.array([tile_dec]),tel_ra,tel_dec,tile_mjd,lst,hexrot_deg=0, use_hardcoded_polmis_rotmat=use_hardcoded_polmis_rotmat)
        xfp_dra,yfp_dra   = tan2fp(xtan,ytan,adc1,adc2) #mm
        xtan,ytan = radec2tan(np.array([tile_ra]),np.array([tile_dec+eps]),tel_ra,tel_dec,tile_mjd,lst,hexrot_deg=0, use_hardcoded_polmis_rotmat=use_hardcoded_polmis_rotmat)
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

    if from_platemaker :
        # apply tranformation from desimeter to platemater
        dm2pm = DM2PM.read(dm2pm_filename())
        xfp,yfp = dm2pm.pm2dm(xfp,yfp)

    xtan,ytan = fp2tan(xfp,yfp,adc1,adc2)
    tmp_ra,tmp_dec    = tan2radec(xtan,ytan,tel_ra,tel_dec,tile_mjd,lst,hexrot_deg=0, use_hardcoded_polmis_rotmat=use_hardcoded_polmis_rotmat)

    # need to deal with field rotation
    # measure field rotation
    tmp_fieldrot = _measure_fieldrot_deg(-tmp_ra,tmp_dec,-tile_ra,tile_dec,xfp,yfp)

    # apply field rotation to match request
    drot = tile_fieldrot-tmp_fieldrot
    s,c  = sincosd(drot)
    tmp_xfp = c * xfp + s * yfp
    tmp_yfp = -s * xfp + c * yfp

    xtan,ytan      = fp2tan(tmp_xfp,tmp_yfp,adc1,adc2)
    ra,dec = tan2radec(xtan,ytan,tel_ra,tel_dec,tile_mjd,lst,hexrot_deg=0, use_hardcoded_polmis_rotmat=use_hardcoded_polmis_rotmat)

    # verify
    realised_fieldrot = _measure_fieldrot_deg(-ra,dec,-tile_ra,tile_dec,xfp,yfp)
    log.info("Requested fieldrot={:3.1f} arcsec delta={:3.1f} arcsec".format(tile_fieldrot*3600.,(tile_fieldrot-realised_fieldrot)*3600.))

    return ra,dec


def fiberassign_radec2xy_flat(ra,dec,tile_ra,tile_dec,tile_mjd,tile_ha,tile_fieldrot,adc1=None,adc2=None,to_platemaker=True, use_hardcoded_polmis_rotmat=False) :
    """Computes X Y focal plane coordinates of targets in the curved, "flat" coordinate system.
    Args:
      ra  : 1D numpy.array with target RA in degrees (need an array of points to compute the field rotation)
      dec : 1D numpy.array with target Dec in degrees, same size as ra
      tile_ra : float, RA of tile center in degrees, such that after pointing and ADC shift corrections, this corresponds to X=0,Y=0
      tile_dec : float, Dec of tile center in degrees, such that after pointing and ADC shift corrections, this corresponds to X=0,Y=0
      tile_mjd : float, days, used for precession (does not need to be super accurate)
      tile_ha :  float, hour angle of observation, in degrees
      tile_fieldrot :  float, design/requested field rotation, in degrees
      adc1 : optional, float, ADC1 angle in degrees: by default it is computed based on the HA and Dec.
      adc2 : optional, float, ADC2 angle in degrees: by default it is computed based on the HA and Dec.
      to_platemaker : assume output coordinates are in platemaker system and use dm2pm transform
      use_hardcoded_polmis_rotmat (optional, defaults to False): use pre-computed polar misalignment matrix, in desimeter.transform.radec2tan.get_hardcoded_polmis_rotmat() (bool)
    Returns xfp,yfp, focal plane coordinates in mm in the curved "flat" coordinate system, 1D numpy arrays of same size as ra,dec
    """

    xfp,yfp = fiberassign_radec2xy_cs5(ra,dec,tile_ra,tile_dec,tile_mjd,tile_ha,tile_fieldrot,adc1,adc2,to_platemaker=to_platemaker, use_hardcoded_polmis_rotmat=use_hardcoded_polmis_rotmat)

    # fiber assign coordinates are on the curved coordinates that follow the curved focal surface
    # the curved coordinates are called 'flat' in the focalplane parlance.
    xflat,yflat = ptl2flat(xfp,yfp)

    return xflat,yflat



def fiberassign_flat_xy2radec(xflat,yflat,tile_ra,tile_dec,tile_mjd,tile_ha,tile_fieldrot,adc1=None,adc2=None,from_platemaker=True, use_hardcoded_polmis_rotmat=False) :
    """Computes RA Dec from focal plane coordinates of targets in curved "flat" coordinate system (inverse of fiberassign_radec2xy_flat)
    Args:
      xflat  : 1D numpy.array with target RA in degrees (need an array of points to compute the field rotation)
      yflat : 1D numpy.array with target Dec in degrees, same size as ra
      tile_ra : float, RA of tile center in degrees, such that after pointing and ADC shift corrections, this corresponds to X=0,Y=0
      tile_dec : float, Dec of tile center in degrees, such that after pointing and ADC shift corrections, this corresponds to X=0,Y=0
      tile_mjd : float, days, used for precession (does not need to be super accurate)
      tile_ha :  float, hour angle of observation, in degrees
      tile_fieldrot :  float, design/requested field rotation, in degrees
      adc1 : optional, float, ADC1 angle in degrees: by default it is computed based on the HA and Dec.
      adc2 : optional, float, ADC2 angle in degrees: by default it is computed based on the HA and Dec.
      from_platemaker : assume input coordinates are from platemaker and use pm2dm transform
      use_hardcoded_polmis_rotmat (optional, defaults to False): use pre-computed polar misalignment matrix, in desimeter.transform.radec2tan.get_hardcoded_polmis_rotmat() (bool)
    Returns RA,Dec, 1D numpy arrays of same size as xflat,yflat
    """

    xfp,yfp = flat2ptl(xflat,yflat)
    ra,dec  = fiberassign_cs5_xy2radec(xfp,yfp,tile_ra,tile_dec,tile_mjd,tile_ha,tile_fieldrot,adc1=adc1,adc2=adc2,from_platemaker=from_platemaker, use_hardcoded_polmis_rotmat=use_hardcoded_polmis_rotmat)
    return ra,dec
