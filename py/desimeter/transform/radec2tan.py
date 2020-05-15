"""desimeter.transform.radec2tan
=============================

Routines to transform RA,Dec coordinates in ICRS to tangent plane coordinates.

Tangent plane coordinates give the direction of incoming photons on
the telescope with respect to the direction defined by the combination
of the telescope pointing and hexapod rotation.

Tangent plane coordinates are :

x_tan = sin(theta)*cos(phi)
y_tan = sin(theta)*sin(phi)

where theta,phi are polar coordinates. theta=0 for along the telescope
pointing. phi=0 for a star with the same Dec as the telescope
pointing but a larger HA (or smaller RA).


This transformation includes :
- precession
- aberration
- atmospheric refraction
- polar misalignement


Inputs are pointing TEL_RA,TEL_DEC,MJD,LST,HEXROT and a list of  TARGET_RA,TARGET_DEC

"""

import numpy as np
#from desimeter.log import get_logger
from desimeter.trig import (sind, tand, put360, arctan2d, arcsind, getXYZ,
                            getNormalized, sincosd)

###################################################################################
# constants
# from http://asa.usno.navy.mil/static/files/2018/Astronomical_Constants_2018.pdf
OBLIQ                      = 23.439279444444445 # 23+26/60.+21.406/3600. , obliquity of the ecliptic, Initial  Values  at  J2000·0
DAYS_PER_YEAR              = 365.256363004
PRECESSION_PERIOD_IN_YEARS = 25771.5753382206 # 360./(5028.796195/100./3600.) , Rates of precession at J2000·0 (IAU 2006) , General precession in longitude
ICRS_MJD                   = 51544.5 # 2451545.0-2400000.5 J2000
LATITUDE                   = 31.96403 # DEG
REFRACTION_AT_30DEG_ELEVATION = 79. # arcsec

# polar mis-alignment
# values based on DESI-5186, fit of commissioning instrument field rotation rotation
MA_ARCSEC = -2. # arcsec; in the northern hemisphere, positive MA means that the pole of the mounting is to the right of due north
ME_ARCSEC = -120. # arcsec; in the northern hemisphere, positive ME means that the pole of the mounting is below the true (unrefracted) north (so elevation=latitude-me)

def getLONLAT(xyz):
    """Convert xyz into its spherical angles"""
    xyz = getNormalized(xyz)  # usually unnecessary
    return  arctan2d(xyz[1],xyz[0]) , arcsind(xyz[2])  # degrees

def vecX(xdeg):  # For positive xdeg=cwRoll: +y=>+z; +z=>-y.
    s,c = sincosd(xdeg)
    return np.array([[1,0,0], [0,c,-s], [0,+s,c]])

def vecY(ydeg):  # For positive ydeg=-elev: +z=>+x; +x=>-z.
    # do not overlook this minus sign: positive ydeg pitches downward.
    s,c = sincosd(ydeg)
    return np.array([[c,0,+s], [0,1,0], [-s,0,c]])

def vecZ(zdeg):  # For positive zdeg=+az: +x=>+y; +y=>-x.
    s,c = sincosd(zdeg)
    return np.array([[c,-s,0], [+s,c,0], [0,0,1]])

def refX(xdeg):  # Rolls reference frame clockwise about X
    return vecX(-xdeg)

def refY(elev):  # Elevates reference frame about Y
    return vecY(+elev)

def refZ(azim):  # Rotates reference frame to new azimuth
    return vecZ(-azim)

def radec2eclip(ra,dec):  # same epoch
    equatorial_xyz = getXYZ(ra,dec)
    ecliptic_xyz = np.dot(refX(OBLIQ), equatorial_xyz)  # roll frame clockwise
    return getLONLAT(ecliptic_xyz)

def eclip2radec(lon,lat):  # same epoch
    ecliptic_xyz = getXYZ(lon,lat)
    equatorial_xyz = np.dot(refX(-OBLIQ), ecliptic_xyz)  # roll frame counterclockwise by obliq
    return getLONLAT(equatorial_xyz)

def apply_precession(ra, dec, years) :
    """
    see DESI-4957
    Equator and zero longitude point glide westward at 0.0139 deg/year, so..
    Star ecliptic longitudes increase +0.0139 deg/year from combined lunar and solar torques on the Earth.
    To precess any sta'’s {RA,DEC}:
    1. Convert to ecliptic coordinates {elon, elat}
    2. Add 0.0139 deg * number of years to elon
    3. Convert back to {RA,DEC}
    """
    deltaELON = 360.* (years / PRECESSION_PERIOD_IN_YEARS) # degrees
    lon,lat=radec2eclip(ra,dec)
    xyz_ecliptic  = getXYZ(lon,lat)
    xyz_precessed = np.dot(vecZ(deltaELON), xyz_ecliptic)
    lon,lat = getLONLAT(xyz_precessed)
    return eclip2radec(lon,lat)

def apply_precession_from_icrs(ra, dec, mjd, use_astropy = False) :
    """ apply precession from ICRS to the date of observation given in MJD

    Args:
        ra: float or 1D np.array RA in degrees
        dec: float or 1D np.array dec in degrees

    Returns: alt, az
        alt: float or 1D np.array refracted altitude in degrees

    """
    if use_astropy : #  using astropy
        from astropy.coordinates import SkyCoord,FK5
        import astropy.time
        fk5  = FK5(equinox=astropy.time.Time(mjd,format="mjd")) # precession
        c1   = SkyCoord(ra,dec, frame='icrs', unit='deg')
        c2   = c1.transform_to(fk5)
        ra2  = c2.ra.value
        dec2 = c2.dec.value
    else :
        # precession
        years = (mjd-ICRS_MJD)/DAYS_PER_YEAR
        ra2, dec2 = apply_precession(ra, dec, years)
    return ra2,dec2

def undo_precession_from_icrs(ra, dec, mjd, use_astropy = False) :
    """ apply precession from ICRS to the date of observation given in MJD

    Args:
        ra: float or 1D np.array RA in degrees
        dec: float or 1D np.array dec in degrees
        mjd: observation date
        use_astropy: obvious

    Returns: alt, az
        alt: float or 1D np.array refracted altitude in degrees

    """
    if use_astropy : #  using astropy
        from astropy.coordinates import SkyCoord,FK5
        import astropy.time
        fk5  = FK5(equinox=astropy.time.Time(mjd,format="mjd")) # precession
        #c1   = SkyCoord(ra,dec, frame='icrs', unit='deg')
        #c2   = c1.transform_to(fk5)
        #ra2  = c2.ra.value
        #dec2 = c2.dec.value

        c1   = SkyCoord(ra,dec, frame=fk5, unit='deg')
        c2   = c1.transform_to('icrs')
        ra2  = c2.ra.value
        dec2 = c2.dec.value

    else :
        # precession
        years = (ICRS_MJD-mjd)/DAYS_PER_YEAR
        ra2, dec2 = apply_precession(ra, dec, years)
    return ra2,dec2

def getSunLon(mjd):  # Given JD.fff, returns sun ecliptic longitude degrees
    # https://en.wikipedia.org/wiki/Position_of_the_Sun
    #days = jd - 2451545.0  # days from Greenwich noon 2000 Jan 01
    days = mjd - 51544.5 # (MJD = JD - 2400000.5)
    mean = 280.460 + 0.9856474*days
    anom = 357.528 + 0.9856003*days
    elon = mean + 1.915*sind(anom) + 0.020*sind(2*anom)
    return put360(elon);

def compute_aberration(ra_eclip,dec_eclip, mjd, magnif):
    # Given target ecliptic lonlat, and JD,
    # returns the fully aberrated & magnified LonLat.
    speed     = 0.99365E-4 # Meeus p.151 in radians
    apexlon   = getSunLon(mjd) - 90.
    if isinstance(ra_eclip, float):
        apexXYZ   = getXYZ(apexlon,0.)
    else :
        apexXYZ   = getXYZ(apexlon*np.ones(ra_eclip.shape), np.zeros(ra_eclip.shape))
    targetXYZ = getXYZ(ra_eclip,dec_eclip)

    if len(targetXYZ.shape) > 1 :
        VxT   = np.zeros(targetXYZ.shape)
        TxVxT = np.zeros(targetXYZ.shape)
        for i in range(targetXYZ.shape[1]) :
            VxT[:,i]=np.cross(apexXYZ[:,i], targetXYZ[:,i])
            TxVxT[:,i]=np.cross(targetXYZ[:,i],VxT[:,i])
    else :
       VxT       = np.cross(apexXYZ, targetXYZ)
       TxVxT     = np.cross(targetXYZ, VxT)

    aberXYZ   = speed * TxVxT
    plotXYZ   = targetXYZ + magnif * aberXYZ
    # DANGER with getNormalized(): use only one triplet per call.
    return getLONLAT(getNormalized(plotXYZ))

def apply_aberration(ra,dec,mjd, use_astropy = False) :
    """ apply aberration given date in MJD

    Args:
        ra: float or 1D np.array RA in degrees
        dec: float or 1D np.array dec in degrees
        mjd: float, mjd

    Returns: alt, az
        alt: float or 1D np.array refracted altitude in degrees

    """
    if use_astropy :
        from astropy.coordinates import SkyCoord,GCRS #FK5
        import astropy.time

        gcrs=GCRS(obstime=astropy.time.Time(mjd,format="mjd")) # precession? + aberration
        #fk5_J2000=FK5(equinox="J2000")
        #fk5=FK5(equinox=astropy.time.Time(mjd,format="mjd")) # precession
        c1 = SkyCoord(ra,dec, frame='icrs', unit='deg')
        c2 = c1.transform_to(gcrs)
        ra_equ  = c2.ra.value
        dec_equ =  c2.dec.value
    else :
        ra_eclip,dec_eclip = radec2eclip(ra,dec)
        ra_eclip,dec_eclip = compute_aberration(ra_eclip,dec_eclip, mjd, 1.)
        ra_equ,dec_equ     = eclip2radec(ra_eclip,dec_eclip)

    return ra_equ,dec_equ

def undo_aberration(ra,dec,mjd, use_astropy = False) :
    """ undo aberration given date in MJD

    Args:
        ra: float or 1D np.array RA in degrees
        dec: float or 1D np.array dec in degrees
        mjd: float, mjd

    Returns: alt, az
        alt: float or 1D np.array refracted altitude in degrees

    """
    # at first order
    ra2,dec2 = apply_aberration(ra,dec,mjd,use_astropy)
    return 2*ra-ra2 , 2*dec-dec2

def apply_refraction(alt) :
    """ apply refraction

    Args:
        alt: float or 1D np.array altitude in degrees outside of atmosphere

    Returns: alt, az
        alt: float or 1D np.array refracted altitude in degrees

    """
    return alt + REFRACTION_AT_30DEG_ELEVATION/3600.*tand(30.)/tand(alt)  # deg , refraction per point in field

def undo_refraction(alt) :
    """ undo refraction

    Args:
        alt: float or 1D np.array altitude in degrees outside of atmosphere

    Returns: alt, az
        alt: float or 1D np.array refracted altitude in degrees

    """
    # need second order inverse
    daprev = 0.
    for _ in range(2) :
        da   = apply_refraction(alt) - alt
        dda  = da-daprev
        alt -= dda
        daprev = da
    return alt

def compute_polar_misalignment_rotation_matrix(me_arcsec,ma_arcsec) :
    """
    compute a rotation matrix to move the polar axis to the north
    vector product
    """
    ha1,dec1=altaz2hadec(alt=LATITUDE-me_arcsec/3600.,az=ma_arcsec/3600.)
    xyz1=getXYZ(ha1,dec1)
    ha2,dec2=altaz2hadec(alt=LATITUDE,az=0.)
    xyz2=getXYZ(ha2,dec2)
    cross = np.cross(xyz1,xyz2)

    norme = np.sqrt(np.sum(cross**2))
    if norme > 0 :
        polar_misalignment_rotation = rotation_matrix(cross/norme,np.arcsin(norme))
    else :
        polar_misalignment_rotation = np.eye(3)

    return polar_misalignment_rotation

def hadec2altaz(ha,dec) :
    """
     Convert HA,Dec to Altitude , Azimuth at Kitt Peak elevation

    Args:
        ha: float or 1D np.array hour angle in degrees
        dec: float or 1D np.array declination in degrees

    Returns: alt, az
        alt: float or 1D np.array altitude in degrees
        az: float or 1D np.array azimuth in degrees
    """
    sha,cha   = sincosd(ha)
    sdec,cdec = sincosd(dec)
    slat,clat = sincosd(LATITUDE)
    x = - cha * cdec * slat + sdec * clat
    y = - sha * cdec
    z = cha * cdec * clat + sdec * slat
    r = np.hypot(x, y)
    az  = arctan2d(y,x)
    alt = arctan2d(z,r)
    return alt,az

def altaz2hadec(alt,az) :
    """
     Convert Altitude, Azimuth to HA, Dec at Kitt Peak elevation

    Args:
        alt: float or 1D np.array altitude in degrees
        az: float or 1D np.array azimuth in degrees

    Returns: ha, dec
        ha: float or 1D np.array Hour Angle in degrees
        dec: float or 1D np.array Hour Angle in degrees
    """
    salt,calt = sincosd(alt)
    saz,caz   = sincosd(az)
    slat,clat = sincosd(LATITUDE)
    ha  = arctan2d( -saz*calt, -caz*slat*calt+salt*clat)
    dec = arcsind(slat*salt+clat*calt*caz)
    return ha,dec


def hadec2xy(ha,dec,tel_ha,tel_dec) :
    """
     Convert HA,Dec to tangent plane x,y given a telescope pointing tel_ha, tel_dec

    Args:
        ha: float or 1D np.array Hour Angle in degrees
        dec: float or 1D np.array Hour Angle in degrees (same size as ha)
        tel_ha: float; telescope pointing Hour Angle in degrees
        tel_dec: float; telescope pointing Dec in degrees

    Returns: x y float or 1D np.array
    """
    xyz = getXYZ(ha,dec)
    sh,ch = sincosd(tel_ha)
    sd,cd = sincosd(tel_dec)
    rh=np.array([[ch,sh,0],[-sh,ch,0],[0,0,1]]) # rotation about HA axis
    rd=np.array([[sd,0,-cd],[0,1,0],[+cd,0,sd]]) # rotation about Dec axis
    tmp = rd.dot(rh.dot(xyz))
    x=tmp[1]
    y=-tmp[0]
    return x,y


def xy2hadec(x,y,tel_ha,tel_dec) :
    """
    Convert tangent plane x,y to HA,Dec given a telescope pointing tel_ha, tel_dec

    Args:
        x: float or 1D np.array with tangent plane coord.
        y: float or 1D np.array with tangent plane coord.
        tel_ha: float; telescope pointing Hour Angle in degrees
        tel_dec: float; telescope pointing Dec in degrees

    Returns: HA, Dec  float or 1D np.array in degrees
    """
    sh,ch = sincosd(tel_ha)
    sd,cd = sincosd(tel_dec)
    rh=np.array([[ch,-sh,0],[sh,ch,0],[0,0,1]])
    rd=np.array([[sd,0,cd],[0,1,0],[-cd,0,sd]])
    z = np.sqrt(1-x**2-y**2)
    xyz = rh.dot(rd.dot(np.array([-y,x,z])))
    return getLONLAT(xyz)

def rotation_matrix(axis,theta) :
    # rotate around axis by an angle
    ct=np.cos(theta)
    st=np.sin(theta)
    u=axis
    uut=np.outer(u,u)
    uxx=np.array([[0,-u[2],u[1]],[u[2],0,-u[0]],[-u[1],u[0],0]])
    return ct*np.eye(axis.size)+st*uxx+(1-ct)*uut

def radec2tan(ra,dec,tel_ra,tel_dec,mjd,lst_deg,hexrot_deg, precession = True, aberration = True, polar_misalignment = True, use_astropy = False) :
    """
    Convert ICRS coordinates to tangent plane coordinates

    Args:
        ra: float or 1D np.array with RA in degrees
        dec: float or 1D np.array with RA in degrees
        tel_ra: float, in degrees, telescope pointing RA
        tel_dec: float, in degrees, telescope pointing Dec
        mjd: float, Modified Julian Date of observation, in days
        lst_deg: float, local sidereal time, in degrees
        hexrot_deg: float, hexapod rotation angle, in degrees

    Returns x_tan,y_tan , tangent plane coordinates:
        x_tan = sin(theta)*cos(phi) : float on np.array (same shape as input ra,dec)
        y_tan = sin(theta)*sin(phi) : float on np.array (same shape as input ra,dec)

        where theta,phi are polar coordinates. theta=0 for along the telescope
        pointing. phi=0 for a star with the same Dec as the telescope
        pointing but a larger HA (or smaller RA).

    Optional arguments:
        aberration: boolean; compute aberration if True
        polar_misalignment: boolean; compute polar misalignment if True
        use_astropy: boolean; use astropy coordinates for precession and aberration if True
    """
    if precession :
        ra,dec = apply_precession_from_icrs(ra, dec, mjd, use_astropy)
        tel_ra,tel_dec = apply_precession_from_icrs(tel_ra, tel_dec, mjd, use_astropy)


    if aberration :
        ra,dec = apply_aberration(ra,dec,mjd, use_astropy)
        tel_ra,tel_dec = apply_aberration(tel_ra,tel_dec,mjd, use_astropy)

    # ha,dec
    ha     = lst_deg - ra
    tel_ha = lst_deg - tel_ra

    if polar_misalignment :
        # rotate
        polar_misalignment_matrix = compute_polar_misalignment_rotation_matrix(me_arcsec=ME_ARCSEC,ma_arcsec=MA_ARCSEC)
        ha,dec = getLONLAT(polar_misalignment_matrix.dot(getXYZ(ha,dec)))
        tel_ha,tel_dec = getLONLAT(polar_misalignment_matrix.dot(getXYZ(tel_ha,tel_dec)))

    # alt,az
    alt,az = hadec2altaz(ha,dec)
    tel_alt,tel_az = hadec2altaz(tel_ha,tel_dec)

    # apply refraction
    alt = apply_refraction(alt)
    tel_alt = apply_refraction(tel_alt)

    # convert back to ha,dec
    ha,dec  = altaz2hadec(alt,az)
    tel_ha,tel_dec  = altaz2hadec(tel_alt,tel_az)

    # tangent plane
    x,y = hadec2xy(ha,dec,tel_ha,tel_dec)

    # hexapod rotation
    shex,chex = sincosd(hexrot_deg)

    return chex*x-shex*y , +shex*x+chex*y

def tan2radec(x_tan,y_tan,tel_ra,tel_dec,mjd,lst_deg,hexrot_deg, precession = True, aberration = True, polar_misalignment = True, use_astropy = False) :
    """
    Convert ICRS coordinates to tangent plane coordinates

    Args:
        xtan: float or 1D np.array with tangent plane coordinates
        ytan: float or 1D np.array with tangent plane coordinates
        tel_ra: float, in degrees, telescope pointing RA
        tel_dec: float, in degrees, telescope pointing Dec
        mjd: float, Modified Julian Date of observation, in days
        lst_deg: float, local sidereal time, in degrees
        hexrot_deg: float, hexapod rotation angle, in degrees

    Returns RA,Dec in ICRS system (hopefully)

    Optional arguments:
        aberration: boolean; compute aberration if True
        polar_misalignment: boolean; compute polar misalignment if True
        use_astropy: boolean; use astropy coordinates for precession and aberration if True
    """

    # undo hexapod rotation
    shex,chex = sincosd(hexrot_deg)
    x =  chex*x_tan+shex*y_tan
    y = -shex*x_tan+chex*y_tan

    # need to apply precession ... etc to telescope pointing to interpret the x,y
    if precession :
        tel_ra,tel_dec = apply_precession_from_icrs(tel_ra, tel_dec, mjd, use_astropy)
    if aberration :
        tel_ra,tel_dec = apply_aberration(tel_ra,tel_dec,mjd, use_astropy)

    tel_ha = lst_deg - tel_ra

    if polar_misalignment :
        polar_misalignment_matrix = compute_polar_misalignment_rotation_matrix(me_arcsec=ME_ARCSEC,ma_arcsec=MA_ARCSEC)
        tel_ha,tel_dec = getLONLAT(polar_misalignment_matrix.dot(getXYZ(tel_ha,tel_dec)))

    # we need to apply refraction for the telescope pointing to interpret the x,y
    tel_alt,tel_az = hadec2altaz(tel_ha,tel_dec)
    # apply refraction
    refracted_tel_alt = apply_refraction(tel_alt)
    # back to ha,dec
    refracted_tel_ha,refracted_tel_dec  = altaz2hadec(refracted_tel_alt,tel_az)
    # now convert x,y to ha,dec
    refracted_ha,refracted_dec = xy2hadec(x,y,refracted_tel_ha,refracted_tel_dec)

    # alt,az
    alt,az = hadec2altaz(refracted_ha,refracted_dec)

    # undo refraction
    alt = undo_refraction(alt)

    # back to ha,dec
    ha,dec  = altaz2hadec(alt,az)

    # now polar mis-alignment
    if polar_misalignment :
        # inverse matrix
        polar_misalignment_matrix = compute_polar_misalignment_rotation_matrix(me_arcsec=-ME_ARCSEC,ma_arcsec=-MA_ARCSEC)
        ha,dec  = getLONLAT(polar_misalignment_matrix.dot(getXYZ(ha,dec)))

    # ra
    ra = lst_deg - ha

    if aberration :
        ra,dec = undo_aberration(ra,dec,mjd, use_astropy)

    if precession :
        ra,dec = undo_precession_from_icrs(ra, dec, mjd, use_astropy)

    return ra,dec
