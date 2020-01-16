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
from desimeter.log import get_logger

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

###################################################################################
D2R = np.pi/180.
R2D = 1./D2R
###################################################################################

def sind(degrees):
    return np.sin(D2R*degrees)

def cosd(degrees):
    return np.cos(D2R*degrees)

def tand(degrees):
    return np.tan(D2R*degrees)

def put360(degrees): # Puts an angle into range 0 to 360.
    return np.fmod(720.+degrees, 360)

def arctan2d(y, x):
    return put360(np.arctan2(y, x)*R2D)

def arcsind(x):
    return np.arcsin(x)*R2D

def getXYZ(lon,lat):  # Convert spherical angles (in degrees) into xyz triplet
    return np.array([cosd(lon)*cosd(lat), sind(lon)*cosd(lat), sind(lat)])

def getNorm(xyz):
    return np.sqrt(np.sum(xyz**2,axis=0))
    
def getNormalized(xyz):
    return xyz/getNorm(xyz)

def getLONLAT(xyz): # Convert xyz into its spherical angles
    xyz = getNormalized(xyz)  # usually unnecessary
    return  arctan2d(xyz[1],xyz[0]) , arcsind(xyz[2])  # degrees

def vecX(xdeg):  # For positive xdeg=cwRoll: +y=>+z; +z=>-y.
    c=np.cos(np.radians(xdeg)); s=np.sin(np.radians(xdeg))
    return np.array([[1,0,0], [0,c,-s], [0,+s,c]])

def vecY(ydeg):  # For positive ydeg=-elev: +z=>+x; +x=>-z.
    # do not overlook this minus sign: positive ydeg pitches downward.
    c=np.cos(np.radians(ydeg)); s=np.sin(np.radians(ydeg))
    return np.array([[c,0,+s], [0,1,0], [-s,0,c]])
    
def vecZ(zdeg):  # For positive zdeg=+az: +x=>+y; +y=>-x.
    c=np.cos(np.radians(zdeg)); s=np.sin(np.radians(zdeg))
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

def precess(ra, dec, years) :
    # see DESI-4957
    # Equator and zero longitude point glide westward at 0.0139 deg/year, so..
    # Star ecliptic longitudes increase +0.0139 deg/year from combined lunar and solar torques on the Earth.
    # To precess any sta'’s {RA,DEC}:
    # 1. Convert to ecliptic coordinates {elon, elat}
    # 2. Add 0.0139 deg * number of years to elon
    # 3. Convert back to {RA,DEC}
    deltaELON = 360.* (years / PRECESSION_PERIOD_IN_YEARS) # degrees
    lon,lat=radec2eclip(ra,dec)
    xyz_ecliptic  = getXYZ(lon,lat)
    xyz_precessed = np.dot(vecZ(deltaELON), xyz_ecliptic) 
    lon,lat = getLONLAT(xyz_precessed)
    return eclip2radec(lon,lat)

def precess_from_icrs(ra, dec, mjd, use_astropy = False) :
    if use_astropy : #  using astropy
        from astropy.coordinates import SkyCoord,FK5,GCRS
        import astropy.time
        fk5  = FK5(equinox=astropy.time.Time(mjd,format="mjd")) # precession
        c1   = SkyCoord(ra,dec, frame='icrs', unit='deg')
        c2   = c1.transform_to(fk5)
        ra2  = c2.ra.value
        dec2 = c2.dec.value
    else :
        # precession    
        years = (mjd-ICRS_MJD)/DAYS_PER_YEAR
        ra2, dec2 = precess(ra, dec, years)
    return ra2,dec2

def apply_refraction(alt) :
    return alt + REFRACTION_AT_30DEG_ELEVATION/3600.*tand(30.)/tand(alt)  # deg , refraction per point in field

def compute_polar_misalignment_rotation_matrix() :
    # define a rotation matrix to move the polar axis to the north
    # vector product
    ha1,dec1=altaz2hadec(alt=LATITUDE-ME_ARCSEC/3600.,az=MA_ARCSEC/3600.)
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
    
    cha  = cosd(ha)
    sha  = sind(ha)
    cdec = cosd(dec)
    sdec = sind(dec)
    clat = cosd(LATITUDE)
    slat = sind(LATITUDE)    
    x = - cha * cdec * slat + sdec * clat
    y = - sha * cdec
    z = cha * cdec * clat + sdec * slat
    r = np.sqrt(x**2 + y**2)
    az  = arctan2d(y,x)
    alt = arctan2d(z,r)
    return alt,az

def altaz2hadec(alt,az) : 
    calt = cosd(alt)
    salt = sind(alt)
    caz  = cosd(az)
    saz  = sind(az)
    clat = cosd(LATITUDE)
    slat = sind(LATITUDE)
    ha  = arctan2d( -saz*calt, -caz*slat*calt+salt*clat)
    dec = arcsind(slat*salt+clat*calt*caz)
    return ha,dec


def hadec2xy(ha,dec,cha,cdec) :
    """
    dec is up (y>0)
    ra is left (x<0)
    ha is right (x>0) (ha=lst-ra)
    """
    xyz = getXYZ(ha,dec)
    ch= cosd(cha)
    sh= sind(cha)
    sd= sind(cdec)
    cd= cosd(cdec)
    rh=np.array([[ch,sh,0],[-sh,ch,0],[0,0,1]]) # rotation about HA axis
    rd=np.array([[sd,0,-cd],[0,1,0],[+cd,0,sd]]) # rotation about Dec axis
    tmp = rd.dot(rh.dot(xyz))
    x=tmp[1]
    y=-tmp[0]
    return x,y
    

def xy2hadec(x,y,cha,cdec) :
    ch= cosd(cha)
    sh= sind(cha)
    sd= sind(cdec)
    cd= cosd(cdec)
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

def radec2tan(ra,dec,tel_ra,tel_dec,mjd,lst_deg,hexrot_deg, aberration = False, polar_misalignment = False, use_astropy = False) :
    
    log = get_logger()
    
    # precession
    ra,dec = precess_from_icrs(ra, dec, mjd, use_astropy)

    
    if aberration :
        log.warning("aberration not implemented")

    # ha,dec
    ha    = lst_deg - ra
    tel_ha = lst_deg - tel_ra
    
    if polar_misalignment :
        # rotate
        polar_misalignment_matrix = compute_polar_misalignment_rotation_matrix()
        ha,dec  = getLONLAT(polar_misalignment_matrix.dot(getXYZ(ha,dec)))
    
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
    if hexrot_deg != 0 :
        log.warning("I am not sure of the hexapod rotation angle")
    
    chex = cosd(hexrot_deg)
    shex = sind(hexrot_deg)
    
    
    return chex*x+shex*y , -shex*x+chex*y
    
