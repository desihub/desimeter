import numpy as np
###################################################################################
D2R = np.pi/180.
R2D = 1./D2R
###################################################################################

def sind(degrees):
    return np.sin(D2R*degrees)

def cosd(degrees):
    return np.cos(D2R*degrees)

def sincosd(deg):
    r = np.deg2rad(deg)
    return np.sin(r), np.cos(r)

def rot2deg(deg):
    '''Returns a 2-d rotation matrix by angle 'deg'.'''
    sa,ca=sincosd(deg)
    return np.array([[ca,-sa],[sa,ca]])

def tand(degrees):
    return np.tan(D2R*degrees)

def put360(degrees):  # Puts an angle into range 0 to 360.
    return degrees % 360

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

def average_angles_deg(a1, a2):
    '''Computes the average of two angles (in degrees),
    handling wrap-around.  That is, if the two angles are thought of as
    directions in the plane, this gives the midpoint between them.'''
    x = cosd(a1) + cosd(a2)
    y = sind(a1) + sind(a2)
    return arctan2d(y, x)
