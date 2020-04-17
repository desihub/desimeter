import numpy as np

def xy2qs(x, y):
    '''Focal tangent plane x,y -> angular q,s on curved focal surface

    Args:
        x, y: cartesian location on focal tangent plane in mm

    Returns (q, s) where q=angle in degrees; s=focal surface radial dist [mm]

    Notes: (x,y) are in the "CS5" DESI coordinate system tangent plane to
    the curved focal surface.  q is the radial angle measured counter-clockwise
    from the x-axis; s is the radial distance along the curved focal surface;
    it is *not* sqrt(x**2 + y**2).  (q,s) are the preferred coordinates for
    the DESI focal plane hardware engineering team.
    '''
    r = np.sqrt(x**2 + y**2)
    s = r2s(r)
    q = (np.degrees(np.arctan2(y, x)) + 360.0) % 360.0
    return q, s

def qs2xy(q, s):
    '''angular q,s on curved focal surface -> focal tangent plane x,y

    Args:
        q: angle in degrees
        s: focal surface radial distance in mm

    Returns (x, y) cartesian location on focal tangent plane in mm

    Notes: (x,y) are in the "CS5" DESI coordinate system tangent plane to
    the curved focal surface.  q is the radial angle measured counter-clockwise
    from the x-axis; s is the radial distance along the curved focal surface;
    it is *not* sqrt(x**2 + y**2).  (q,s) are the preferred coordinates for
    the DESI focal plane hardware engineering team.
    '''
    r = s2r(s)
    x = r*np.cos(np.radians(q))
    y = r*np.sin(np.radians(q))
    return x, y

def r2s(r):
    '''Convert polar coordinate r to radial distance along focal surface q.'''
    # fitted on desimodel/data/focalplane/fiberpos.ecsv
    # residuals are < 0.4 microns
    # JHS -- I find this not true for 0 < r < 25
    #        For example r2s(0) == 0.0036, clearly not quite right
    #        Suggest comparison to posconstants.R2S_lookup()
    r = r if isinstance(r, np.ndarray) else np.array(r)
    c = np.array([-3.01232440e-03,  1.45324708e-02, -2.55244612e-02,  2.15885180e-02, -8.05287872e-03,  2.05529419e-03,  9.99773920e-01,  9.12275165e-06])
    pol = np.poly1d(c)    
    s = 400.*pol(r/400.)
    return s

def s2r(s):
    '''Convert radial distance along focal surface to polar coordinate r.'''
    # fitted on desimodel/data/focalplane/fiberpos.ecsv
    # residuals are < 0.4 microns
    s = s if isinstance(s, np.ndarray) else np.array(s)
    c = np.array([-2.60833797e-03,  6.40671681e-03, -5.64913181e-03,  6.99354170e-04, -2.13171265e-04,  1.00000009e+00,  9.75790364e-07])
    pol = np.poly1d(c)
    r = 400.*pol(s/400.)
    return r

#- Transform between CS5 x,y and curved focal surface
def xy2uv(x, y):
    q, s = xy2qs(x, y)
    qrad = np.radians(q)
    u = s*np.cos(qrad)
    v = s*np.sin(qrad)
    return u, v

def uv2xy(u, v):
    s = np.sqrt(u**2 + v**2)
    qrad = np.arctan2(v, u)
    q = np.degrees(qrad)
    x, y = qs2xy(q, s)
    return x, y
