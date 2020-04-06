"""
Utility functions to fit circles from a set of cartesian coordinates
"""

import numpy as np
from scipy import optimize

#- Least squares fit to circle; adapted from "Method 2b" in
#- https://scipy-cookbook.readthedocs.io/items/Least_Squares_Circle.html

def fit_circle(x, y):
    """
    fit a circle from a set of 2D cartesian coordinates
    Args:
        x : float numpy array of coordinates along first axis of cartesian coordinate system
        y : float numpy array of coordinates along second axis in same system

    returns:
        xc : float, coordinates along first axis of center of circle 
        yc : float, coordinates along second axis of center of circle
        r  : float, radius of circle
    """
    
    x_m, y_m, r = _fast_fit_circle(x, y)
    
    #- If r is too small or too big then either this positioner wasn't moving
    #- or the points are mismatched for a bad fit.
    if (r < 1.0) or (r > 5.0):
        raise ValueError('Bad circle fit')

    def calc_R(xc, yc):
        """ calculate the distance of each data points from the center (xc, yc) """
        return np.sqrt((x-xc)**2 + (y-yc)**2)

    def f_2b(c):
        """ calculate the algebraic distance between the 2D points and the mean circle centered at c=(xc, yc) """
        Ri = calc_R(*c)
        return Ri - Ri.mean()

    def Df_2b(c):
        """ Jacobian of f_2b
        The axis corresponding to derivatives must be coherent with the col_deriv option of leastsq"""
        xc, yc     = c
        df2b_dc    = np.empty((len(c), x.size))

        Ri = calc_R(xc, yc)
        df2b_dc[0] = (xc - x)/Ri                   # dR/dxc
        df2b_dc[1] = (yc - y)/Ri                   # dR/dyc
        df2b_dc    = df2b_dc - df2b_dc.mean(axis=1)[:, np.newaxis]

        return df2b_dc

    
    center_estimate = x_m, y_m
    center_2b, ier = optimize.leastsq(f_2b, center_estimate, Dfun=Df_2b, col_deriv=True)

    xc_2b, yc_2b = center_2b
    Ri_2b        = calc_R(*center_2b)
    R_2b         = Ri_2b.mean()
    residu_2b    = sum((Ri_2b - R_2b)**2)

    if (R_2b < 1.0) or (R_2b > 5.0):
        raise ValueError('Bad circle fit')

    return xc_2b, yc_2b, R_2b

def _fast_fit_circle(x,y) :
    """
    fast and approximate method to fit a circle from a set of 2D cartesian coordinates
    (better use fit_circle)
    Args:
        xin : float numpy array of coordinates along first axis of cartesian coordinate system
        yin : float numpy array of coordinates along second axis in same system

    returns:
        xc : float, coordinates along first axis of center of circle 
        yc : float, coordinates along second axis of center of circle
        r  : float, radius of circle
    """
    # init
    nn=len(x)
    i1=np.arange(nn)
    ### i2=(i1+1)%nn
    i2=(i1+nn//2-1)%nn
    
    # midpoints
    mx=((x[i1]+x[i2])/2.)
    my=((y[i1]+y[i2])/2.)
    nx=(y[i2]-y[i1])
    ny=-(x[i2]-x[i1])

    # solve for intersection of perpendicular bisectors
    # with s1,s2 are affine parameters of 2 adjacent perpendicular bisectors
    # 2 equations:
    # mx1 + nx1*s1 = mx2 + nx2*s2
    # my1 + ny1*s1 = my2 + ny2*s2
    s1 = (ny[i2]*mx[i2]-nx[i2]*my[i2]-ny[i2]*mx[i1]+nx[i2]*my[i1])/(ny[i2]*nx[i1]-nx[i2]*ny[i1])

    # coordinates of intersections are estimates of center of circle
    xc=mx[i1]+nx[i1]*s1[i1]
    yc=my[i1]+ny[i1]*s1[i1]
    
    # first estimate of center is mean of all intersections
    xc=np.mean(xc)
    yc=np.mean(yc)
    r=np.mean(np.sqrt((x-xc)**2+(y-yc)**2))
   
    return xc,yc,r

