"""
Transformations between tangent plane and focal plane xy based upon the
Echo22 optics model in DESI-0329v18 Echo22Platescale.txt
"""

import numpy as np
from astropy.table import Table
from pkg_resources import resource_filename

_r2t_coeff = None
_t2r_coeff = None
_rscale = 415.0

def tan2fp(xtan, ytan):
    """
    Convert tangent plane coordinates to focal plane mm
    
    Args:
        xtan: sin(theta)*cos(phi)
        ytan: sin(theta)*sin(phi)

    Returns xfp, yfp in CS5 coordinates on the focal plane
    """
    phi = np.arctan2(ytan, xtan)
    rtan = np.sqrt(xtan**2 + ytan**2)
    theta = np.degrees(np.arcsin(rtan))
    r = theta2radius(theta)
    x = r*np.cos(phi)          
    y = r*np.sin(phi)
    return x, y

def fp2tan(xfp, yfp):
    """
    Convert focal plane to tangent plane coordinates
    
    Args:
        xfp, yfp: CS5 focal plane coordinates in mm

    return xtan, ytan where xtan=sin(theta)*cos(phi), ytan=sin(theta)*sin(phi)
    """
    #- phi=0 aligned with +xtan = -RA = +HA = +xfp
    phi = np.arctan2(yfp, xfp)
    r = np.sqrt(xfp**2 + yfp**2)
    theta = radius2theta(r)
    
    rtan = np.sin(np.radians(theta))
    
    xtan = rtan * np.cos(phi)
    ytan = rtan * np.sin(phi)

    return xtan, ytan

def radius2theta(radius):
    """
    Convert radius on focal plane [mm] to radial angle [deg]
    """
    global _r2t_coeff, _t2r_coeff, _rscale
    
    if _r2t_coeff is None:
        _r2t_coeff, _t2r_coeff = _fit_r2theta()

    return np.polyval(_r2t_coeff, radius/_rscale)

def theta2radius(theta):
    """
    Convert radial angle theta [deg] to radius on the focal plane [mm]
    """
    global _r2t_coeff, _t2r_coeff, _rscale
    
    if _t2r_coeff is None:
        _r2t_coeff, _t2r_coeff = _fit_r2theta()

    return np.polyval(_t2r_coeff, theta)

def _fit_r2theta(deg=5):
    """
    Fit radius vs. radial angle using Echo22 design
    
    Args:
        deg (int): degree of polynomials to fit
    
    Returns r2t_coeff, t2r_coeff: coefficients for np.polyval for transforming
        radius -> theta and theta -> radius.
    
    Note: use radius/_rscale to normalize radius.
    """
    global _rscale

    #- Format of Echo22Platescale.txt file:
    """
    # Radius	Theta	MF/#	SF/#	MFL	SFL	MPS	SPS
    0	0	3.678707383	3.678707457	13.92433205	13.92433205	67.50706676	67.50706676
    1	0.004114797	3.678711788	3.678708901	13.92436237	13.92435007	67.50721377	67.50715416
    ...
    """
    echo22file = resource_filename('desimeter.transform.tan2fp', 'data/Echo22Platescale.txt')
    dtype = [('radius', float), ('theta', float)]
    echo22 = np.loadtxt(echo22file, dtype=dtype, usecols=[0,1])
    
    #- Maximum fiber reach is 413.5 mm, so don't try to fit well beyond that
    echo22 = echo22[echo22['radius'] <= _rscale]

    r2t_coeff = np.polyfit(echo22['radius']/_rscale, echo22['theta'], deg)
    t2r_coeff = np.polyfit(echo22['theta'], echo22['radius'], deg)
    
    return r2t_coeff, t2r_coeff
    
    # print(' n  rms_dt max_dt rms_dr max_dr')
    # for n in range(3,15):
    #     r2t_coeff = np.polyfit(echo22['radius']/_rscale, echo22['theta'], n)
    #     t2r_coeff = np.polyfit(echo22['theta'], echo22['radius'], n)
    #     theta = np.polyval(r2t_coeff, echo22['radius']/_rscale)
    #     radius = np.polyval(t2r_coeff, echo22['theta'])
    #
    #     dtheta = 3600*(theta - echo22['theta'])     #- arcsec
    #     dradius = 1000*(radius - echo22['radius'])  #- microns
    #     dtheta_rms = np.sqrt(np.mean(dtheta**2))
    #     dradius_rms = np.sqrt(np.mean(dradius**2))
    #
    #     print('{:2d}  {:.3f}  {:.3f}  {:.3f}  {:.3f}'.format(
    #         n, dtheta_rms, np.max(np.abs(dtheta)),
    #         dradius_rms, np.max(np.abs(dradius))))
    
