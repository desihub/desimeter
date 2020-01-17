
def mjd2lst(mjd) :
    
    """
    Converts a floating point MJD in days to local sideral time in degrees for the Mayall telescope

    All we care about here is consistency with what is done at the telescope. So, this is less accurate 
    than using values in the headers provided by the telescope, (like LST = HA + RA), but, unfortunatly, 
    we don't always have this.
    """

    
    # Exposure /global/homes/j/jguy/data/20191101/00023498/gfa-00023498.fits.fz
    # TCSST   = '22:32:04.726'       / Local Sidereal time reported by TCS (HH:MM:SS) 
    # TCSMJD  =         58789.135177 / MJD reported by TCS
    # (22+32/60.+4.726/3600.)*360./24. = 338.01969167
    # 
    # A few month later,
    # Exposure /global/homes/j/jguy/data/20200116/00041433/gfa-00041433.fits.fz
    # MJD-OBS =       58865.54008097 / Modified Julian Date of observation            
    # ST      = '13:17:13.130'       / Local Sidereal time at observation start (HH:MM
    # TCSST   = '13:16:57.957'       / Local Sidereal time reported by TCS (HH:MM:SS) 
    # TCSMJD  =         58865.540489 / MJD reported by TCS                            
    # giving (13+16/60.+57.957/3600.)*360./24. = 199.24148750
    #
    # I measure a drift of 2.8 arcsec over 2.5 month ...
    #
    # Base on routine ct2lst.pro in IDL astrolib.
    # But with a correction term for Mayall based on exposure 23498.
    
    
    LONGITUDE = -111.59989 # Mayall 
    jd = mjd + 2400000.5
    c = [280.46061837, 360.98564736629, 0.000387933, 38710000.0 ]
    jd2000 = 2451545.0
    t0 = jd - jd2000
    t = t0/36525.
    theta = c[0] + (c[1] * t0) + t**2*(c[2] - t/ c[3] )
    # correction term calibrated on exposure 23498 = -0.16035388
    lst = (theta + LONGITUDE - 0.16035388) % 360.
    return lst
