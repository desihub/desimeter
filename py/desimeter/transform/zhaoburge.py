# code extracted from DESI-3798 IterFitter8.py


#-------------------IMPORT ZONE-----------------------

import numpy as np
import time

#-----------ZERNIKE NORMALIZED EDGE=1 PACKAGE-----------
#-----------Using BORN-WOLF {n,m} definitions-----------

def factorial(n):
    if n > 1:
       return int(n*factorial(n-1))
    else:
       return 1
       
def convertNolltoBW(noll):
    # converts a Noll Zernike index to the B&W {n,m} doublet
    n = int(-0.5+np.sqrt(2*noll-1.75))
    m = 0
    base = int(0.5*n**2 +0.5*n + 1)
    diff = noll - base
    if n%2==0:
        m = 2*int(0.5*diff + 0.7)
    else:
        m = 2*int(0.5*diff + 1.2) - 1
    if noll%2>0:  
        m = -m
    return np.array([n, m])

def convertWyanttoBW(wyant):
    # converts a Wyant Zernike index to the B&W {n,m}   
    halfsum = int(np.sqrt(wyant))
    idif = int(wyant - halfsum**2)
    halfdif = int (idif/2)
    n = halfsum + halfdif
    m = halfsum - halfdif
    if idif%2 >0:
        m = -m
    return np.array([n, m]) 
    
def getZernFuncXY(nm, xnorm, ynorm):   # nm is the BornWolf index; BIG NINE:  #1
    # Here, xnorm and ynorm must lie within the unit circle
    rnorm = np.sqrt(xnorm*xnorm + ynorm*ynorm)
    angle = np.arctan2(ynorm,xnorm)
    return getZernRadial(nm,rnorm) * getZernAngular(nm,angle)

def getZernRadial(nm, rnorm):    # BIG NINE: #4
    n = nm[0]             # B&W
    m = np.abs(nm[1])     # B&W
    halfsum = (n+m)/2
    idif = n-m
    halfdif = int(idif/2)
    # n = halfsum + halfdif   # or, halfsum = (n+m)/2
    # m = halfsum - halfdif   # or, halfdif = (n-m)/2
    # loop through the polynomial
    result = 0.
    for i in range(0, halfdif+1):
        expon = int(n-2*i)
        sign = 1 if i%2 == 0 else -1
        numer = sign * factorial(n-i)
        denom = factorial(i) * factorial(halfsum-i) * factorial(halfdif-i)
        coef = numer / denom
        term = coef*np.power(rnorm, expon)
        result = result + term
    return result  
    
def getZernAngular(nm, theta): 
    m = nm[1]    # B&W
    if m==0:
        return 1.
    if m>0:
        return np.cos(m*theta)
    m = np.abs(m)         # note this abs() function
    return np.sin(m*theta)  

def zernFormulaText(nm):           # BIG NINE: #8
    #---generates a text representation of the specified Zernike function--
    n = nm[0]  # B&W
    m = nm[1]  # B&W
    # print 'New zernFormulaText() is using n, m: ', n, m
    needsine = True if m<0 else False
    m = np.abs(m)
    halfsum = (n+m)/2
    idif = n-m
    halfdif = int(idif/2)
    
    #--first do the radial part----
    nterms = 0
    s = ''
    #--evaluate the radial polynomial-----
    for i in range(0, halfdif+1):
        nterms = nterms + 1
        # print "Starting with n,  m, i, nterms = ", n, m, i, nterms
        expon = int(n-2*i)                # start with highest exponent
        # print "  expon = ", expon
        sign = 1 if i%2 == 0 else -1      # alternating signs in Zernike series
        strsign = '+' if i%2==0 else '-'
        numer = sign * factorial(n-i)
        denom = factorial(i) * factorial(halfsum-i) * factorial(halfdif-i)
        coef = numer / denom
        scoef = str(coef)
        if coef==1 and expon>0:          # suppress showing coef=1
            scoef = ''
        if coef > 0 and nterms > 1:
            scoef = '+' + scoef
        s = s + scoef
        if expon > 0:
            s = s + 'r'
        if expon > 1:
            s = s + '^' + str(expon)
    if nterms>1 and m!=0:
        s = '('+ s + ')'
        
    #--then do the azimuthal part, if any--------
    if m==0:
        return s
    strm = ''
    if m>1:    
        strm = str(m)
    if needsine:
        s = s + '*sin(' + strm + 't)'
    else:
        s = s + '*cos(' + strm + 't)'
    return s

#----------End Zernike package with {n,m} B&W indexing---------------------


#-----Zhao-Burge functions built from Zernikes Package---------------------

rh = np.sqrt(0.5)
rt = np.sqrt(2.0)

NCOEFS = 33  # available item coefs

#LUT = [2,  5,  6,   9,  20,  28, 29,  30]    # 8 polynomials

#--------Zernike-Noll Renormalization for AreaIntegral = Pi----------------

squares = np.array([0, 1, 4, 4, 3, 6, 6, 8, 8, 8, 8,  5, 10, 10, 10, 10, 12, 12, 12, 12, 12, 12, 7, 14, 14, 14, 14, 14, 14, 16, 16, 16, 16, 16, 16, 16, 16, 9])
# for Noll index =  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,12, 13, 14, 15, 16, 17, 18, 19, 20, 21,22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37 

def normalizeArea(noll):
    # converts an edge=1 Zernike to an area=Pi normalized Zernike
    return np.sqrt(squares[noll])
  
def getZ(noll, x, y):
    return normalizeArea(noll) * getZernFuncXY(convertNolltoBW(noll), x, y)

def getZhaoBurgeXY(polids, coeffs, x, y):
    """
    Args:
        coeffs: array of coefficients
        x,y : locations at which to evaluate

    returns dx, dy arrays
    """
    dx = np.zeros(len(x))
    dy = np.zeros(len(y))
    for polid, coeff in zip(polids,coeffs):
        zbx, zby, name = getZhaoBurgeTerm(polid, x, y)
        dx += coeff*zbx
        dy += coeff*zby

    return dx, dy

def getZhaoBurgeTerm(polid, x, y):
    
    # Cartesian input x,y; cartesian output x, y, and label.
    # Calls getZ() and delivers area normalized ZB terms
    if polid==0:   # case "S2"  r^0, keep; X translate; 
        result =  getZ(1,x,y), 0.0,            'S2'   # parm=0/13
    elif polid==1:   # case "S3"  r^0, keep; Y translate; 
        result =  0.0, getZ(1,x,y),             'S3'   # parm=1/13
    elif polid==2:   # case "S4"  magnify
        result =  rh*getZ(2,x,y), rh*getZ(3,x,y), 'S4'  # parm=2/13
    elif polid==3:   # case "S5"
        result =  rh*getZ(3,x,y), rh*getZ(2,x,y), 'S5'  # n/a
    elif polid==4:   # case "S6"
        result =  rh*getZ(2,x,y), -rh*getZ(3,x,y), 'S6'  # n/a
    elif polid==5:   # case "S7" 
        result =  0.5*getZ(5,x,y), rh*getZ(4,x,y)-0.5*getZ(6,x,y), 'S7' # parm=3/13
    elif polid==6:   # case "S8"
        result =  rh*getZ(4,x,y)+0.5*getZ(6,x,y), 0.5*getZ(5,x,y), 'S8'  # parm=4/13
    elif polid==7:   # case "S9"
        result =  rh*getZ(5,x,y), rh*getZ(6,x,y), 'S9'      # parm=5/13 but TREFOIL: DROPPING THIS
    elif polid==8:   # case "S10"
        result =  rh*getZ(6,x,y), -rh*getZ(5,x,y), 'S10'    # n/a but is SISTER TO S9 !!!
    elif polid==9:   # case "S11" 
        result =  rh*getZ(8,x,y), rh*getZ(7,x,y), 'S11'   # parm=6/13
    elif polid==10:  # case "S12"
        result =  0.5*getZ(8,x,y)+0.5*getZ(10,x,y), -0.5*getZ(7,x,y)+0.5*getZ(9,x,y), 'S12' # n/a
    elif polid==11:  # case "S13" 
        result =  0.5*getZ(7,x,y)+0.5*getZ(9,x,y), 0.5*getZ(8,x,y)-0.5*getZ(10,x,y), 'S13'  # n/a
    elif polid==12:  # case "S14"
        result =  rh*getZ(10,x,y), -rh*getZ(9,x,y),  'S14'   # n/a
    elif polid==13:  # case "S15"
        result =  rh*getZ(9,x,y), rh*getZ(10,x,y),  'S15'   # n/a
    elif polid==14:  # Case "S16"
        result =  rh*getZ(11,x,y)+0.5*getZ(12,x,y), 0.5*getZ(13,x,y), 'S16'  # n/a; TYPO! 13 not 3
    elif polid==15:  # Case "S17"
        result =  0.5*getZ(13,x,y), rh*getZ(11,x,y)-0.5*getZ(12,x,y), 'S17'  # parm=7/13  TYPO! 13 not 3
    elif polid==16:  # Case "S18"
        result =  0.5*getZ(12,x,y)+0.5*getZ(14,x,y), 0.5*getZ(15,x,y)-0.5*getZ(13,x,y), 'S18' # n/a
    elif polid==17:  # Case "S19"
        result =  0.5*getZ(13,x,y)+0.5*getZ(15,x,y), 0.5*getZ(12,x,y)-0.5*getZ(14,x,y), 'S19'  # n/a
    elif polid==18:  # Case "S20"
        result =  rh*getZ(14,x,y), -rh*getZ(15,x,y), 'S20'  # n/a
    elif polid==19:  # Case "S21"
        result =  rh*getZ(15,x,y), rh*getZ(14,x,y), 'S21', # n/a
    elif polid==20:  # Case "S22"
        result =  rh*getZ(16,x,y), rh*getZ(17,x,y), 'S22'  # parm=8/13
    elif polid==21:  # Case "S23"
        result =  0.5*getZ(17,x,y)+0.5*getZ(19,x,y), 0.5*getZ(16,x,y)-0.5*getZ(18,x,y), 'S23'  # n/a
    elif polid==22:  # Case "S24"
        result =  0.5*getZ(16,x,y)+0.5*getZ(18,x,y), -0.5*getZ(17,x,y)+0.5*getZ(19,x,y), 'S24' # n/a
    elif polid==23:  # Case "S25" 
        result =  0.5*getZ(19,x,y)+0.5*getZ(21,x,y), 0.5*getZ(18,x,y)-0.5*getZ(20,x,y), 'S25'  # n/a
    elif polid==24:  # Case "S26"
        result =  0.5*getZ(18,x,y)+0.5*getZ(20,x,y), -0.5*getZ(19,x,y)+0.5*getZ(21,x,y), 'S26'  # n/a
    elif polid==25:  # Case "S27"
        result =  rh*getZ(21,x,y), rh*getZ(20,x,y), 'S27'                 # n/a
    elif polid==26:  # Case "S28"
        result =  rh*getZ(20,x,y), -rh*getZ(21,x,y), 'S28'                # n/a
    elif polid==27:  #  case "T4"
        result =  rh*getZ(3,x,y), -rh*getZ(2,x,y), 'T4'                   # parm=9/13
    elif polid==28:  # case "T7" 
        result =  rh*getZ(4,x,y)-0.5*getZ(6,x,y), -0.5*getZ(5,x,y), 'T7'  # parm=10/13
    elif polid==29:  # case "T8"
        result =  0.5*getZ(5,x,y), -rh*getZ(4,x,y)-0.5*getZ(6,x,y), 'T8'  # parm=11/13; TWO MINUS SIGNS
    elif polid==30:  # case "T11" 
        result =  rh*getZ(7,x,y), -rh*getZ(8,x,y), 'T11'                  # parm=12/13
    elif polid==31:  # case "T12"
        result =  -0.5*getZ(7,x,y)+0.5*getZ(9,x,y), -0.5*getZ(8,x,y)-0.5*getZ(10,x,y), 'T12' # n/a
    elif polid==32:  # case "T13"
        result =  0.5*getZ(8,x,y)-0.5*getZ(10,x,y), -0.5*getZ(7,x,y)-0.5*getZ(9,x,y), 'T13' # n/a
    else: 
        print("ZhaoBurgeTerm() is exitting because unsatisfied polid = ", polid) 
        quit()
    return result   

#---------end of Zhao-Burge evaluator------------------------


#- Transformation and fit routines

from scipy.optimize import minimize

#- Tranformation in reduced coordinates space, to be used by minimizer
def transform(x, y, scale, rotation, offset_x, offset_y, zbpolids=None, zbcoeffs=None):
    """
    TODO: document
    """
    x = x + offset_x
    y = y + offset_y
    xx = (x*np.cos(rotation) - y*np.sin(rotation))*scale # + offset_x
    yy = (x*np.sin(rotation) + y*np.cos(rotation))*scale # + offset_y

    if zbpolids is not None:
        if zbcoeffs is None:
            raise RuntimeError("need none or both zppolids and zbcoeffs")
        dx, dy = getZhaoBurgeXY(zbpolids,zbcoeffs, xx, yy)
        xx += dx
        yy += dy

    return xx, yy

def fit_scale_rotation_offset(x, y, xp, yp, fitzb=False):
    """
    Fit scale, rotation, offset plus optional Zhao-Burge corrections
    for x,y -> xp,yp.

    TODO: document details
    """

    def func(params, x, y, xp, yp, fitzb):
        scale, rotation, offset_x, offset_y = params[0:4]
        xx, yy = transform(x, y, scale, rotation, offset_x, offset_y)

        if fitzb:
            polids, coeffs, zbx, zby = fitZhaoBurge(xx, yy, xp, yp)
            xx += zbx
            yy += zby

        dr2 = np.sum((xx-xp)**2 + (yy-yp)**2)
        return dr2

    p0 = np.array([1.0, 0.0, 0.0, 0.0])
    p = minimize(func, p0, args=(x, y, xp, yp, fitzb)) #, method='Nelder-Mead')

    scale, rotation, offset_x, offset_y = p.x
   
    if fitzb:
        #- including ZB in every iteration is ~10x better than fitting
        #- scale,rotation,offset, then separately fitting ZB
        xx, yy = transform(x, y, scale, rotation, offset_x, offset_y)
        zbpolids, zbcoeffs, zbx, zby = fitZhaoBurge(xx, yy, xp, yp)
        return scale, rotation, offset_x, offset_y, zbpolids, zbcoeffs
    else:
        return scale, rotation, offset_x, offset_y

def fitZhaoBurge(x, y, xp, yp):
    dx = xp-x
    dy = yp-y

    nx = len(x)

    # here we choose the polynomials
    # 0 = translation along X
    # 1 = translation along Y
    # 2 = magnification
    
    polids = np.array([2, 5,  6,   9,  20,  28, 29,  30],dtype=int)
    
    nzb = polids.size
    H = np.zeros((2*nx, nzb))
    for i,polid in enumerate(polids) :
        zbx, zby, name = getZhaoBurgeTerm(polid, x, y)
        H[0:nx, i] = zbx
        H[nx:, i] = zby

    A = H.T.dot(H)
    b = H.T.dot(np.concatenate([dx, dy]))
    coeffs = np.linalg.solve(A, b)

    zbx, zby = getZhaoBurgeXY(polids, coeffs, x, y)

    return polids, coeffs, zbx, zby
