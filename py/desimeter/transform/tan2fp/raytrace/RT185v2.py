#  RT185v2.py  Updated 22 July 2020 MLL to use DESI-ADC-2 files. 
#  RT185.py  Implementing the precorrection angles for ADCs.
#  RT184.py  Going base zero and preparing for single star single wave probe trace
#  RT183.py  Fraunhofer lines (line 1030) are now base=0 but so far unused. 
#     Keeping base=1 Nsurfs, Nrays, Nglasses
#     Repaired plot color error by calling PLOTCOLORS[], line 1519.
#  RT182.py  Examining the max ADC angular dispersion over the field. Externally called. 
#    set up to use Fraunhofer "i" line 365nm and "t" line 1013nm.    
#  RT181.py  Looking at through-focus for KAP-50100 36x49mm sensor, FVC
#  RT180.py  Studying field distortion of a simple doublet FVC lens
#  RT179-square.py  Like RT179-both.py but calls DESI-SQUARE.RAY.CSV square pupil grid
#  RT179-both.py: outputs both text files of individual rays, and spot diagrams.
#  RT179-text.py: Added text file outputs for each ray, fifty files. 
#  RT179-thru.py: Preparing monochromatic spot diagrams with full obstruction and well sampled pupils.
#     goal is to predict precise dither coordinates at 638nm (DECam Rband LambdaEff).
#     Linear-interpolated DESI-SPOT.MED has this new wavelength. 
#  RT178.py: Implementing Steve kent's ADC roll directions: + is RH about +Z5.
#  RT177.py: Adapting Fraunhofer letters to base=1: 7 wavels but 8 array elements; zero is dummy.
#  RT177.py: commenting out the many print(xxx) statements to run silently, except for trouble
#  RT176.py: rewriting main() to have doTrace() callable as module "import RT176.py"
#  RT175.py: Doing a chroma demo rotating ADCs and making spot diagrams.  
#  RT174.py: Simplified edition, no zernikes, no sagfiles; validated Chroma. 
#  Second half of the EndToEnd study:  Receives a task list, traces rays.
#  Set up to use Tasks.txt, DESI-E2E.OPT, DESI-E2E-RAY, DESI-E2E.MED
#
#  MLL UCB SSL October November 2019
#
#
# OUTLINE:
#   Library imports:    line 50 onward
#   Macro definitions:   line 70...
#   List definitions:    lines 180...
#   Table parsing tools:  lines 220..
#   Math helpers:          lines 460...  
#   Surface generators:     lines 540...
#   Slope generators:       lines 580..
#   Interceptors:           lines 640...
#   Validator:               lines 750...
#   Redirectors:              lines 790...
#   Coordinate changers:       lines 860...
#   Text output formatters:     lines 925...
#   Input CSV file unpackers:    lines 1050...
#   Ray tracing:                lines 1350...
#   Output helpers:              lines 1450...
#   Externally called methods:    lines 1620...
#  
# 
#
#




#-------------------RAY TRACING IMPORT ZONE-----------------------

import math
import numpy as np
import matplotlib.pyplot as plt       # for showing PSF output
import csv                            # for interpreting CSV files; 'rU' = universal EOL recognizer
import matplotlib.colors as mcolors   # for colorizing pixels with specified gamma


    
def printArray(name, array):
    print(name, end='')
    # if len(array) == 1:
    for item in array:
        print('{:14.6f}'.format(item), end='')
    print()   
    
PROGNAME   = 'RT185'



#--------RAY TRACING SET UP MACROS FOR INDEXING THE ARRAYS-------------
#----------------group 0 "OBASIC"-----------
OABSENT     = -1
OINDEX      = 0   # user-supplied refraction index
OACTIONTYPE = 1   # user-supplied surface  action type: lens, mirror, etc
OODIAM      = 2   # user-supplied limiting Outer Diameter
OIDIAM      = 3   # user-supplied inner diameter
OCURVE      = 4   # user-supplied curvature
OASPH       = 5   # user-supplied asphericity

#-------------------group 1 "OXYZ"---------
OX          = 6   # user-supplied X coordinate of surface vertex, lab frame
OY          = 7   # user-supplied Y coordinate of surface vertex, lab frame
OZ          = 8   # user-supplied Z coordinate of surface vertex, lab frame

#-------------------group 2 "OTPR" be sure to call setEulers()---------------
OTILT     = 9   # user-supplied tilt angle about lab X axis, degrees
OPITCH    = 10  # user-supplied pitch angle about tilted Y axis, degrees
OROLL     = 11  # user-supplied roll angle about pitched and tilted Z axis, degrees

#------------------group 3 "OPOLY"------------
OA1       = 12  # user-supplied polynomial coefficients
OA2       = 13 
OA3       = 14
OA4       = 15 
OA5       = 16
OA6       = 17
OA7       = 18
OA8       = 19
OA9       = 20
OA10      = 21 
OA11      = 22
OA12      = 23
OA13      = 24
OA14      = 25

#------------------spider macro definitions------------
ONSPIDER   = 30  # number of equally spaced spider legs
OWSPIDER   = 31  # width of each spider leg
OMAXINPUT  = 32

#-----------------group 4 "EULER" internal access only---------
OE11      = 41  # internal Euler matrix element
OE12      = 42  # internal Euler matrix element
OE13      = 43  # internal Euler matrix element
OE21      = 44  # internal Euler matrix element
OE22      = 45  # internal Euler matrix element
OE23      = 46  # internal Euler matrix element
OE31      = 47  # internal Euler matrix element
OE32      = 48  # internal Euler matrix element
OE33      = 49  # internal Euler matrix element
OFINAL    = 50

#------------RAY TRACING optics action code values---------

OLENSACTION    = 0
OMIRRORACTION  = 1
ORETROACTION   = 2
OIRISACTION    = 3
OCBINACTION    = 4
OCBOUTACTION   = 5
actionLookup   = ['OLENSACTION', 'OMIRRORACTION', 'ORETROACTION', 'OIRISACTION', 'OCBINACTION', 'OCBOUTACTION']



#-----------RAY TRACING FIELD HEADER IDENTIFIERS----------------
#-----these apply at every surface: 0=Raystarts; F=final surface.
#-------------parser must supply the surface number.-------------

RX      = 0  #uppercase: lab frame
RY      = 1
RZ      = 2
RU      = 3
RV      = 4
RW      = 5
RWAVE   = 6
RXG     = 7   # goal field
RYG     = 8   # goal field
RFINALINPUT = 8

Rx      = 7  # lowercase = vertex frame
Ry      = 8
Rz      = 9
Ru      = 10
Rv      = 11
Rw      = 12
RFINAL  = 13

#--------RAY TRACING ray failure classifiers--------------------------    

OK       = 0   # no failure
BACK     = 1   # intercept failure; only roots are d<=0, i.e. backwards
MISS     = 2   # intercept failure
ODIAM    = 3   # outer diameter validate failure
IDIAM    = 4   # inner diameter validate failure
SPI      = 5   # spider leg hit killed ray
TIR      = 6   # refraction failure
EDGE     = 7   # interpolation beyond safe boundary
UNKNOWN  = 8   # no valid action 

failures = ['OK', 'BACK', 'MISS', 'ODIAM', 'IDIAM', 'SPI', 'TIR', 'EDGE', 'UNKNOWN']







#------------------RAY TRACING LIST AND INDEX DEFINITIONS-----------------
#-------------Helpers that just read these need not declare them global---
#-----Helpers that create or modify these MUST declare them global--------
#--------Here, the .CSV files have no ruler line; parsing is by field-----

#----RAY TRACING Optics params and lists------------------------
Nsurfs          = 0                 # how many optical surfaces, modified by guide number.
Onfields        = 0                 # how many data fields in the table

Ostrings        = []                # list of strings starting with title string
Otable          = []                # list of lists of fields: 2D, starting with row=2
Oheaders        = []                # list of strings, one per field
OhasZern        = []                # one-based list of booleans
OhasPoly        = []                # one-based list of booleans
OglassNames     = []                # one-based list of refractive indices and glass Names. 
OneedsMedia     = False             # if True, nonnumerical glass names will require a validated .MED table
Oarray          = np.empty([0,0])   # Working, possibly modified, optical specification [allSurfaces, allParms]

#----RAY TRACING Ray arrays--------------------------
Nrays         = 0                    # how many ray starts, modified by guide number
Rnfields      = 0                    # how many data fields in the table
Rstrings      = []                   # list of strings, 1D, starting with title string
Rtable        = []                   # list of lists of fields: 2D, starting with row=2
Rheaders      = []                   # list of strings
RItoF         = []                   # lookup index-to-field list ONLY FOR RAY INPUT FIELDS
RFtoI         = []                   # lookup field-to-index list.
Raystarts     = np.zeros([0,0])      # [nrays,nparms]; updated in getRaysCSV()
Rarray        = np.zeros([0,0,0])    # [nlines,nsurfs,nparms]; updated in getRaysCSV()
RwaveNames    = []                   # list of wavelengths from "@wave" field

#----RAY TRACING Media lists----------------------------
Nglasses    = 0                      # how many glasses, modified by guide number
Mnfields    = 0                      # how many data fields in the table
Mstrings    = []                     # list of strings starting with title string
Mtable      = []                     # list of lists of fields: 2D, starting with row=2
Mheaders    = []                     # list of strings = wavelength Names
MwaveNames  = []                     # list of wavelength Names from headers > 0
MglassNames = []                     # list of glass Names from column zero
Marray      = np.empty([0,0])        # 2D numpy array




#------RAY TRACING TABLE PARSING TOOLS------------------

def cuberoot(x):
    # allows negative arguments
    return math.copysign(math.pow(abs(x), 1./3.), x)

def makeTen(s):
    while len(s)<10:
        s += ' '
    return s

def getFloatValue(s):
    try:
        x = float(s)
    except ValueError:
        x = -0.0
    return x

def fixup(u,v,w):
    arg = 1 - u*u - v*v
    if arg < 0.:
        print('fixup() has u,v too large. Quitting.', u, v)
        quit()
    if w>= 0:  # includes -0.0
        return np.sqrt(arg)
    else:
        return -np.sqrt(arg)

def suckInt(s):
    # extracts a positive integer from a messy string
    slen = len(s)
    if slen<1:
        return 0
    numstring = ""
    i = 0
    while i<slen and not s[i].isdigit():
        i += 1
    while i<slen and s[i].isdigit():
        numstring += s[i]
        i += 1
    if len(numstring) ==0:
        return 0
    try:
        result = int(numstring)
        return result
    except ValueError:
        return 0
    
def getActionType(snippet):
    # this number can be put right into Oarray[OACTIONTYPE field]
    length = len(snippet)
    if length < 1:
        return OLENSACTION
    c0 = snippet[0].upper()
    c2 = ' '
    if length > 2:
        c2 = snippet[2].upper()
    if c0 == 'M':
        return OMIRRORACTION
    if c0 == 'R':
        return ORETROACTION
    if c0 == 'I':     # includes spider
        return OIRISACTION
    if c0 == 'C':    #coordinate break
        if c2 == 'I':
            return OCBINACTION
        if c2 == 'O':
            return OCBOUTACTION
    return OLENSACTION

def getOpticsAttribute(header):
    # From a column header, return its math array column index.
    # Surface numbers are row numbers, they do not come from headers.
    # Includes Zernike coefficients, SagMaps.
    guidenum = suckInt(header)
    header += "      "
    c0 = header[0]
    c1 = header[1]
    c2 = header[2]
    c3 = header[3]
    c0up = c0.upper()
    c1up = c1.upper()
    c2up = c2.upper()
    c3up = c3.upper()
    
    if c0=='D':
        return OODIAM
    elif c0=='d':
        return OIDIAM
    elif (c0up=='T' and c1up=='Y') or (c0up=='L' and c1up=='E') or (c0up=='M' and c1up=='I'):
        # rint 'getOpticsIndex() header and OACTIONTYPE:', header, OACTIONTYPE
        return OACTIONTYPE
    elif c0up=='I' or c0up=='G':
        return OINDEX
    elif c0up=='X':
        return OX
    elif c0up=='Y':
        return OY
    elif c0up=='Z' and not c1up=='E':
        return OZ
    elif c0up=='P':
        return OPITCH
    elif c0up=='T' and c1up=='I':
        return OTILT
    elif c0up=='R':
        return OROLL
    elif c0up=='C':
        return OCURVE
    elif c0up=='A' and c1up=='S':
        return OASPH
    elif c0up=='A' and guidenum>=0 and guidenum<15:
        sum = OA1 - 1 + guidenum
        # print 'getOpticsIndex() finds polynomial field whose index = ', sum
        return sum
    elif c0up=='Z' and c1up=='E' and guidenum<36:
        sum = OZ0 + guidenum
        # print 'getOpticsIndex() finds Zernike field whose index = ', sum
        return sum
    elif c0up=='S' and c3up=='F':   # sag header
        return OSAGFILE
    elif c0up=='S' and c3up=='M':
        return OSAGMULT
    elif c0up=='S' and c3up=='S':
        return OSAGSTEP
    elif c0up=='N' and c1up=='S':
        return ONSPIDER
    elif c0up=='W' and c1up=='S':
        return OWSPIDER
    else:
        return -1
    
def getRayStartAttribute(header):
    # Returns the input field parameter, or -1 if not an input field.
    # print 'getRayStart(header) has header = ', header
    if len(header) < 1:
        return -1
    header = header.upper()
    c0 = header[0]
    if c0=='@':
        return RWAVE
    c1 = ' '
    if len(header) > 1:
        c1 = header[1]
    # For AutoRay, X0... are ray inputs, XG... are ray goals   
    # but sometimes I want to accept any Xxxxx or Yxxxx as a goal. 
    # if len(header) < 2:
    #    return -1
    # if header[1] != '0':  
    #     return -1
    
    if c0=='X':
        if c1 == 'G':
            # print '......XG detected; returning RXG = ', RXG
            return RXG
        if c1 == '0':
            return RX
        return -1
    if c0=='Y':
        if c1 == 'G':
            return RYG
        if c1 == '0':
            return RY
        return -1
    if c0=='Z':
        if c1 == '0':
            return RZ
        return -1
    if c0=='U':
        if c1 == '0':
            return RU
        return -1
    if c0=='V':
        if c1 == '0':
            return RV
        return -1
    if c0=='W':
        if c1 == '0':
            return RW
        return -1
    return -1
    

def findGlassRow(glassname):
    # Search MglassNames to find a given glassname.
    # Return -1 if not found. 
    for kglass in range(1, len(MglassNames)):
        if glassname == MglassNames[kglass]:
            # print 'findGlassRow has glassname, kglass = ', glassname, kglass
            return kglass
    print('RT: findGlassRow() not found. Quitting')
    quit()
    
    
def findWaveColumn(wavename):
    # Search MwaveNames trying to locate a given wavename.
    # Skip column zero: it is the glass name header.
    # quit() if not found.
    # print('findWaveColumn() is searching for wavename = ', wavename)
    for col in range(1, len(MwaveNames)):  # ignore column zero.
        if wavename == MwaveNames[col]:  
            # print 'findWaveColumn has wavename, column = ',wavename, col
            return col
    print('RT: findWaveColumn() for wavename = ' + wavename +' not found. Quitting.')
    quit()
    

def findRefraction(iray, jsurf):
    glassname = OglassNames[jsurf]      # numbering 1...Nsurfs
    if len(glassname) < 1:
        return 1.0                      # assumes blank = air = 1.0    
    try:
        result = float(glassname)
        return result
    except ValueError:                  # need a Media table     
        if Nglasses < 1:
            print('Media table is needed for glassname = ', glassname)
            quit()
        wavename = RwaveNames[iray]         # numbering 1...Nrays
        # print('findRefraction() is searching for .RAY wavename = ', wavename)
        result = 1.0
        mediarow = findGlassRow(glassname)
        # print('findRefraction() is accessing mediarow = ', mediarow)
        if mediarow < 0:
            print('findRefraction() is quitting since no GlassName = ', glassname)
            quit()
        mediacol = findWaveColumn(wavename)
        if mediacol < 1:
            print('findRefraction() is quitting since no WaveName = ', wavename)
            quit()
        result = float(Marray[mediarow][mediacol])
    return result
    



#--------RAY TRACING MATH HELPERS----------------

def deg(radians):
    return math.degrees(radians)
    
def rad(degrees):
    return math.radians(degrees)

def isMinusZero(x):
    return x==0. and np.signbit(x)==True

def getBothRoots(A, B, C):
    # Solves for the real roots of a quadratic function.
    # Returns 0, 1, or two roots as a tuple: "-0.0" means failed.
    # Method is Press et al 'Numerical Recipes' 2nd edition p.183
    if A==0.0 and B!= 0.0:
        return -C/B, -0.0
    if B==0.0 and A!=0.0:
        if C/A>0:
            return np.sqrt(C/A), -np.sqrt(C/A)
        else:
            return -0.0, -0.0
    if C==0 and A!= 0.0:
        return -B/A, -0.0
    D = B*B-4*A*C
    if D <= 0.0:
        return -0.0, -0.0
    Q = -0.5*(B + np.sign(B)*np.sqrt(D))
    return Q/A, C/Q

def setEulers():  # call this after any change in OTILT, OPITCH, or OROLL
    for j in range(1, Nsurfs+1):      
        ct = np.cos(np.radians(Oarray[j, OTILT]))
        st = np.sin(np.radians(Oarray[j, OTILT]))
        cp = np.cos(np.radians(Oarray[j, OPITCH])) 
        sp = np.sin(np.radians(Oarray[j, OPITCH])) 
        cr = np.cos(np.radians(Oarray[j, OROLL])) 
        sr = np.sin(np.radians(Oarray[j, OROLL]))  
        Oarray[j,OE11] = cr*cp;               # X <- x; M11
        Oarray[j,OE12] = -sr*cp;              # X <- y; M12
        Oarray[j,OE13] = sp;                  # X <- z; M13
        Oarray[j,OE21] = cr*sp*st + sr*ct;    # Y <- x; M21
        Oarray[j,OE22] = cr*ct - sr*sp*st;    # Y <- y; M22
        Oarray[j,OE23] = -cp*st;              # Y <- z; M23
        Oarray[j,OE31] = -cr*sp*ct + sr*st;   # Z <- x; M31
        Oarray[j,OE32] = sr*sp*ct + cr*st;    # Z <- y; M32
        Oarray[j,OE33] = cp*ct;               # Z <- z; M33   
 
def dotproduct(abc, xyz):
    # returns the dot product of two triplets
    return abc[0]*xyz[0] + abc[1]*xyz[1] + abc[2]*xyz[2]
    
    
def crossproduct(abc, xyz):
    # returns the cross product of two triplets
    product = np.zeros(3)
    product[0] = abc[1]*xyz[2] - abc[2]*xyz[1]
    product[1] = abc[2]*xyz[0] - abc[0]*xyz[2]
    product[2] = abc[0]*xyz[1] - abc[1]*xyz[0]
    return product    
    
def normalize(norm):
    # modifies given host array.
    len = np.sqrt(norm[0]**2 + norm[1]**2 + norm[2]**2)
    if len==0:
        print("cannot normalize a zero vector")
        return
    norm[0] /= len
    norm[1] /= len
    norm[2] /= len

def testUVW(iray, jsurf):
    err = Rarray[iray, jsurf, RU]**2  \
        + Rarray[iray, jsurf, RV]**2  \
        + Rarray[iray, jsurf, RW]**2  \
        - 1.0
    if math.fabs(err) > 1E-14:
        print('UVW normalization error at iray, surf = ', iray, jsurf, err)

def isNegZero(x):
    return x==0. and np.signbit(x)==True

    

#-----RAY TRACING SURFACE GENERATORS--------------------

def getZtotal(iray, jsurf, d):
    # "d" is a positive trial distance along current ray
    z = getZconic(iray, jsurf, d)
    if OhasPoly[jsurf]:
        z += getZpoly(iray, jsurf, d)
    if OhasZern[jsurf]:
        z += getZzern(iray, jsurf, d)
    if OhasSag[jsurf]:
        z += getZsag(iray, jsurf, d)
    return z   

    
def getZconic(iray, jsurf, d):
    #  coordinates here are local "vertex frame" values.
    # "d" is a positive trial distance along current ray being tested here.
    x = Rarray[iray, jsurf, Rx] + d * Rarray[iray, jsurf, Ru]
    y = Rarray[iray, jsurf, Ry] + d * Rarray[iray, jsurf, Rv]
    r2 = x*x + y*y
    s = Oarray[jsurf, OASPH] + 1.0
    c = Oarray[jsurf, OCURVE]
    numer = c*r2
    arg = 1 - s*c*c*r2
    if arg < 0:
        print('negative argument found by getZconic(); flange case failure code -0.0')
        return -0.0, MISS  # failure code
    denom = 1 + np.sqrt(arg)
    zconic = numer/denom
    return zconic, OK

def getZpoly(iray, jsurf, d):
    #  coordinates here are local vertex frame values.
    x = Rarray[iray, jsurf, Rx] + d * Rarray[iray, jsurf, Ru]
    y = Rarray[iray, jsurf, Ry] + d * Rarray[iray, jsurf, Rv]
    r = np.sqrt(x*x + y*y)
    product = 1.0
    sum = 0.0
    for attrib in range(OA1, OA14+1):   # OA1=11 ... OA14=24
        product *= r
        sum += Oarray[jsurf, attrib] * product
    return sum, OK
    

#----RAY TRACING SURFACE SLOPES AND NORMALS---------    
    
def getNormal(iray, jsurf):
    # There are two surface normals. Should not matter. I always use the one with Nz>0.
    gx, gy = gradTotal(iray, jsurf)
    normal = np.array([-gx, -gy, 1.0])
    normalize(normal)
    return normal

def gradTotal(iray, jsurf):
    gx, gy = gradConic(iray, jsurf)
    if OhasPoly[jsurf]:
        px, py = gradPoly(iray, jsurf)
        gx += px
        gy += py
    # print 'gradTotal() is returning gx, gy = ', gx, gy
    return gx, gy
    
def gradConic(iray, jsurf):
    s = Oarray[jsurf, OASPH] + 1
    c = Oarray[jsurf, OCURVE]
    if c==0:               # plano case
        return 0.0, 0.0
    x = Rarray[iray, jsurf, Rx]
    y = Rarray[iray, jsurf, Ry]
    r2 = x*x + y*y
    arg = 1.0 - s*c**2*r2
    if arg <= 0:           # flange case
        # print '   gradConic() is returning the flange case.'
        return -0.0, -0.0
    coef = c/np.sqrt(arg)  # conic case
    gx = x*coef
    return x*coef, y*coef
    
def gradPoly(iray, jsurf):
    x = Rarray[iray, jsurf, Rx]
    y = Rarray[iray, jsurf, Ry]
    r = np.sqrt(x*x + y*y)
    if r==0:
        return 0.0, 0.0
    product = 1.0
    dzdr = 0.0
    for index in range(OA1, OA14+1):   # OA1=11 ... OA14=24
        coef = 1 + index - OA1
        dzdr += coef * Oarray[jsurf, index] * product
        product *= r
    return (x/r)*dzdr, (y/r)*dzdr 








#----------RAY TRACING INTERCEPTOR TOOLS------------------------

def intercept(iray, jsurf):
    # Works entirely in local coordinate frame Rx,Ry,Rz,Ru,Rv,Rw.
    # These numbers come entirely out of labToVx() routine.
    # zlocal = Rarray[iray, jsurf, Rz]   #
    # wlocal = Rarray[iray, jsurf, Rw]   # always negative, as planned.
    if OhasPoly[jsurf]:
        return higherIntercept(iray, jsurf)
    return conicIntercept(iray, jsurf)

def func(iray, jsurf, d):  
    # Function whose root is to be found using Newton's method.
    zc, code = getZconic(iray, jsurf, d)
    if code!=OK:
        return -0.0, code
    # print 'zc = ', zc 
    zp, code = getZpoly(iray, jsurf, d)
    if code!=OK:
        return -0.0, code
    # print 'zp = ', zp
    z0 = Rarray[iray, jsurf, Rz]
    w = Rarray[iray, jsurf, Rw]
    sum = zc + zp - (z0 + w*d)
    # print 'func() gets total sum = ', sum
    return sum, OK
    
def deriv(iray, jsurf, d):
    # Estimator of the derivative of func() for Newton's method.
    DELTA = 0.00001   # should be made adaptive
    fplus, status = func(iray, jsurf, d+DELTA)
    fminus, status = func(iray, jsurf, d-DELTA)
    return (fplus - fminus)/(2.*DELTA)

def higherIntercept(iray, jsurf):
    # First do a conic intercept to get close to the correct root. 
    # print '\nHigherIntercept starting its conic intercept...'
    status  = conicIntercept(iray, jsurf)
    if status !=OK:
       return status
    # Set up for using the Newton method.
    # print '\nHigherIntercept() is setting up for Newton rootfinder'
    d = 0
    niters = 0
    while True:   # Newton rootfinder
        niters += 1
        f, status = func(iray, jsurf, d)
        if status!=OK:
            return status
        slope = deriv(iray, jsurf, d)
        d -= f/slope
        # print '    Newton rootfinder: niters, d, f: ', niters, d, f
        if abs(f) < 1E-12 or niters > 8:
            break;
    Rarray[iray, jsurf, Rx] = Rarray[iray, jsurf, Rx] + d*Rarray[iray, jsurf, Ru]
    Rarray[iray, jsurf, Ry] = Rarray[iray, jsurf, Ry] + d*Rarray[iray, jsurf, Rv]
    Rarray[iray, jsurf, Rz] = Rarray[iray, jsurf, Rz] + d*Rarray[iray, jsurf, Rw]
    return status
    
def conicIntercept(iray, jsurf):   
    s = Oarray[jsurf, OASPH] + 1.0
    c = Oarray[jsurf, OCURVE]
    
    # Note: labtovx() will have already set the vertex-frame ray starts.
    
    x = Rarray[iray, jsurf, Rx]
    y = Rarray[iray, jsurf, Ry]
    z = Rarray[iray, jsurf, Rz]
    u = Rarray[iray, jsurf, Ru]
    v = Rarray[iray, jsurf, Rv]
    w = Rarray[iray, jsurf, Rw]
    sos = u*u + v*v + w*w
    err = sos - 1.0
    if abs(err) > 1E-12:
        print('Yikes, faulty direction cosines at jsurf =', jsurf)
    
    A = c*(u*u + v*v + s*w*w)
    B = 2*c*(x*u + y*v + s*z*w) - 2*w
    C = c*(x*x + y*y + s*z*z) - 2*z
    r1, r2 = getBothRoots(A, B, C)
    rLessPositive = min(r1, r2)
    rMorePositive  = max(r1, r2)
    if rMorePositive<=0. and rLessPositive<0.:
        return BACK
        
    d=-0.0
    sheetLessPositive = s*c*(z + w * rLessPositive)
    sheetMorePositive  = s*c*(z + w * rMorePositive)
    
    # Always try the shorter path first...
    if sheetLessPositive<1.0 and rLessPositive>0.0:  # if OK, this is our winner.
        d=rLessPositive
    elif sheetMorePositive<1.0 and rMorePositive>0.0:  # else maybe this is our winner. 
        d=rMorePositive
    if d<=0.0:                             # Neither? well then... failed.
        # print 'intercept failed: d1, d2 = ', rLessPositive, rMorePositive
        return MISS
        
    # Now we have a good ray segment length "d"  -- SO PROPAGATE.
    # print 'interceptor has d = ', d
    Rarray[iray, jsurf, Rx] = Rarray[iray, jsurf, Rx] + d*u
    Rarray[iray, jsurf, Ry] = Rarray[iray, jsurf, Ry] + d*v
    Rarray[iray, jsurf, Rz] = Rarray[iray, jsurf, Rz] + d*w  
    return OK
    
    
    
    
    
    
#---------RAY TRACE VALIDATOR with spider----------------------

def validate(iray, jsurf):
    # Is intercept within {Diameter, diameter} pupil annulus?
    # print 'Starting validation for jsurf = ', jsurf
    x = Rarray[iray, jsurf, Rx]       # Rx is index for local frame ray x
    y = Rarray[iray, jsurf, Ry]       # Ry is index for local frame ray y
    r = np.sqrt(x*x + y*y)            # ray's local radius off axis
    Diam = Oarray[jsurf, OODIAM]      # outer Diameter; -0.0 if absent
    if Diam > 0.0:
        if r > 0.5*Diam:
            # print("clobbered by OODIAM.")
            return ODIAM
    if r < 0.5*Oarray[jsurf, OIDIAM]:  # inner diameter; -0.0 if absent
        # print("clobbered by IDIAM")
        return IDIAM

    nlegs = int(Oarray[jsurf, ONSPIDER])
    halfwidth = 0.5*Oarray[jsurf, OWSPIDER]
    if nlegs < 1 or halfwidth <= 0.0:   # no spider declared
        return OK
        
    # Test for rectangular spider leg blockage. Don't rotate the spider: rotate the ray. 
    rolldeg = Oarray[jsurf, OROLL]                      # spider roll pattern, degrees CCW from +X
    raydeg = deg(np.arctan2(y,x))                       # degrees CCW from +X;  0/0 is OK here
    # print('Spider at iray, jsurf: nlegs, halfwidth, rolldeg, raydeg = {:6d}{:3d}{:3d}{:8.2f}{:8.2f}{:8.2f}'.format(iray,jsurf, nlegs, halfwidth, rolldeg, raydeg))
    for ileg in range(0, nlegs):
        diff = ileg*360.0/nlegs + rolldeg - raydeg      # degrees between ray and leg
        xtemp = math.fabs(r * math.sin(rad(diff)))
        ytemp = r * math.cos(rad(diff))
        # print('...spider ileg, rolldeg, raydeg, diffdeg, xtemp, ytemp = {:3d}{:8.3f}{:8.3f}{:8.3f}{:8.1f}{:8.1f}'.format(ileg, rolldeg, raydeg, diff, xtemp, ytemp))
        if xtemp < halfwidth and ytemp > 0.:
            # print('....Finding that xtemp<halfwidth and ytemp>0, so declaring ray failure of type SPI')
            return SPI
    return OK





#-------RAY TRACING REDIRECTORS-------------------

def redirect(iray, jsurf):
    # This is the switchyard to the detail redirectors..
    if jsurf == Nsurfs:
        return OK                 # no redirection at final surface.
        
    u = Rarray[iray, jsurf, Ru]   # input directions to be modified
    v = Rarray[iray, jsurf, Rv]
    w = Rarray[iray, jsurf, Rw]
    action = int(Oarray[jsurf, OACTIONTYPE])
    # print('====edirect starting with u,v,w,action = {:12.6f}{:12.6f}{:12.6f}{:4d}'.format(u,v,w,action)+"  "+actionLookup[action])
    if action == OIRISACTION:     # cannot fail.
        return OK

    if action == OMIRRORACTION:   # cannot fail.
        normal = getNormal(iray, jsurf)
        dotp = u*normal[0] + v*normal[1] + w*normal[2]
        u = Rarray[iray, jsurf, Ru] = Rarray[iray, jsurf, Ru] - 2*dotp*normal[0]
        v = Rarray[iray, jsurf, Rv] = Rarray[iray, jsurf, Rv] - 2*dotp*normal[1]
        w = Rarray[iray, jsurf, Rw] = Rarray[iray, jsurf, Rw] - 2*dotp*normal[2]
        return OK

    if action == OLENSACTION:   # can fail via TIR.    
        numer = findRefraction(iray, jsurf)
        if numer==0.0:
            numer = 1.0
        denom = findRefraction(iray, jsurf+1)
        if denom==0.0:
            denom = 1.0
        # print 'numer, denom = ', numer, denom
        mu = numer/denom
        normal = getNormal(iray, jsurf)
        kinput = np.array([u, v, w])
        kparallel = crossproduct(normal, crossproduct(kinput, normal))
        
        kparallel = mu * kparallel                         # vector equation 
        kparallelSQ = dotproduct(kparallel, kparallel)     # scalar equation 
        kperpSQ = 1.- kparallelSQ                          # Pythagoras
        if kperpSQ <= 0.0:
            return TIR
        kperpmagnitude = np.sqrt(kperpSQ)                  # scalar equation 
        perpsign = np.sign(dotproduct(normal, kinput))
        kperp = perpsign*kperpmagnitude * normal           
        kout = kparallel + kperp                           # vector equation 
        Rarray[iray, jsurf, Ru] = kout[0]
        Rarray[iray, jsurf, Rv] = kout[1]
        Rarray[iray, jsurf, Rw] = kout[2]
        return OK
        
    if action == ORETROACTION:    # cannot fail
        Rarray[iray, jsurf, Ru] *= -1.0
        Rarray[iray, jsurf, Rv] *= -1.0
        Rarray[iray, jsurf, Rw] *= -1.0
        return OK
        
    if action == OCBINACTION or action == OCBOUTACTION:
         return OK
         
    return UNKNOWN
    
    
    
    
    
    
    

#---------RAY TRACING COORDINATE CHANGERS-------------------

def labtovx(iray, jsurf):
    global Rarray
    # Assumes that raystarts have been loaded into Rarray[iray, 0].
    # Coordinate frame changer, moves from jsurf-1 LAB to current jsurf LOCAL.
    # Matrix OE converts local to lab coordinates; must transpose it here.
    # M.Lampton STELLAR SOFTWARE (C) 1989, 2003, 2017

    Xprev = Rarray[iray, jsurf-1, RX]
    Yprev = Rarray[iray, jsurf-1, RY]
    Zprev = Rarray[iray, jsurf-1, RZ]
    Uprev = Rarray[iray, jsurf-1, RU]
    Vprev = Rarray[iray, jsurf-1, RV]
    Wprev = Rarray[iray, jsurf-1, RW] # if forward: Rw is positive.  Reverse: Rw is negative.

    xLocal = Xprev - Oarray[jsurf, OX]
    yLocal = Yprev - Oarray[jsurf, OY]
    zLocal = Zprev - Oarray[jsurf, OZ]  # if forward: Rz is negative.  Reverse: Rz is positive.  
    
    x = Rarray[iray,jsurf,Rx] = xLocal*Oarray[jsurf,OE11] + yLocal*Oarray[jsurf,OE21] + zLocal*Oarray[jsurf,OE31]
    y = Rarray[iray,jsurf,Ry] = xLocal*Oarray[jsurf,OE12] + yLocal*Oarray[jsurf,OE22] + zLocal*Oarray[jsurf,OE32]
    z = Rarray[iray,jsurf,Rz] = xLocal*Oarray[jsurf,OE13] + yLocal*Oarray[jsurf,OE23] + zLocal*Oarray[jsurf,OE33]

    u = Rarray[iray,jsurf,Ru] = Uprev*Oarray[jsurf,OE11] + Vprev*Oarray[jsurf,OE21] + Wprev*Oarray[jsurf,OE31]
    v = Rarray[iray,jsurf,Rv] = Uprev*Oarray[jsurf,OE12] + Vprev*Oarray[jsurf,OE22] + Wprev*Oarray[jsurf,OE32]
    w = Rarray[iray,jsurf,Rw] = Uprev*Oarray[jsurf,OE13] + Vprev*Oarray[jsurf,OE23] + Wprev*Oarray[jsurf,OE33]

    return OK

def vxtovx(iray, jsurf):
    # used only by CBout coordinate break. No math; it just copies locals.
    Rarray[iray, jsurf, Rx] = Rarray[iray, jsurf-1, Rx]
    Rarray[iray, jsurf, Ry] = Rarray[iray, jsurf-1, Ry]
    Rarray[iray, jsurf, Rz] = Rarray[iray, jsurf-1, Rz]    
    Rarray[iray, jsurf, Ru] = Rarray[iray, jsurf-1, Ru]
    Rarray[iray, jsurf, Rv] = Rarray[iray, jsurf-1, Rv]
    Rarray[iray, jsurf, Rw] = Rarray[iray, jsurf-1, Rw]    

def vxtolab(iray, jsurf):
    # Coordinate frame changer at a single surface.
    # Here the Euler matrix is used directly, local to lab conversion. 
    # M.Lampton STELLAR SOFTWARE (C) 1989, 2003, 2017

    x = Rarray[iray, jsurf, Rx]
    y = Rarray[iray, jsurf, Ry]
    z = Rarray[iray, jsurf, Rz]    
    u = Rarray[iray, jsurf, Ru]
    v = Rarray[iray, jsurf, Rv]
    w = Rarray[iray, jsurf, Rw]
    
    Rarray[iray, jsurf, RU] = u*Oarray[jsurf,OE11] + v*Oarray[jsurf,OE12] + w*Oarray[jsurf,OE13]
    Rarray[iray, jsurf, RV] = u*Oarray[jsurf,OE21] + v*Oarray[jsurf,OE22] + w*Oarray[jsurf,OE23]
    Rarray[iray, jsurf, RW] = u*Oarray[jsurf,OE31] + v*Oarray[jsurf,OE32] + w*Oarray[jsurf,OE33]

    Rarray[iray, jsurf, RX] = x*Oarray[jsurf,OE11] + y*Oarray[jsurf,OE12] + z*Oarray[jsurf,OE13]
    Rarray[iray, jsurf, RY] = x*Oarray[jsurf,OE21] + y*Oarray[jsurf,OE22] + z*Oarray[jsurf,OE23]
    Rarray[iray, jsurf, RZ] = x*Oarray[jsurf,OE31] + y*Oarray[jsurf,OE32] + z*Oarray[jsurf,OE33]
    
    Rarray[iray, jsurf, RX] = Rarray[iray, jsurf, RX] + Oarray[jsurf, OX]
    Rarray[iray, jsurf, RY] = Rarray[iray, jsurf, RY] + Oarray[jsurf, OY]
    Rarray[iray, jsurf, RZ] = Rarray[iray, jsurf, RZ] + Oarray[jsurf, OZ]
    return OK




""" #---------RAY TRACING NUMERICAL TEXT DISPLAY TOOLS----------------

def showEuler(jsurf):
    eulerlist = [[Oarray[jsurf, OE11], Oarray[jsurf, OE12], Oarray[jsurf, OE13]],
                 [Oarray[jsurf, OE21], Oarray[jsurf, OE22], Oarray[jsurf, OE23]],
                 [Oarray[jsurf, OE31], Oarray[jsurf, OE32], Oarray[jsurf, OE33]]]
    euler = np.array(eulerlist)
    print(euler)
    
def displayInput(iray):
    X = Raystarts[iray, RX]
    Y = Raystarts[iray, RY]
    Z = Raystarts[iray, RZ]    
    U = Raystarts[iray, RU]
    V = Raystarts[iray, RV]
    W = Raystarts[iray, RW] 
    print('+++ Input:   iray,   XYZUVW: {:3d}{:16.8f}{:16.8f}{:16.8f}{:16.8f}{:16.8f}{:16.8f}'.format(iray,X,Y,Z,U,V,W))

def displayLocal(iray, jsurf):
    x = Rarray[iray, jsurf, Rx]
    y = Rarray[iray, jsurf, Ry]
    z = Rarray[iray, jsurf, Rz]    
    u = Rarray[iray, jsurf, Ru]
    v = Rarray[iray, jsurf, Rv]
    w = Rarray[iray, jsurf, Rw] 
    print('*** Output: howfar,  xyzuvw:{:3d}{:16.8f}{:16.8f}{:16.8f}{:16.8f}{:16.8f}{:16.8f}'.format(jsurf,x,y,z,u,v,w))

def displayLabs(iray, jsurf):
    X = Rarray[iray, jsurf, RX]
    Y = Rarray[iray, jsurf, RY]
    Z = Rarray[iray, jsurf, RZ]    
    U = Rarray[iray, jsurf, RU]
    V = Rarray[iray, jsurf, RV]
    W = Rarray[iray, jsurf, RW] 
    print('*** Output: howfar,  XYZUVW:{:3d}{:16.8f}{:16.8f}{:16.8f}{:16.8f}{:16.8f}{:16.8f}'.format(jsurf,X,Y,Z,U,V,W))

def displayLongOutput(iray, jsurf):
    X = Rarray[iray, jsurf, RX]
    Y = Rarray[iray, jsurf, RY]
    Z = Rarray[iray, jsurf, RZ]
    print('*** Output: howfar, X,Y,Z:{:4d}{:24.12f}{:24.12f}{:24.12f}'.format(jsurf,X,Y,Z))

def displayXYUV(iray, jsurf):
    X = Rarray[iray, jsurf, RX]
    Y = Rarray[iray, jsurf, RY]
    U = Rarray[iray, jsurf, RU]
    V = Rarray[iray, jsurf, RV]
    print(' X,Y,U,V{:22.14f}{:22.14f}{:22.14f}{:22.14f}'.format(X,Y,U,V))

def displayHXYZ(iray, jsurf):
    X = Rarray[iray, jsurf, RX]
    Y = Rarray[iray, jsurf, RY]
    Z = Rarray[iray, jsurf, RZ]    
    print('*** Output: howfar,  XYZ:{:3d}{:18.12f}{:18.12f}{:20.12f}'.format(jsurf,X,Y,Z))

def displayXYZ(iray, jsurf):
    X = Rarray[iray, jsurf, RX]
    Y = Rarray[iray, jsurf, RY]
    Z = Rarray[iray, jsurf, RZ]
    print(' X,Y,Z: {:24.16f}{:24.16f}{:24.16f}'.format(X,Y,Z))

def displayXY(iray, jsurf):
    X = Rarray[iray, jsurf, RX]
    Y = Rarray[iray, jsurf, RY]
    print(' X, Y: {:24.16f}{:24.16f}'.format(X,Y))

def doMonsterListing():       # List user's ray table results
    prepTableRayGroup(False)  # True=wantListing
    runTableRayGroup(True)    # True=wantListing
    ngood = len(xpoints)
    print("\nNgoodRays = {:6d}".format(ngood) + '  ' + xRMSstring + '  ' + yRMSstring)
    print('Successive ray listing...')
    for iray in range(1, Nrays+1):
        X0 = Rarray[iray, 0, RX]
        U0 = Rarray[iray, 0, RU]
        X1 = Rarray[iray, 1, RX]
        U1 = Rarray[iray, 1, RU]
        X2 = Rarray[iray, 2, RX]
        U2 = Rarray[iray, 2, RU]
        Xf = Rarray[iray, Nsurfs, RX]
        Uf = Rarray[iray, Nsurfs, RU] 
        print('iray, X0, U0, X1, U1, X2, U2, Xf, Uf = {:4d}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:21.15f}{:12.6f}'.format(iray, X0, U0, X1, U1, X2, U2, Xf, Uf))
"""    

  



#-------RAY TRACING INPUT CSV FILE READER--------------

def isEmpty(anyList):
    if len(anyList) < 1:
        return True
    width = 0
    for i in range(0, len(anyList)):
        width = max(width, len(anyList[i]))
    return True if width==0 else False
    
    
#------Wavelength labelling-base=0----------------

FraunhoferIndex      = [   0,     1,     2,     3,     4,     5,     6,     7 ]
FraunhoferLetters    = [  'i',   'g',   'F',   'd',   'R',   'C',   's',   't']
FraunhoferMicrons    = [0.365, 0.436, 0.486, 0.588, 0.638, 0.656, 0.852, 1.014]
FraunhoferNanometers = [  365,   436,   486,   588,   638,   656,   852,  1014]
PLOTCOLORS           = [  'b',   'c',   'g',   'y',   'r',   'c',   'm',   'k']
NWAVES               = 8

def letter2index(letter):
    letter = letter.strip()
    if letter in FraunhoferLetters:
        return FraunhoferLetters.index(letter)
    else:
        print('invalid Fraunhofer letter ' + letter +'   Quitting.')
        quit()
        
def letter2microns(letter):
    index = letter2index(letter)
    return FraunhoferMicrons[index]
    
def letter2nanometers(letter):
    index = letter2index(letter)
    return FraunhoferNanometers[index]
    
def letter2plotColor(letter):
    index = letter2index(letter)
    return PLOTCOLORS[index]
    

#----------FILE GENERAL UNPACKER-----------------------

def unpackCSV(fname):   # returns a 2D rectangular list of all the fields in fname. 
    data = list()       # initially empty.
    # print('\nunpackCSV() Trying: ', fname)
    try:
        data = list(csv.reader(open(fname)))  # 2D list of snippets
    except IOError:
        print("Quitting since could not open "+fname)
        quit()
    if len(data) < 3:
        print("Fewer than three CSV records are found. Quitting this file.")
        return data
    for irow in range(0, len(data)):    
        for jcol in range(0, len(data[irow])):
            data[irow][jcol] = data[irow][jcol].strip()  # unnecessary from Excel
    # print "Initial nrows = ", len(data)

    # delete all empty rows: crucial if empty rows are mixed into the .csv file
    initialLength = len(data)
    for irow in range(initialLength-1, 0, -1):  # don't eliminate title row even if empty
        if isEmpty(data[irow]):
            del data[irow]
    # print "After removing empties, nrows = ", len(data)
    if len(data) < 1:
        print("Nothing left to return. Quitting this file.")
        return data

    # get number of fields from longest row that holds any actual data
    nfields = 0
    for irow in range(0, len(data)):
       for jcol in range(0, len(data[irow])):
           if len(data[irow][jcol]) > 0:
               nfields = max(nfields, jcol+1)
    # print "Nfields = ", nfields
    if nfields < 1:
        print("Nfields is < 1.  Quitting this file.")
        return data

    # make all rows have nfields by appending empty fields where needed.
    for irow in range(len(data)):
        data[irow] = data[irow][:nfields]  # truncate beyond nfields
        while len(data[irow]) < nfields:   # append empty fields
            data[irow].append("")

    return data
    

#-------RAY TRACING SPECIFIC FILE TYPE UNPACKERS----------

def getOpticsCSV(optname):
    # puts user CSV data into a global list "Odata" 
    global Onfields, Nsurfs, Oarray, OglassNames, OneedsMedia, OhasAnySags, Oheaders, Ojmirror, Ojfocal
    # print(optname)
    if len(optname) < 4:
        print("Optics table was not found, but is mandatory.  Quitting.")
        quit()
    data = unpackCSV(optname)    
    
    nrows = len(data)
    ncols = len(data[0])
    # print('getOpticsCSV() has nrows, ncols = ', nrows, ncols)
    # for item in data:
    #     print(item)
    if nrows < 3 or ncols < 2:
        print(myMedFileName, " has returned no usable data. Quitting.")
        quit()    
    
    Onfields = ncols
    Nsurfs = nrows - 2
    guideNumber = suckInt(data[0][0])
    # print "Optics guideNumber = ", guideNumber
    if guideNumber>0:
        Nsurfs = min(Nsurfs, guideNumber)
    Ojfocal = Nsurfs
    # print("  getOpticsCSV() is setting Ojfocal = Nsurfs = ", Ojfocal)
    Oheaders = data[1]
    if Onfields<1 or Nsurfs < 1:
        print(optname, ' has no data available.  Quittimg.')
        quit()
    Oarray = np.zeros([Nsurfs+1, OFINAL])        # rebuild host Oarray
    Oarray.fill(-0.0)
    
    #---set up complete empty lists------------
    OneedsMedia = False
    del OglassNames[:]
    OglassNames.append("base=1")
    for k in range(1, Nsurfs+1):
        OglassNames.append("")

    #-----set literal and numerical OpticsData() field by field--------
    for ifield in range(0, Onfields):
        header = data[1][ifield]
        attrib = getOpticsAttribute(header)
        if attrib == OINDEX:                      # literal not numerical
            for jsurf in range(1, Nsurfs+1):      # jsurf=1, 2, ...Nsurfs
                snippet = data[jsurf+1][ifield]
                OglassNames[jsurf] = snippet 
                if len(snippet) > 0:
                    try:
                        x = float(snippet)        # if not numeric,
                    except:                       # will need a .MED lookup table.
                        OneedsMedia = True
        elif attrib == OACTIONTYPE:
            for jsurf in range(1, Nsurfs+1):
                snippet = data[jsurf+1][ifield]
                iaction = getActionType(snippet)  # returns an action code number
                Oarray[jsurf, OACTIONTYPE] = iaction
                if iaction == OMIRRORACTION:
                    # print('  getOpticsCSV() is setting Ojmirror = ', jsurf)
                    Ojmirror = jsurf
        elif attrib>=0 and attrib<OMAXINPUT:      # numerical data fields including SAGMULT etc
            # if attrib in (26,27,28,29):           # sag attribute range
            #     print "++++++++Parsing attribute number ", attrib
            for jsurf in range(1, Nsurfs+1):      # jsurf=1, 2, ...Nsurfs
                snippet = data[jsurf+1][ifield]
                x = -0.0
                if len(snippet) > 0:
                    try:
                        x = float(snippet)
                    except ValueError:
                        x = -0.0
                Oarray[jsurf, attrib] = x 

    #---For each surface: any need polynomials?-----
    del OhasPoly[:]                           # empty the host list
    OhasPoly.append("base=1")                 # 1-based to match jsurf
    for jsurf in range(1, Nsurfs+1):          # jsurf = 1, 2, ...Nsurfs
        numPoly = 0                           # poly search
        for index in range(OA1, OA14+1):
            if Oarray[jsurf, index] != 0.0: 
                numPoly += 1   
        if numPoly>0:
            OhasPoly.append(True)             # 1-based list like jsurf
        else:
            OhasPoly.append(False)

    #----evaluate all the Euler matrices----
    setEulers()
    
        
def getRaysCSV(rayname, maxrays):  # sets up Raystarts[], Rarray[] etc
    # Be sure to getOpticsCSV() first! since this needs to know Nrays. ?????
    global Rnfields, Nrays, Rarray, Raystarts, Rheaders, RwaveNames
    # print(rayname)
    if len(rayname) < 4:
        print("  Ray table was not found, but is mandatory.  Quitting.")
        quit()
        
    data = unpackCSV(rayname)
    nrows = len(data)
    ncols = len(data[0])
    # print('getRaysCSV() has nrows, ncols = ', nrows, ncols)
    # for item in data:
    #     print(item)
    if nrows < 3 or ncols < 2:
        print(myMedFileName, " has returned no usable data. Quitting.")
        quit()    
    
    Rnfields = ncols
    Rheaders = data[1]
    Nrays = len(data) - 2 
    Nrays = min(Nrays, maxrays)  # limit the table to get monochromatic ray starts 
    guideNumber = suckInt(data[0][0])
    # print("  Rays guideNumber = ", guideNumber)
    if guideNumber>0:
        Nrays = min(Nrays, guideNumber)

    if Rnfields<1 or Nrays < 1:
        print(rayname, ' has no data available.  Quitting.')
        quit()
        
    #---set Raystarts[iray,attrib] field by field-------
    # print 'Creating and Setting Raystarts[iray,attrib]'
    Raystarts = np.zeros([Nrays+1, RFINALINPUT+1])   # base=1, base=0
    Raystarts.fill(-0.0)
    del RwaveNames[:]
    del RFtoI[:]
    del RItoF[:]
    RwaveNames.append("base=1")  # numbering 1...Nrays inclusive
    for iray in range(1, Nrays+1):
        RwaveNames.append(" ")  # start with all blanks
    
    for ifield in range(0, Rnfields):
        header = data[1][ifield]
        attrib = getRayStartAttribute(header)
        RFtoI.append(attrib)
        if attrib == RWAVE:                     # this is the "@wave" field
            for iray in range(1, Nrays+1):
                snippet = data[iray+1][ifield]  # for example "i" for 365nm
                RwaveNames[iray] = snippet      # put snippet int RwaveNames.
        if attrib>=0 and attrib<=RFINALINPUT:   # X0,Y0,Z0,U0,V0,W0,XG,YG
            # found a ray start column in ray table.
            for iray in range(1, Nrays+1):
                snippet = data[iray+1][ifield]
                Raystarts[iray, attrib] = getFloatValue(snippet)
                
    # -----Next set up Rarray[]-------            
    Rarray = np.zeros([Nrays+1, Nsurfs+1, RFINAL])
    # print('  Rarray.shape = ', Rarray.shape)  # expect {Nrays+1, Nsurfs+1,  13}
    Rarray.fill(-0.0)                           # -0.0 means empty field
    for iray in range(1, Nrays+1): 
        #----Fill Rarray[0] based on table RayStarts----------
        # printArray('Raystarts ', Raystarts[iray])
        Rarray[iray, 0, RWAVE]  = Raystarts[iray, RWAVE]           
        X = Rarray[iray, 0, RX] = Raystarts[iray, RX]
        Y = Rarray[iray, 0, RY] = Raystarts[iray, RY]
        Z = Rarray[iray, 0, RZ] = Raystarts[iray, RZ]
        U = Rarray[iray, 0, RU] = Raystarts[iray, RU]
        V = Rarray[iray, 0, RV] = Raystarts[iray, RV]
        W = Rarray[iray, 0, RW] = Raystarts[iray, RW]
        W = Rarray[iray, 0, RW] = fixup(U,V,W)
    
    #---evaluate lookup table RItoF[]-----------
    del RItoF[:]
    for attrib in range(RX, RFINALINPUT+1):
        try:
            field = RFtoI.index(attrib)
        except ValueError:
            field = -1
        RItoF.append(field)

    
def getMediaCSV(myMedFileName):
    global Marray, Nglasses, MwaveNames, MglassNames, OglassNames, RwaveNames    
    # print(myMedFileName)
    if len(myMedFileName) < 4:
        print('valid media table was not found.  Quitting.')
        quit()
    data = unpackCSV(myMedFileName)
    nrows = len(data)
    ncols = len(data[0])
    # print('getMediaCSV() has nrows, ncols = ', nrows, ncols)
    # for item in data:
    #    print(item)
    if nrows < 3 or ncols < 2:
        print(myMedFileName, " has returned no usable data. Quitting.")
        quit()
    guideNumber = suckInt(data[0][0])
    # print('guideNumber = ', guideNumber)    
    nrows = nrows - 1   # avoids title row, includes wavelength header row
    Nglasses = min(nrows-1, guideNumber)
    if Nglasses < 1:
        print('getMediaCSV() finds no glasses.  Quitting.')
        quit()    
    # print('now, after min(), Nglasses including dummy = ', Nglasses)
    Mnfields = len(data[1])
    # print('Mnfields from HEADERROW is... ', Mnfields)
    if Mnfields<1 or Nglasses<1:
        print(medianame, " has no data available.  Quitting.")
        quit()
    # print('Mnfields = ', Mnfields)
    nwaves = Mnfields - 1                    # column zero is glassname, not index
    # print('nwaves   = ', nwaves)
    #---get the glass names first. These are in table column zero----------------
    del MglassNames[:]                       # empty the host list
    MglassNames.append("base=1")             # using base=1 so this is a dummy
    for kglass in range(1, Nglasses+1):      # glasses are 1....Nglasses
        row = kglass + 1                     # glass #1 is in row 2. 
        snippet = data[row][0]               # rows are 0....Nglasses
        MglassNames.append(snippet)
    # print('Here are the MEDIA glassnames including dummy: ', MglassNames)    
    
    #-----next, gather the wavelength names from header row 1, columns 1...Nwaves----
    del MwaveNames[:]                        # empty the host list
    MwaveNames.append("base=1")              # base=1
    for ifield in range(1, nwaves+1):       
        wavename = data[1][ifield]   # HEADERROW is row=1
        MwaveNames.append(wavename)        
    # print('Here are the MEDIA wave names including dummy: ', MwaveNames)      
    
    #----Now get the refractive indices data[][] field by field---------
    Marray = np.zeros([Nglasses+1, Mnfields+1])  # both have base=1
    Marray.fill(-0.0)            
    for kglass in range(1, Nglasses+1):      # row by row, base=1
        row = kglass + 1                     # glass #1 is in row 2
        for ifield in range(1, nwaves+1):    # column by column, base=1
            snippet = data[row][ifield]
            x = -0.0
            if len(snippet) > 0:
                try:
                    x = float(snippet)
                except ValueError:
                    x = -0.0
            Marray[kglass, ifield] = x
          
    #-------check to see that all needed .OPT glass names are present-------------
    for iglass in range(1, Nsurfs+1):   # skip base=0
        glassname = OglassNames[iglass]
        if len(glassname)>0 and not glassname in MglassNames:
            print('  Failed to find .OPT glass name among MEDIA: ', glassname, '  Quitting.')
            quit()
            
    #-------check to see that all needed .RAY wavelength names are present-----------
    for i in range(1, Nrays+1):   #skip base=0
        wavename = RwaveNames[i]
        if not wavename in MwaveNames:
            print('  Failed to find .RAY wave name among MEDIA: ', wavename, '  Quitting.')
            quit()
        

  
    
    
    
    
    
#----RAY TRACE METHODS----------------    
#----RAY TRACE METHODS----------------    
#----RAY TRACE METHODS----------------    
    
def modifyADC(adc1, adc2):  
    global Oarray
    # for use with DESI-ADC.OPT: 26 surfaces like DESI-E2E.OPT 
    # Surfaces are numbered base=1: there is no surface zero.
    # CBout row 11 then ADC1 then CBout row 15 
    # CBout row 17 then ADC2 then CBout row 21.
    # Call this only after .OPT table has been loaded & parsed.
    # get no disperion if adc1 = adc2: same roll angle.
    # New: Southern sky max dispersion if adc1=-90, adc2=+90
    # New: Northern sky max dispersion if adc1=+90, adc2=-90
    # these indices are for DESI-ADC.OPT
    Oarray[11, OROLL] = -adc1  # CBout
    Oarray[15, OROLL] = adc1   # CBout
    Oarray[17, OROLL] = -adc2  # Cbout
    Oarray[21, OROLL] = adc2   # CBout
    setEulers()
    

def modifyWaveName(wavename):
    global RwaveNames
    # Replaces RayTable Rwavenames with this user supplied wavename
    # Often these are Fraunhofer letters 'i' ... 't'
    for iray in range(1, Nrays+1):   # ray table is base=1
        RwaveNames[iray] = wavename

def prepForcedRayGroup(wuv):  
    # {w,u,v} = {wavelength number, Uincoming, Vincoming}
    # Steers all table rays that were installed in getRayCSV().
    # Afffects Rarray[] but not Raystarts[]
    # print('prepForcedRayGroup has received Nrays = ', Nrays)
    # printArray('prepForcedRayGroup wuv = ', wuv)
    waveID = int(wuv[0])                         # base=0
    waveName = FraunhoferLetters[waveID]
    modifyWaveName(waveName)
    for iray in range(1, Nrays+1):               # base=1
        U = Rarray[iray, 0, RU] = wuv[1]
        V = Rarray[iray, 0, RV] = wuv[2] 
        W = Rarray[iray, 0, RW]
        # print('UVW = {:12.6f}{:12.6f}{:12.6f}'.format(U,V,W))
        Rarray[iray, 0, RW] = fixup(U,V,W)
            
def runOneRay(iray):
    # THIS IS THE INNER LOOP given a single ray start 1 <= iray <= Nrays
    # uses Rarray[iray, 0] to get starting coordinates, 
    # so be sure to prep Rarray[] and prepForcedRayGroup() before calling this method.
    howfar = 0 
    code = OK
    for jtarget in range(1, Nsurfs+1):                # base=1
        isCBout = OCBOUTACTION == int(Oarray[jtarget, OACTIONTYPE])
        if isCBout:
            vxtovx(iray, jtarget) 
            vxtolab(iray, jtarget)
            howfar = jtarget
            continue
        labtovx(iray, jtarget)
        code = intercept(iray, jtarget)
        if code == OK:
            code = validate(iray, jtarget)
        if code == OK:
            code = redirect(iray, jtarget)
        if code == OK:
            howfar = jtarget
        vxtolab(iray, jtarget) 
        if code != OK:
            break
        testUVW(iray, jtarget)
    return howfar

    
def runAllTableRays():  
    # THIS IS THE OUTER LOOP that builds [Xfinal, Yfinal, Zfinal, waveNum]
    # Return each table ray result [X0, Y0, Xf, Yf, c] global CS-5 coordinates
    # Assumes all tables have been set up.
    xyzcList = []                                      # list to gather outputs
    for iray in range(1, Nrays+1):                     # base=1
        howfar  = runOneRay(iray)                      # calls inner loop
        if howfar == Nsurfs:                           # this is a "good" ray
            xf   = Rarray[iray, Nsurfs, RX]            # gather ray trace result
            yf   = Rarray[iray, Nsurfs, RY]            # gather ray trace result
            zf   = Rarray[iray, Nsurfs, RZ]            # gather ray trace result
            
            letter = RwaveNames[iray]                  # get waveNum, 0..7 inclusive
            cint = letter2index(letter)
            
            cf   = float(cint)                         # float to allow array()
            xyzc = [xf,yf,zf,cf]                       # List of four numbers
            xyzcList.append(xyzc)                      # stack it into output
    nrows = len(xyzcList)
    if nrows < 2:
        print('RT finds fewer than two good rays.  Quitting.')
        quit()
    xyzcArray = np.array((xyzcList))  
    return xyzcArray
    
    








#-------RAY TRACING OUTPUT GENERATORS-------------------------



def doSpotDiagram(xyzcArray, title, u0, v0, adc1, adc2): 
    # This plot ignores zfinal; plots xf,yf with specified color
    # Here, c = float color ID: 1.0, 2.0, ... 8.0
    nrays = len(xyzcArray)
    if nrays < 1:
        print('Yikes, spotDiagram() has received no rays.  Quitting.')
        quit()
    ngood = nrays
    # print('RT184 doSpotDiagram() has received nrays = ', nrays)
    
    xvals = xyzcArray[:,0]   # all rows column 0 = xfinal
    yvals = xyzcArray[:,1]   # all rows column 1 = yfinal
    cvals = xyzcArray[:,3]   # all rows column 3 = wavenums: 0.0, 1.0, 2.0, ...
    
    xave = np.average(xvals)
    yave = np.average(yvals)
    xrms = np.std(xvals)
    yrms = np.std(yvals)
    xmax  = np.amax(xvals)
    xmin  = np.amin(xvals)
    xmid  = 0.5*(xmax+xmin)
    xspan = xmax - xmin
    ymax  = np.amax(yvals)
    ymin  = np.amin(yvals)
    ymid  = 0.5*(ymax + ymin)
    yspan = ymax - ymin
    MINSPAN = 0.1  # 0.1mm or 100um
    EXTRA = 1.1
    span = max(MINSPAN, EXTRA*max(xspan, yspan))
    half = 0.5*span    
    
    colors = list()
    for iray in range(nrays):
        icolor = int(cvals[iray])
        plotcolor = PLOTCOLORS[icolor]
        colors.append(plotcolor)

    fig, ax = plt.subplots(figsize=(6,6))     # yes plural even for one plot
    ax.tick_params(direction='in')          
    ax.scatter(xvals, yvals, c=colors)        # plot the dots
    
    ax.set_xlabel('Xfp, mm, +eastward')
    ax.set_ylabel('Yfp, mm, +southward')
    ax.axis('equal')  # equal x and y display scales, or see span generator below
    # title = 'U0, V0 = {:+12.6f}{:+12.6f}'.format(u0, v0)
    plt.title(title)   

    # print('SpotDiagram scaled span = {:9.3f}'.format(span))
    ax.set_xlim(xmid-half, xmid+half)
    ax.set_ylim(ymid-half, ymid+half) 
    
    for i in range(NWAVES):  # base=0
        item = FraunhoferLetters[i] + '  ' + str(FraunhoferNanometers[i]) + 'nm'
        ax.text(0.86, 0.96-0.03*i, item, transform=ax.transAxes, fontsize=8, color=PLOTCOLORS[i])  

    u0string     = 'U0 = {:9.6f}'.format(u0)
    v0string     = 'V0 = {:9.6f}'.format(v0)
    ngoodstring  = 'Ngood = ' + str(ngood)
    xAVEstring   = 'Xave = {:9.3f}'.format(xave)
    yAVEstring   = 'Yave = {:9.3f}'.format(yave)
    xRMSstring   = 'Xrms = {:9.3f}'.format(xrms)
    yRMSstring   = 'Yrms = {:9.3f}'.format(yrms)
    intspanum    = int(1000*span)
    spanstring   = 'PlotSpan = '+str(intspanum) + r'$\mu$m'

    ax.text(0.02, 0.96, ngoodstring,  transform=ax.transAxes, fontsize=8) 
    ax.text(0.02, 0.93, xAVEstring,   transform=ax.transAxes, fontsize=8)  
    ax.text(0.02, 0.90, yAVEstring,   transform=ax.transAxes, fontsize=8) 
    ax.text(0.02, 0.87, xRMSstring,   transform=ax.transAxes, fontsize=8)  
    ax.text(0.02, 0.84, yRMSstring,   transform=ax.transAxes, fontsize=8) 
    ax.text(0.02, 0.81, spanstring,   transform=ax.transAxes, fontsize=8) 
     
    ax.text(0.87, 0.02, PROGNAME, transform=ax.transAxes, fontsize=8)  
    fig.tight_layout()
    
    # CAUTION LATEX FIGURE FILENAMES FORBIDS MOST PUNCTUATION AND SPACES
    # period is OK before suffix but nowhere else: NO DECIMAL POINTS darn it. 
    # However, plus, minus, equals, underscores are OK.
    # To eliminate this trouble I use tiny unit integers: no decimal points needed.  

    ustring   = '{:+10.0f}'.format(1000000*u0).strip()   # micro radians
    vstring   = '{:+10.0f}'.format(1000000*v0).strip()   # micro radians
    figfilename = title + '.png' 
    print('Saving... ' +figfilename)
    fig.savefig(figfilename, dpi=300)
    plt.close(fig)
    
    
"""
def showStatistics(xyzcArray, u0, v0):    
    ngood = len(xyzcArray)
    xvals = xyzcArray[:,0]   # all rows column 0 = xfinal
    yvals = xyzcArray[:,1]   # all rows column 1 = yfinal
    zvals = xyzcArray[:,2]   # all rows column 2 = zfinal
    xave = np.average(xvals)
    yave = np.average(yvals)
    zave = np.average(zvals)
    xrms = np.std(xvals)
    yrms = np.std(yvals)
    zrms = np.std(zvals)
    result = '{:+10.6f}{:+10.6f}{:8d}{:+10.3f}{:+10.3f}{:+10.3f}{:+10.3f}{:+10.3f}{:+10.3f}' \
        .format(  u0,      v0,   ngood,  xave,   yave,    zave,   xrms,    yrms,    zrms)
    print(result) 
  
def getSevenStatistics(xyzcArray):    
    print('getStatistics has received xyzcArray.shape = ', xyzcArray.shape) 
    ngood = 1.0*len(xyzcArray)
    xvals = xyzcArray[:,0]   # all rows column 0 = xfinal
    yvals = xyzcArray[:,1]   # all rows column 1 = yfinal
    zvals = xyzcArray[:,2]   # all rows column 2 = zfinal
    print('getStatistics has extracted xvals.shape = ', xvals.shape)
    print('these xvals are...', xvals)
    xave = np.average(xvals)
    yave = np.average(yvals)
    zave = np.average(zvals)
    xrms = np.std(xvals)
    yrms = np.std(yvals)
    zrms = np.std(zvals)
    result = [ngood, xave, yave, zave, xrms, yrms, zrms]
    return np.array((result))
"""    
     
     
     
     
     
     
     
     
     
#------FUNCTION FOR EXTERNAL (END-TO-END) CALLS in e2e-25.py ----------------


def getNine(wuv12s):
    global Nsurfs, Nrays, Nglasses
    # for one star, one wuv, specified wavelength
    # Monochromatic: just one row of {waveNo, U0, V0, adc1, adc2}
    # Polychromatic: 8 rows of 5 columns
    
    #--------set up the optics, rays, and media files--------------
    myOptFileName = 'DESI-ADC-2.OPT.CSV'     # Hexapod roll row 4; ADC rolls CBouts at row 11,15,17,21
    myRayFileName = 'DESI-ADC-2.RAY.CSV'     # on axis; 84 rays per pupil, only 84 rays used.
    myMedFileName = 'DESI-ADC-2.MED.CSV'     # 8 glasses, 8 wavels' i'...'t'
    # print("Loading files: " + myOptFileName + '  ' + myRayFileName + '  ' + myMedFileName)
    getOpticsCSV(myOptFileName)            # set up optical prescription
    if Nsurfs < 1:
       print('  Nsurfs < 1; quitting.')
       quit()
    getRaysCSV(myRayFileName, 84)          # set up Raystarts[] and Rarray[]; 84 forces monochromatic
    if Nrays < 1:
        print('  Nrays < 1; quitting.')
        quit()
    getMediaCSV(myMedFileName)             # set up internal glass refraction table
    if Nglasses < 1:
        print('  Nglasses < 1; quitting.')
        quit()
        
    # Next study our wuv12s: monochromatic, or polychromatic?
    # print('RT184 has received these wuv12s: \n', wuv12s)    
    size = wuv12s.size
    nwavels = size // 5
    print('RT185:getNine is given wuv12s.size =   ', size)
    print('so RT185:getNine has adoptec nwavels = ', nwavels)
    if size < 5:
        print('RT185:getNine quitting; wuv12s.size = ', size)
        quit()
    if size == 5:   # enlarge to 2D array for indexing
        wuv12s = np.vstack((wuv12s, np.array([0.,0.,0.,0.,0.])))
    adc1 = wuv12s[0,3]                    # get the 4th element
    adc2 = wuv12s[0,4]                    # get the 5th element
    # print('RT185:getNine is passing the two ADC angles: ', adc1, adc2)
    modifyADC(adc1, adc2)                 # use them
    
    # Next work through the specified wavelengths
    xlist = []
    ylist = []
    zlist = []
    ngood = 0.0
    for iwave in range(nwavels):
        wuv  = wuv12s[iwave,:3]            # first three elements 0, 1, 2        
        prepForcedRayGroup(wuv)            # use these {wavels, U0, V0}    
        xyzcArray = runAllTableRays()      # run this wavelength on this target; get 84rows, 4cols

        #-------statistics() calls here-----------
        # print('getNine has received xyzcArray.shape = ', xyzcArray.shape) 
        ngood += 1.0*len(xyzcArray)
        xlist.append(xyzcArray[:,0])   # all rows column 0 = xfinal
        ylist.append(xyzcArray[:,1])   # all rows column 1 = yfinal
        zlist.append(xyzcArray[:,2])   # all rows column 2 = zfinal
    print('RT185:getNine has extracted len(xlist) = ', len(xlist))
    xave = np.average(xlist)
    yave = np.average(ylist)
    zave = np.average(zlist)
    xrms = np.std(xlist)
    yrms = np.std(ylist)
    zrms = np.std(zlist)
    resultstr = '{:9.3f}{:9.3f}{:6.0f}{:9.3f}{:9.3f}{:9.3f}{:9.3f}{:9.3f}{:9.3f}'.format(adc1, adc2, ngood, xave, yave, zave, xrms, yrms, zrms)
    print('RT185:getNine: ', resultstr)
    resultArray = np.array([adc1, adc2, ngood, xave, yave, zave, xrms, yrms, zrms])
    return resultArray

    
    
"""
def getSeven(starnumber, wuv12s):
    ndim = np.asarray(wuv12s).ndim
    if ndim == 0:
        print("Quitting, no rays. ")
        quit()
    if ndim == 1:
        return getSevenMono(starnumber, wuv12s)
    else:
        return getSevenPoly(starnumber, wuv12s)
        

def getSevenPoly(starnumber, wuv12s):   
    global Nsurfs, Nrays, Nglasses
    # star number is just an output key
    # Polychromatic: Each row of wuv12s array is:  {waveNo, U0, V0, adc1, adc2}
    # so, 8 rows of 5 columns
    nwaves = len(wuv12s)                  # 8 rows   
    print('nwaves = ', nwaves)
    if nwaves < 1:
        print('getResultsOneStar(): no rays. Quitting.')
        quit()
    #--------set up the optics, rays, and media files--------------
    myOptFileName = 'DESI-ADC.OPT.CSV'     # Hexapod roll row 4; ADC rolls CBouts at row 11,15,17,21
    myRayFileName = 'DESI-ADC.RAY.CSV'     # on axis; 84 rays per pupil, only 84 rays used.
    myMedFileName = 'DESI-ADC.MED.CSV'     # 8 glasses, 8 wavels' i'...'t'
    # print("Loading files: " + myOptFileName + '  ' + myRayFileName + '  ' + myMedFileName)
    getOpticsCSV(myOptFileName)            # set up optical prescription
    if Nsurfs < 1:
       print('  Nsurfs < 1; quitting.')
       quit()
    getRaysCSV(myRayFileName, 84)          # set up Raystarts[] and Rarray[]; 84 forces monochromatic
    if Nrays < 1:
        print('  Nrays < 1; quitting.')
        quit()
    getMediaCSV(myMedFileName)             # set up internal glass refraction table
    if Nglasses < 1:
        print('  Nglasses < 1; quitting.')
        quit()
    # now specialize the table for 8 monochromatic runs with 8 {u,v} pairs, one per wavelength
    # Trick here is to concatenate eight 84x4 subLists into a bigList

    xyzcBigList = []                       # Start with empty polychromatic list
    for iwave in range(NWAVES):            # trace each wavelength from dispersed sky; base=0
        wuv12 = wuv12s[iwave]              # {waveno, U0, V0} triplet; wuvs are base=1
        adc1 = wuv12[3]                    # get the 4th element
        adc2 = wuv12[4]                    # get the 5th element
        modifyADC(adc1, adc2)              # use them
        wuv  = wuv12[:3]                   # first three elements 0, 1, 2        
        prepForcedRayGroup(wuv)            # use these {wavels, U0, V0}    
    
        xyzcArray = runAllTableRays()      # run this wavelength on this target; get 84rows, 4cols
        xyzcSubList = xyzcArray.tolist()   
        xyzcBigList += xyzcSubList         # extending not appending: makes the list taller.
    xyzcArray = np.array(xyzcBigList)      # convert to an array
    #--------SpotDiagram() and getStatistics() calls here-----------
    adc1str = '{:+6.0f}'.format(adc1).strip()
    adc2str = '{:+6.0f}'.format(adc2).strip()
    u0 = wuv12s[4][1]   # R band wavelength
    v0 = wuv12s[4][2]   # R band wavelength
    title = PROGNAME+'_'+str(starnumber) + '_' + adc1str+'_'+adc2str
    print('RT184 is writing.... ' + title + '.png')
    doSpotDiagram(xyzcArray, title, u0, v0, adc1, adc2)
    return getSevenStatistics(xyzcArray)


def getSevenMono(starnumber, wuv12):   
    global Nsurfs, Nrays, Nglasses
    # star number is just an output key
    # Monochromatic: just one row of {waveNo, U0, V0, adc1, adc2}
    # so, 8 rows of 5 columns
    printArray('RT184 has received this wuv12: " ', wuv12)

    #--------set up the optics, rays, and media files--------------
    myOptFileName = 'DESI-ADC.OPT.CSV'     # Hexapod roll row 4; ADC rolls CBouts at row 11,15,17,21
    myRayFileName = 'DESI-ADC.RAY.CSV'     # on axis; 84 rays per pupil, only 84 rays used.
    myMedFileName = 'DESI-ADC.MED.CSV'     # 8 glasses, 8 wavels' i'...'t'
    # print("Loading files: " + myOptFileName + '  ' + myRayFileName + '  ' + myMedFileName)
    getOpticsCSV(myOptFileName)            # set up optical prescription
    if Nsurfs < 1:
       print('  Nsurfs < 1; quitting.')
       quit()
    getRaysCSV(myRayFileName, 84)          # set up Raystarts[] and Rarray[]; 84 forces monochromatic
    if Nrays < 1:
        print('  Nrays < 1; quitting.')
        quit()
    getMediaCSV(myMedFileName)             # set up internal glass refraction table
    if Nglasses < 1:
        print('  Nglasses < 1; quitting.')
        quit()
    adc1 = wuv12[3]                    # get the 4th element
    adc2 = wuv12[4]                    # get the 5th element
    modifyADC(adc1, adc2)              # use them
    wuv  = wuv12[:3]                   # first three elements 0, 1, 2        
    prepForcedRayGroup(wuv)            # use these {wavels, U0, V0}    
    xyzcArray = runAllTableRays()      # run this wavelength on this target; get 84rows, 4cols
    #--------SpotDiagram() and getStatistics() calls here-----------
    adc1str = '{:+6.0f}'.format(adc1).strip()
    adc2str = '{:+6.0f}'.format(adc2).strip()
    u0 = wuv12[1]
    v0 = wuv12[2]
    title = PROGNAME+'_'+str(starnumber) + '_' + adc1str+'_'+adc2str
    print('RT184 is writing.... ' + title + '.png')
    doSpotDiagram(xyzcArray, title, u0, v0, adc1, adc2)
    return getSevenStatistics(xyzcArray)
"""
