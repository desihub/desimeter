#!/usr/bin/env python


import argparse
import sys
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Plot FVC spots""")
parser.add_argument('-i','--infile', type = str, nargs = "*", default = None, required = True,
                    help = 'path to one or several FVC spots tables in CSV format')


args  = parser.parse_args()

timeindex=dict()
xpix=dict()
ypix=dict()

for number,filename in enumerate(args.infile) :
    
    table=Table.read(filename,format="csv")
    
    if "XPIX" in table.dtype.names :
        x=table["XPIX"]
        y=table["YPIX"]
#    elif "Q" in  table.dtype.names :
#        x=table["S"]*np.cos(table["Q"]*np.pi/180)
#        y=table["S"]*np.sin(table["Q"]*np.pi/180)
    else :
        print("don't know what to do")
        sys.exit(12)
    
    ii = np.where(table["LOCATION"]>0)[0]
    plt.plot(x[ii],y[ii],".")
    
   
    for i in ii :
        loc=table["LOCATION"][i]
        if int(loc) not in timeindex.keys() :
            timeindex[loc] = [ number ,]
            xpix[loc] = [ x[i], ]
            ypix[loc] = [ y[i], ]
        else :
            timeindex[loc].append(number)
            xpix[loc].append(x[i])
            ypix[loc].append(y[i])
    
    
if "DOTID" in table.dtype.names and number == 0 :
    ii=table["DOTID"]==4
    if np.sum(ii)>0 :
        plt.plot(x[ii],y[ii],"X")

if len(args.infile)>1 :
    plt.figure("XPIX")
    for loc in xpix.keys() :
        plt.plot(timeindex[loc],xpix[loc]-np.mean(xpix[loc]),"o-")
    plt.figure("YPIX")
    for loc in ypix.keys() :
        plt.plot(timeindex[loc],ypix[loc]-np.mean(ypix[loc]),"o-")
    

plt.show()

