#!/usr/bin/env python

import argparse
import sys
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Plot FVC spots, showing residuals with metrology""")
parser.add_argument('-i','--infile', type = str, default = None, required = True,
                    help = 'path to a FVC spots table in CSV format (with XFP,YFP,XMETRO,YMETRO columns)')

args  = parser.parse_args()

timeindex=dict()
xpix=dict()
ypix=dict()

filename = args.infile
    
table=Table.read(filename,format="csv")

ii = np.where(table["LOCATION"]>0)[0]
x = table["XFP"][ii]
y = table["YFP"][ii]
xm = table["XMETRO"][ii]
ym = table["YMETRO"][ii]
    
plt.plot(x,y,".")
dx=xm-x
dy=ym-y
plt.quiver(x,y,dx,dy)
plt.show()

