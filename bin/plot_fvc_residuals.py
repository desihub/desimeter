#!/usr/bin/env python

import argparse
import os.path
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

ii = np.where((table["LOCATION"]>0)&(table["XMETRO"]!=0))[0]
x = table["XFP"][ii]
y = table["YFP"][ii]
xm = table["XMETRO"][ii]
ym = table["YMETRO"][ii]

a = plt.subplot(1,1,1)
a.set_title(os.path.basename(filename))

plt.plot(x,y,".")
plt.plot(xm,ym,".")

for i in range(x.size) :
    plt.plot([x[i],xm[i]],[y[i],ym[i]],"-",c="gray",alpha=0.5)

dx=xm-x
dy=ym-y
plt.quiver(x,y,dx,dy)
plt.xlabel("XFP (mm)")
plt.ylabel("YFP (mm)")
plt.show()

