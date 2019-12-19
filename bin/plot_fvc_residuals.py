#!/usr/bin/env python

import argparse
import os.path
import sys
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from pkg_resources import resource_filename

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


# plotting all of FIF and GIF
filename = resource_filename('desimeter',"data/fp-metrology.csv")
metrology = Table.read(filename,format="csv")
selection=(metrology["Device Type"]=="FIF")|(metrology["Device Type"]=="GIF")
plt.plot(metrology["XFP"][selection],metrology["YFP"][selection],"o",c="gray")

plt.plot(x,y,"o")
plt.plot(xm,ym,"X")

dx=xm-x
dy=ym-y

jj=np.where((np.abs(dx)>0.1)|(np.abs(dy)>0.1))[0]
if jj.size > 0 : # sadly
    print("Large residuals")
    for j in jj :
        i=ii[j]
        line = "dx={:4.3f}mm dy={:4.3f}mm".format(dx[j],dy[j])
        for k in table.dtype.names :
            if k=="XERR" or k=="YERR" : continue
            line += " {}={:4.3f}".format(k,table[k][i])
        print(line)
        plt.plot(x[j],y[j],"+",color="red")


plt.quiver(x,y,dx,dy)
plt.xlabel("XFP (mm)")
plt.ylabel("YFP (mm)")
plt.show()

