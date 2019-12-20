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

fig = plt.figure()
a = plt.subplot(1,1,1)
a.set_title(os.path.basename(filename))

a.plot(table["XFP"],table["YFP"],".",alpha=0.5,label="all spots")

# plotting all of FIF and GIF
filename = resource_filename('desimeter',"data/fp-metrology.csv")
metrology = Table.read(filename,format="csv")
selection=(metrology["Device Type"]=="FIF")|(metrology["Device Type"]=="GIF")
a.scatter(metrology["XFP"][selection],metrology["YFP"][selection],marker="o",edgecolors="gray",alpha=1.,facecolors="none",label="all FIF and GIF metrology")

a.plot(x,y,".",color="purple",label="matched measured spots")
a.scatter(xm,ym,marker="o",edgecolors="orange",facecolors="none",label="matched metrology")

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
        label=None
        if j==jj[0]: label="large residual"
        a.plot(x[j],y[j],"+",color="red",markersize=12,label=label)


a.quiver(x,y,dx,dy)
a.set_xlabel("XFP (mm)")
a.set_ylabel("YFP (mm)")
a.legend(loc="upper left")

dist=np.sqrt(dx**2+dy**2)
a2=fig.add_subplot(666)
a2.hist(dist[dist<0.1]*1000.)
a2.set_xlabel("dist. (um)")


if False :
    print("Write coordinates of missing or incorrect fiducials")
    x=table["XFP"]
    y=table["YFP"]

    xc=-289
    yc=217
    device_type="GIF"
    petal_id=11
    petal_loc=6
    device_loc=999 # made up
    location=1000*petal_loc+device_loc
    print(table.dtype.names)
    
    # find nearests
    ii=np.where((x-xc)**2+(y-yc)**2<5**2)[0]
    #col0,Petal ID,Petal Loc ID,Device Loc ID,Device Type,DOTID,XPL,YPL,ZPL,X MNT,Y MNT,Z MNT,FCL Source,Proj Distance,Provenance,Notes,XFP,YFP,ZFP,LOCATION
    j=np.argmin((x-xc)**2+(y-yc)**2)
    line=10000
    c=1
    for i in ii :
        if i==j :
            dot=4 # center
        else :
            dot=c
            c+=1
            if c==4 : c+= 1
            
        print(x[i],y[i])
        print("{line},{petal_id},{petal_loc},{device_loc},{device_type},{dot},0,0,0,0,0,0,'interpolated',0,'none','none',{x},{y},0,{location}".format(line=line,petal_id=petal_id,petal_loc=petal_loc,device_loc=device_loc,device_type=device_type,dot=dot,x=x[i],y=y[i],location=location))
        line += 1

    xc=194.7
    yc=248.4
    # find nearests
    ii=np.where((x-xc)**2+(y-yc)**2<2**2)[0]
    for i in ii :
        print(table[:][i])
    
    device_type="FIF"
    petal_id=8
    petal_loc=4
    device_loc=321 # made up
    location=4321
    dot=5 # new
    print(table.dtype.names)
    j=np.argmin((x-xc)**2+(y-yc)**2)
    print("{line},{petal_id},{petal_loc},{device_loc},{device_type},{dot},0,0,0,0,0,0,'interpolated',0,'patch','none',{x},{y},0,{location}".format(line=line,petal_id=petal_id,petal_loc=petal_loc,device_loc=device_loc,device_type=device_type,dot=dot,x=x[j],y=y[j],location=location))
    line += 1


    # replacements
    for (xc,yc) in [ (-47.4,41.7) , ( 367.6,176.7) ] :
        
        ii=np.where((x-xc)**2+(y-yc)**2<2**2)[0]

        for i in ii :
            selection=(metrology["LOCATION"]==table["LOCATION"][i])
            for k in ["DOTID"] :
                selection &= (metrology[k]==table[k][i])
                jj=np.where(selection)[0]
                if jj.size == 0 :
                    print("no match???????")
                    continue
                j=jj[0]
                print("X {} -> {}".format(metrology["XFP"][j],table["XFP"][i]))
                print("Y {} -> {}".format(metrology["YFP"][j],table["YFP"][i]))
                metrology["XFP"][j] = table["XFP"][i]
                metrology["YFP"][j] = table["YFP"][i]
                line="{} : ".format(i)
                for k in metrology.dtype.names :
                    line += ",{}".format(metrology[k][j])
            print(line)
    
    

plt.show()


