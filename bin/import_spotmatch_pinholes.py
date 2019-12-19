#!/usr/bin/env python


"""
Code used (once) to import spotmatch pinhole references.
"""


import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from pkg_resources import resource_filename

filename = resource_filename('desimeter',"data/fp-metrology.csv")
spots    = Table.read(filename,format="csv")
spots_fidid = spots["Petal ID"]*1000+spots["Device Loc ID"]
spots_identifier = np.array(spots_fidid)*100 + np.array(spots["DOTID"])

print(spots_identifier)


#echo "FID_ID,PIN_ID,XPIX,YPIX" > tmp.csv
#cat spotmatch-pinhole-references.dat | grep pinhole | awk '{print $1","$5","$2","$3}' >> tmp.csv

pinholes_filename = resource_filename('desimeter',"tmp.csv")
pinholes = Table.read(pinholes_filename,format="csv")

print(pinholes.dtype.names)

pinholes_dotid = np.zeros(pinholes["FID_ID"].size)
pinholes_dotid[pinholes["PIN_ID"]==259]  = 1
pinholes_dotid[pinholes["PIN_ID"]==515]  = 2
pinholes_dotid[pinholes["PIN_ID"]==771]  = 3
pinholes_dotid[pinholes["PIN_ID"]==1027] = 4
pinholes_identifier = np.array(pinholes["FID_ID"])*100 + pinholes_dotid


dico_spots={f:index for index,f in enumerate(spots_identifier)}

indices_pinholes=[]
indices_spots=[]

for i,f in enumerate(pinholes_identifier) :
    if f in dico_spots :
        indices_spots.append(dico_spots[f])
        indices_pinholes.append(i)



newpinholes = spots[:][indices_spots]
newpinholes["XPIX"] = pinholes["XPIX"][indices_pinholes]
newpinholes["YPIX"] = pinholes["YPIX"][indices_pinholes]

newpinholes.write("pinholes-fvxy.csv",format="csv",overwrite=True)

    
