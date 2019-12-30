#!/usr/bin/env python


import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from pkg_resources import resource_filename

filename = resource_filename('desimeter',"data/fp-metrology.csv")
spots    = Table.read(filename,format="csv")

print(spots.dtype.names)
print(np.unique(spots["DEVICE_TYPE"]))

plt.figure(figsize=(6,6))
for petal in np.unique(spots['PETAL_LOC']) :
    selection=(spots['PETAL_LOC']==petal)&(spots['DEVICE_TYPE']=="POS")
    petalid = spots['PETAL_ID'][selection][0]
    plt.plot(spots["X_FP"][selection],spots["Y_FP"][selection],".",alpha=0.5,label="PETAL LOC={} ID={}".format(petal,petalid))

selection=(spots['DEVICE_TYPE']=="FIF")
plt.plot(spots["X_FP"][selection],spots["Y_FP"][selection],".",label="FIF")
selection=(spots['DEVICE_TYPE']=="GIF")
plt.plot(spots["X_FP"][selection],spots["Y_FP"][selection],".",label="GIF")
    
plt.legend()
plt.show()


