#!/usr/bin/env python


import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from pkg_resources import resource_filename

filename = resource_filename('desicoord',"data/fp-metrology.csv")
spots    = Table.read(filename,format="csv")

print(spots.dtype.names)
print(np.unique(spots["Device Type"]))

for petal in np.unique(spots['Petal Loc ID']) :
    selection=(spots['Petal Loc ID']==petal)&(spots['Device Type']=="POS")
    petalid = spots['Petal ID'][selection][0]
    plt.plot(spots["XFP"][selection],spots["YFP"][selection],".",alpha=0.5,label="PETAL LOC={} ID={}".format(petal,petalid))

selection=(spots['Device Type']=="FIF")
plt.plot(spots["XFP"][selection],spots["YFP"][selection],".",label="FIF")
selection=(spots['Device Type']=="GIF")
plt.plot(spots["XFP"][selection],spots["YFP"][selection],".",label="GIF")
    
plt.legend()
plt.show()


