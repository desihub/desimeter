#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import desimeter.transform.zhaoburge

x1d=np.linspace(-400,400,11)
x=np.tile(x1d,(x1d.size,1))
y=x.copy().T

for p in range(8) :
    a = plt.subplot(2,4,p+1)
    xp,yp,label = desimeter.transform.zhaoburge.getZhaoBurgeTerm(p, x, y)
    a.set_title(label)
    a.quiver(x,y,xp-x,yp-y,label=label)
    
plt.show()
