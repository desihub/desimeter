import numpy as np
from scipy.spatial import cKDTree as KDTree
from desimeter.match import match_same_system
from desimeter.simplecorr import SimpleCorr

def average_coordinates(tables,xkey,ykey) :
    """
    Average x,y coordinates given by xkey and ykey from a list of astropy tables
    and return an astropy table with the average coordinates.
    This function includes a match and a transformation per table.
    The tables do not necessarily have the same number of entries (in case of false detections).
    Args
      tables: list of astropy.Table objects with same columns
    Returns astropy.Table with same columns as input
    """
    table1=None
    x1=None
    y1=None
    indices=None
    xx=[]
    yy=[]
    for table in tables :
        x2=np.array(table[xkey])
        y2=np.array(table[ykey])
        if x1 is None :
            table1=table
            x1=x2
            y1=y2
            indices=np.arange(len(x1),dtype=int)
        else :
            # match the two sets of spots
            indices_1 = np.arange(len(x1),dtype=int)
            indices_2, distances = match_same_system(x1,y1,x2,y2)
            ok=np.where((indices_2>=0)&(distances<5.))[0]
            indices_1 = indices_1[ok]
            indices_2 = indices_2[ok]
            distances = distances[ok]
            x1=x1[indices_1]
            y1=y1[indices_1]
            indices=indices[indices_1]
            x2=x2[indices_2]
            y2=y2[indices_2]
            for i in range(len(xx)) :
                xx[i]=xx[i][indices_1]
                yy[i]=yy[i][indices_1]
            # adjust a possible transfo between the two FVC images
            corr=SimpleCorr()
            corr.fit(x2,y2,x1,y1)
            x2,y2=corr.apply(x2,y2)

        xx.append(x2)
        yy.append(y2)
    xx=np.vstack(xx)
    yy=np.vstack(yy)

    xrms=np.std(xx,axis=0)
    yrms=np.std(yy,axis=0)
    mx=np.mean(xx,axis=0)
    my=np.mean(yy,axis=0)
    table1[xkey][indices]=mx
    table1[ykey][indices]=my
    
    print("number of entries found in all tables= {}".format(indices.size))
    print("rms({})= {:4.3f}".format(xkey,np.median(xrms)))
    print("rms({})= {:4.3f}".format(ykey,np.median(yrms)))

    return table1
