"""
Utility functions to match two catalogs of 2D coordinates
"""


import numpy as np
from scipy.spatial import cKDTree as KDTree
from desimeter.log import get_logger


def compute_triangles_with_arbitrary_orientation(x,y) :
    """
    triangle vertices are ordered using the side length
     so it is independent of a possible rotation
    """
    
    tk=[] # indices 
    tr=[] # max side length ratio
    tc=[] # cosine of first vertex (after sorting according to side length)
    ts=[] # orientation
    
    nn=len(x)
    for i in range(nn) :
        for j in range(i+1,nn) :
            for k in range(j+1,nn) :
                # x y of vertices
                ijk=np.array([i,j,k])
                tx=x[ijk]
                ty=y[ijk]

                # sorting according to length (square)
                # is not working well for our case because
                # a lot of triangles are isosceles
                tl2=np.array([(tx[1]-tx[0])**2+(ty[1]-ty[0])**2,(tx[2]-tx[1])**2+(ty[2]-ty[1])**2,(tx[0]-tx[2])**2+(ty[0]-ty[2])**2])
                pairs=np.array([[0,1],[1,2],[0,2]])
                ii=np.argsort(tl2)
                ordering = np.zeros(3).astype(int)
                ordering[0] = np.intersect1d(pairs[ii[0]],pairs[ii[2]]) # vertex connected to shortest and longest side  
                ordering[1] = np.intersect1d(pairs[ii[0]],pairs[ii[1]]) # vertex connected to shortest and intermediate side  
                ordering[2] = np.intersect1d(pairs[ii[1]],pairs[ii[2]]) # vertex connected to intermediate and longest side
                
                ijk=ijk[ordering]
                tx=tx[ordering]
                ty=ty[ordering]
                
                length2=np.array([(tx[1]-tx[0])**2+(ty[1]-ty[0])**2,(tx[2]-tx[1])**2+(ty[2]-ty[1])**2,(tx[0]-tx[2])**2+(ty[0]-ty[2])**2])
                r=np.sqrt(np.max(length2)/np.min(length2)) # ratio of longest to shortest side
                
                c=((tx[1]-tx[0])*(tx[2]-tx[0])+(ty[1]-ty[0])*(ty[2]-ty[0]))/np.sqrt( ((tx[1]-tx[0])**2+(ty[1]-ty[0])**2)*((tx[2]-tx[0])**2+(ty[2]-ty[0])**2)) # cos of angle of first vertex

                s=((tx[1]-tx[0])*(ty[2]-ty[0])-(tx[2]-tx[0])*(ty[1]-ty[0]))/np.sqrt( ((tx[1]-tx[0])**2+(ty[1]-ty[0])**2)*((tx[2]-tx[0])**2+(ty[2]-ty[0])**2)) # orientation whether vertices are traversed in a clockwise or counterclock-wise sense

                                
                tk.append(ijk)
                tr.append(r)
                tc.append(c)
                ts.append(s)
    return np.array(tk),np.array(tr),np.array(tc),np.array(ts)



def compute_triangles_with_fixed_orientation(x,y) :
    """
    Triangle vertices are ordered using one of the coordinates.
    This makes this algorithm sensitive to the orientation of the coordinate frame,
    but more robust if we do not expect large rotations.
    It would however fail for situations where we don't know the rotation
    """
    
    tk=[] # indices 
    tr=[] # max side length ratio
    tc=[] # cosine of first vertex (after sorting according to side length)
    
    nn=len(x)
    for i in range(nn) :
        for j in range(i+1,nn) :
            for k in range(j+1,nn) :
                # x y of vertices
                ijk=np.array([i,j,k])
                tx=x[ijk]
                ty=y[ijk]

                # sorting according to length (square)
                # is not working well for our case because
                # a lot of triangles are isosceles
                #tl2=np.array([(tx[1]-tx[0])**2+(ty[1]-ty[0])**2,(tx[2]-tx[1])**2+(ty[2]-ty[1])**2,(tx[0]-tx[2])**2+(ty[0]-ty[2])**2])
                #pairs=np.array([[0,1],[1,2],[0,2]])
                
                #ii=np.argsort(tl2)
                #ordering = np.zeros(3).astype(int)
                #ordering[0] = np.intersect1d(pairs[ii[0]],pairs[ii[2]]) # vertex connected to shortest and longest side  
                #ordering[1] = np.intersect1d(pairs[ii[0]],pairs[ii[1]]) # vertex connected to shortest and intermediate side  
                #ordering[2] = np.intersect1d(pairs[ii[1]],pairs[ii[2]]) # vertex connected to intermediate and longest side
                
                # I sort the vertices with one arbitrary axis, x, because the rotation
                # of the field is not large enough to alter this
                # (in fact I could choose a different ordering per triangle
                # if differences of tx are too small, but empirically,
                # I have not find it necessary so far)
                ordering = np.argsort(tx)                
                ijk=ijk[ordering]
                tx=tx[ordering]
                ty=ty[ordering]
                
                #r=np.sqrt(tl2[ii[2]]/tl2[ii[0]]) # ratio of longest to shortest side
                length2=np.array([(tx[1]-tx[0])**2+(ty[1]-ty[0])**2,(tx[2]-tx[1])**2+(ty[2]-ty[1])**2,(tx[0]-tx[2])**2+(ty[0]-ty[2])**2])
                r=np.sqrt(np.max(length2)/np.min(length2)) # ratio of longest to shortest side
                
                c=((tx[1]-tx[0])*(tx[2]-tx[0])+(ty[1]-ty[0])*(ty[2]-ty[0]))/np.sqrt( ((tx[1]-tx[0])**2+(ty[1]-ty[0])**2)*((tx[2]-tx[0])**2+(ty[2]-ty[0])**2)) # cos of angle of first vertex

                # orientation does not help here because many isosceles triangles, so I don't compute that
                #s=((tx[1]-tx[0])*(ty[2]-ty[0])-(tx[2]-tx[0])*(ty[1]-ty[0]))/np.sqrt( ((tx[1]-tx[0])**2+(ty[1]-ty[0])**2)*((tx[2]-tx[0])**2+(ty[2]-ty[0])**2)) # orientation whether vertices are traversed in a clockwise or counterclock-wise sense

                                
                tk.append(ijk)
                tr.append(r)
                tc.append(c)
                
    return np.array(tk),np.array(tr),np.array(tc)

def match_same_system(x1,y1,x2,y2) :
    """
    match two catalogs, assuming the coordinates are in the same coordinate system (no transfo)
    Args:
        x1 : float numpy array of coordinates along first axis of cartesian coordinate system
        y1 : float numpy array of coordinates along second axis in same system
        x2 : float numpy array of coordinates along first axis in same system
        y2 : float numpy array of coordinates along second axis in same system
    
    returns:
        indices_2 : integer numpy array. if ii is a index array for entries in the first catalog, 
                            indices_2[ii] is the index array of best matching entries in the second catalog.
                            (one should compare x1[ii] with x2[indices_2[ii]])
        distances : distances between pairs. It can be used to discard bad matches. 

    """
    xy1=np.array([x1,y1]).T
    xy2=np.array([x2,y2]).T
    tree2 = KDTree(xy2)
    distances,indices_2 = tree2.query(xy1,k=1)
    return indices_2,distances

def match_arbitrary_translation_dilatation(x1,y1,x2,y2) :
    """
    Match two catalogs in different coordinate systems, 1 and 2, related by a translation, a dilatation, and possibly a "small" rotation
    The orientation of triangles is used for the match so the rotation has to be small.
    Inspired from http://articles.adsabs.harvard.edu/pdf/1986AJ.....91.1244G
    
    Args:
        x1 : float numpy array of coordinates along first axis of cartesian coordinate system 1
        y1 : float numpy array of coordinates along second axis of cartesian coordinate system 1
        x2 : float numpy array of coordinates along first axis of cartesian coordinate system 2
        y2 : float numpy array of coordinates along second axis of cartesian coordinate system 2
    
    returns:
        indices_2 : integer numpy array. if ii is a index array for entries in the first catalog, 
                            indices_2[ii] is the index array of best matching entries in the second catalog.
                            (one should compare x1[ii] with x2[indices_2[ii]])
                            negative values for unmatched entries.
        distance : distance between pairs of triangles. It can be used to discard bad matches. 

    """

    log = get_logger()
    
    # compute all possible triangles in both data sets
    # 'tk' is index in x,y array of triangle vertices
    # 'tr' is a side length ratio
    # 'tc' a vertex cosine
    tk1,tr1,tc1 = compute_triangles_with_fixed_orientation(x1,y1)
    tk2,tr2,tc2 = compute_triangles_with_fixed_orientation(x2,y2)
    
    # we also need to use the orientation of the triangles ...
    tdu1 = np.array([x1[tk1[:,1]]-x1[tk1[:,0]],y1[tk1[:,1]]-y1[tk1[:,0]]])
    tdu1 /= np.sqrt(np.sum(tdu1**2,axis=0))
    tdu2 = np.array([x2[tk2[:,1]]-x2[tk2[:,0]],y2[tk2[:,1]]-y2[tk2[:,0]]])
    tdu2 /= np.sqrt(np.sum(tdu2**2,axis=0))
    cos12 = tdu1.T.dot(tdu2) # cosine between first side of both triangles

    # distance defined as difference of ratio of side length, cosine of primary vertex, and cosine between triangles
    # the following is the distances between all pairs of triangles
    dist2 = (tr1[:,None] - tr2.T)**2 + (tc1[:,None] - tc2.T)**2 + (np.abs(cos12)-1)**2
    matched = np.argmin(dist2,axis=1) # this is the best match
    dist = np.sqrt(np.min(dist2,axis=1)) # keep distance values

    # now that we have match of triangles , need to match back catalog entries
    ranked_pairs = np.argsort(dist)

    
    indices_2 = -1*np.ones(x1.size,dtype=int)
    distances = np.zeros(x1.size)
    
    all_matched = False
        
    for p in ranked_pairs :

        k1=tk1[p] # incides (in x1,y1) of vertices of this triangle (size=3)
        k2=tk2[matched[p]] # incides (in x2,y2) of vertices of other triangle

        #if dist2[p] > 1.e-2 : break # NEED A BETTER CRITERION THAN THIS

        # check unmatched or equal
        if np.any((indices_2[k1]>=0)&(indices_2[k1]!=k2)) :
            log.warning("skip {} <=> {}".format(k1,k2))
            continue
        indices_2[k1]=k2
        distances[k1]=dist[p]
        all_matched = (np.sum(indices_2>=0)==x1.size)
        if all_matched :
            log.debug("all matched")
            break

    # check duplicates
    for i2 in np.unique(indices_2[indices_2>=0]) :
        ii=(indices_2==i2)
        if np.sum(ii) > 1 :
            log.warning("{} duplicates for i2={}".format(np.sum(ii),i2))
            indices_2[ii]=-1
    
    return indices_2 , distances
