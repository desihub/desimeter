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
    
    nn=len(x)
    ntri=(nn*(nn-1)*(nn-2))//6
    tk=np.zeros((ntri,3),dtype=int) # indices 
    txyz=np.zeros((ntri,3),dtype=float) # properties of the shape and orientation of triangles
        
    # I can do this ordering once, outside of the loop on triangles, to got faster
    ordering = np.argsort(x) 
    x=x[ordering]
    y=y[ordering]
    tx=np.zeros(3)
    ty=np.zeros(3)

    triangle_index=0
    for i in range(nn) :
        tx[0]=x[i]
        ty[0]=y[i]
        for j in range(i+1,nn) :
            tx[1]=x[j]
            ty[1]=y[j]
            for k in range(j+1,nn) :
                # x y of vertices
                tx[2]=x[k]
                ty[2]=y[k]
                length2=np.array([(tx[1]-tx[0])**2+(ty[1]-ty[0])**2,(tx[2]-tx[1])**2+(ty[2]-ty[1])**2,(tx[0]-tx[2])**2+(ty[0]-ty[2])**2])
                r=np.sqrt(np.max(length2)/np.min(length2)) # ratio of longest to shortest side
                tk[triangle_index]=ordering[[i,j,k]]
                txyz[triangle_index,0]=np.sqrt(np.max(length2)/np.min(length2)) # ratio of longest to shortest side
                txyz[triangle_index,1]=((tx[1]-tx[0])*(tx[2]-tx[0])+(ty[1]-ty[0])*(ty[2]-ty[0]))/np.sqrt( ((tx[1]-tx[0])**2+(ty[1]-ty[0])**2)*((tx[2]-tx[0])**2+(ty[2]-ty[0])**2)) # cos of angle of first vertex
                txyz[triangle_index,2]=(tx[1]-tx[0])/np.sqrt((tx[1]-tx[0])**2+(ty[1]-ty[0])**2) # cos of angle of first side
                triangle_index += 1
    return tk,txyz

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
    # txyz are properties of the shape and orientation of the triangles
    log.debug("compute triangles")
    tk1,txyz1 = compute_triangles_with_fixed_orientation(x1,y1)
    tk2,txyz2 = compute_triangles_with_fixed_orientation(x2,y2)
    
    log.debug("match triangles")
    # match with kdtree triangles with same shape and orientation
    tree2=KDTree(txyz2)
    triangle_distances,triangle_indices_2 = tree2.query(txyz1,k=1)
    
    # now that we have match of triangles , need to match back catalog entries
    ranked_pairs = np.argsort(triangle_distances)
    
    indices_2 = -1*np.ones(x1.size,dtype=int)
    distances = np.zeros(x1.size)
    
    all_matched = False
    log.debug("match catalogs using pairs of triangles")
    for p in ranked_pairs :

        k1=tk1[p] # incides (in x1,y1) of vertices of this triangle (size=3)
        k2=tk2[triangle_indices_2[p]] # incides (in x2,y2) of vertices of other triangle
        
        # check unmatched or equal
        if np.any((indices_2[k1]>=0)&(indices_2[k1]!=k2)) :
            log.warning("skip {} <=> {}".format(k1,k2))
            continue
        indices_2[k1]=k2
        distances[k1]=triangle_distances[p]
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
