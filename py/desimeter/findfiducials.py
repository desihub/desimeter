"""
Utility functions to find fiducials in a list of spots given a know pattern of pinholes
"""

import os,sys
import numpy as np
from desimeter.log import get_logger
from astropy.table import Table,Column
from pkg_resources import resource_filename
from scipy.spatial import cKDTree as KDTree
from desimeter.transform.fvc2fp.poly2d import FVCFP_Polynomial
from desimeter.transform.fvc2fp.zb import FVCFP_ZhaoBurge


metrology_pinholes_table = None
metrology_fiducials_table = None

def compute_triangles(x,y) :
    
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


def findfiducials(spots,input_transform=None,separation=7.) :
    
    
    global metrology_pinholes_table
    global metrology_fiducials_table
    log = get_logger()
    
    log.debug("load input tranformation we will use to go from FP to FVC pixels")
    if input_transform is None :
        input_transform = resource_filename('desimeter',"data/default-fvc2fp.json")
    log.info("loading input tranform from {}".format(input_transform))
    try :
        input_tx = FVCFP_ZhaoBurge.read_jsonfile(input_transform)
    except AssertionError as e :
        log.warning("Failed to read input transfo as Zhao Burge, try polynomial...")
        input_tx = FVCFP_Polynomial.read_jsonfile(input_transform)
    
    
    if metrology_pinholes_table is None :
        
        filename = resource_filename('desimeter',"data/fp-metrology.csv")
        if not os.path.isfile(filename) :
            log.error("cannot find {}".format(filename))
            raise IOError("cannot find {}".format(filename))
        log.info("reading metrology in {}".format(filename)) 
        metrology_table = Table.read(filename,format="csv")

        log.debug("keep only the pinholes")
        metrology_pinholes_table = metrology_table[:][metrology_table["PINHOLE_ID"]>0]
        
        # use input transform to convert X_FP,Y_FP to XPIX,YPIX
        xpix,ypix = input_tx.fp2fvc(metrology_pinholes_table["X_FP"],metrology_pinholes_table["Y_FP"])
        metrology_pinholes_table["XPIX"]=xpix
        metrology_pinholes_table["YPIX"]=ypix

        log.debug("define fiducial location as central dot")
        metrology_fiducials_table = metrology_pinholes_table[:][metrology_pinholes_table["PINHOLE_ID"]==4]
    
    # find fiducials candidates  
    # select spots with at least two close neighbors (in pixel units)
    xy   = np.array([spots["XPIX"],spots["YPIX"]]).T
    tree = KDTree(xy)
    measured_spots_distances,measured_spots_indices = tree.query(xy,k=4,distance_upper_bound=separation)
    number_of_neighbors = np.sum( measured_spots_distances<separation,axis=1)
    fiducials_candidates_indices = np.where(number_of_neighbors>=3)[0]  # including self, so at least 3 pinholes
    
    # match candidates to fiducials from metrology 
    # using nearest neighbor
    
    spots_tree  = KDTree(np.array([spots["XPIX"][fiducials_candidates_indices],spots["YPIX"][fiducials_candidates_indices]]).T)
    metrology_xy   = np.array([metrology_fiducials_table["XPIX"],metrology_fiducials_table["YPIX"]]).T

    distances,indices = spots_tree.query(metrology_xy,k=1)
    log.debug("med. distance = {:4.2f} pixels for {} candidates and {} known fiducials".format(np.median(distances),fiducials_candidates_indices.size,metrology_fiducials_table["XPIX"].size))

    for loop in range(2) :
        # fit offset
        dx = np.median(metrology_xy[:,0]-spots["XPIX"][fiducials_candidates_indices][indices])
        dy = np.median(metrology_xy[:,1]-spots["YPIX"][fiducials_candidates_indices][indices])
        log.debug("offset dx={:3.1f} dy={:3.1f}".format(dx,dy))
        
        # rematch
        metrology_xy[:,0] -= dx
        metrology_xy[:,1] -= dy
        distances,indices = spots_tree.query(metrology_xy,k=1)
        log.debug("med. distance = {:4.2f} pixels for {} candidates and {} known fiducials".format(np.median(distances),fiducials_candidates_indices.size,metrology_fiducials_table["XPIX"].size))
    
    maxdistance = 10.
    selection = np.where(distances<maxdistance)[0]
    
    fiducials_candidates_indices     = fiducials_candidates_indices[indices][selection]
    matching_known_fiducials_indices = selection
    
    log.debug("mean distance = {:4.2f} pixels for {} matched and {} known fiducials".format(np.mean(distances[distances<maxdistance]),fiducials_candidates_indices.size,metrology_fiducials_table["XPIX"].size))
        
    #import matplotlib.pyplot as plt
    #plt.hist(distances[distances<maxdistance],bins=100)
    #plt.figure()
    #plt.plot(spots["XPIX"],spots["YPIX"],".")
    #plt.plot(spots["XPIX"][fiducials_candidates],spots["YPIX"][fiducials_candidates],"o")
    #plt.plot(metrology_table["XPIX"]-dx,metrology_table["YPIX"]-dy,"X",color="red")
    #plt.show()

    
    log.debug("now matching pinholes ...")
    
    nspots=spots["XPIX"].size
    if 'LOCATION' not in spots.dtype.names :
        spots.add_column(Column(np.zeros(nspots,dtype=int)),name='LOCATION')
    if 'PINHOLE_ID' not in spots.dtype.names :
        spots.add_column(Column(np.zeros(nspots,dtype=int)),name='PINHOLE_ID')
    
    
    for index1,index2 in zip ( fiducials_candidates_indices , matching_known_fiducials_indices ) :
        location = metrology_fiducials_table["LOCATION"][index2]
        
        # get indices of all pinholes for this matched fiducial
        # note we now use the full pinholes metrology table
        pi1 = measured_spots_indices[index1][measured_spots_distances[index1]<separation]
        pi2 = np.where(metrology_pinholes_table["LOCATION"]==location)[0]
        
        x1 = spots["XPIX"][pi1]
        y1 = spots["YPIX"][pi1]

        x2 = metrology_pinholes_table["XPIX"][pi2]
        y2 = metrology_pinholes_table["YPIX"][pi2]

        # inspired from http://articles.adsabs.harvard.edu/pdf/1986AJ.....91.1244G
        # compute all possible triangles in both data sets
        # 'tk' is index in x,y array of triangle vertices
        # 'tr' is a side length ratio
        # 'tc' a vertex cosine
        tk1,tr1,tc1 = compute_triangles(x1,y1)
        tk2,tr2,tc2 = compute_triangles(x2,y2)

        # we also need to use the orientation of the triangles ...
        tdu1 = np.array([x1[tk1[:,1]]-x1[tk1[:,0]],y1[tk1[:,1]]-y1[tk1[:,0]]])
        tdu1 /= np.sqrt(np.sum(tdu1**2,axis=0))
        tdu2 = np.array([x2[tk2[:,1]]-x2[tk2[:,0]],y2[tk2[:,1]]-y2[tk2[:,0]]])
        tdu2 /= np.sqrt(np.sum(tdu2**2,axis=0))
        cos12 = tdu1.T.dot(tdu2) # cosine between first side of both triangles
        
        # distance defined as difference of ratio and cosine (and cosine between triangles)
        # the following is the distances between all pairs of triangles
        dist2 = (tr1[:,None] - tr2.T)**2 + (tc1[:,None] - tc2.T)**2 + (np.abs(cos12)-1)**2
        
        matched = np.argmin(dist2,axis=1) # this is the best match
        dist2 = np.min(dist2,axis=1) # keep distance values
        
        ranked_pairs = np.argsort(dist2)
        
        metrology_pinhole_ids = metrology_pinholes_table["PINHOLE_ID"][pi2]
        
        pinhole_ids = np.zeros(x1.size,dtype=int)
        all_matched = False
        
        for p in ranked_pairs :

            k1=tk1[p] # incides (in x1,y1) of vertices of this triangle (size=3)
            k2=tk2[matched[p]] # incides (in x2,y2) of vertices of other triangle
            #print(p,k1,"<->",k2,tc1[p],dist2[p])

            #if tc1[p]>0.9 : continue # don't use ambiguous flat triangle
            if dist2[p] > 1.e-2 : break # bad pairs now
            
            # check unmatched or equal
            if np.any((pinhole_ids[k1]>0)&(pinhole_ids[k1]!=metrology_pinhole_ids[k2])) :
                log.warning("skip {} <=> {}".format(k1,k2))
                continue
            
            pinhole_ids[k1] = metrology_pinhole_ids[k2]
            spots["LOCATION"][pi1[k1]]   = location
            spots["PINHOLE_ID"][pi1[k1]] = pinhole_ids[k1]

            all_matched = (np.sum(pinhole_ids>0)==x1.size)
            if all_matched :
                log.debug("all pinholes matched for loc={}".format(location))
                allgood = True
                break
        if not all_matched :
            log.warning("only matched pinholes {} for {} detected at LOCATION {:04d} xpix~{} ypix~{}".format(pinhole_ids[pinhole_ids>0],x1.size,location,int(np.mean(x1)),int(np.mean(y1))))
            """
            import matplotlib.pyplot as plt
            plt.plot(x1,y1,"o",color="b",label="x1")
            for i in range(x1.size) :
                plt.text(x1[i]+0.3,y1[i]-0.15,"{}".format(i),color="b")
            dx=-np.mean(x2)+np.mean(x1)
            dy=-np.mean(y2)+np.mean(y1)
            plt.plot(x2+dx,y2+dy,"o",label="x2+offset",color="r")
            for i in range(x2.size) :
                plt.text(x2[i]+dx+0.3,y2[i]+dy+0.15,"{}".format(i),color="r")
            plt.plot(x1[pinhole_ids>0],y1[pinhole_ids>0],"x",c="k",label="matched")
            plt.legend()
            plt.show()
            """
        # check duplicates
        if np.unique(pinhole_ids[pinhole_ids>0]).size != np.sum(pinhole_ids>0) :
            xfp=np.mean(metrology_pinholes_table[pi2]["X_FP"])
            yfp=np.mean(metrology_pinholes_table[pi2]["Y_FP"])
            log.warning("duplicate(s) pinhole ids in {} at LOCATION={:04d} xpix~{} ypix~{} xfp~{} yfp~{}".format(pinhole_ids,location,int(np.mean(x1)),int(np.mean(y1)),int(xfp),int(yfp)))
            bc=np.bincount(pinhole_ids[pinhole_ids>0])
            duplicates = np.where(bc>1)[0]
            for duplicate in duplicates :
                log.warning("Unmatch ambiguous pinhole id = {}".format(duplicate))
                selection=(spots["LOCATION"]==location)&(spots["PINHOLE_ID"]==duplicate)
                spots["PINHOLE_ID"][selection]=0
                
            """
            import matplotlib.pyplot as plt
            plt.plot(x1,y1,"o",color="b",label="x1")
            for i in range(x1.size) :
                plt.text(x1[i]+0.3,y1[i]-0.15,"{}".format(i),color="b")
            dx=-np.mean(x2)+np.mean(x1)
            dy=-np.mean(y2)+np.mean(y1)
            plt.plot(x2+dx,y2+dy,"o",label="x2+offset",color="r")
            for i in range(x2.size) :
                plt.text(x2[i]+dx+0.3,y2[i]+dy+0.15,"{}".format(i),color="r")
            plt.plot(x1[pinhole_ids>0],y1[pinhole_ids>0],"x",c="k",label="matched")
            plt.legend()
            plt.show()
            sys.exit(12)
            """
    
    n_matched_pinholes  = np.sum(spots["PINHOLE_ID"]>0)
    n_matched_fiducials = np.sum(spots["PINHOLE_ID"]==4)
    log.info("matched {} pinholes from {} fiducials".format(n_matched_pinholes,n_matched_fiducials))
    return spots
