"""
Utility functions to find fiducials in a list of spots given a know pattern of pinholes
"""

import os,sys
import numpy as np
from desiutil.log import get_logger
from astropy.table import Table,Column
from pkg_resources import resource_filename
from scipy.spatial import cKDTree as KDTree
from desimeter.transform.fvc2fp.poly2d import FVCFP_Polynomial


metrology_table = None


def compute_triangles(x,y) :
    
    tti=[] # indices
    ttx=[]
    tty=[]
    tr=[]
    tc=[]
    
    nn=len(x)
    for i in range(nn) :
        for j in range(i+1,nn) :
            for k in range(j+1,nn) :
                # x y of vertices
                ti=[i,j,k]
                tx=x[[i,j,k]]
                ty=y[[i,j,k]]

                # sort according to length
                # length of sides square
                tl2=np.array([(tx[1]-tx[0])**2+(ty[1]-ty[0])**2,(tx[2]-tx[1])**2+(ty[2]-ty[1])**2,(tx[0]-tx[2])**2+(ty[0]-ty[2])**2])
                pairs=np.array([[0,1],[1,2],[0,2]])
                
                ii=np.argsort(tl2)
                ordering = np.zeros(3).astype(int)
                ordering[0] = np.intersect1d(pairs[ii[0]],pairs[ii[2]]) # vertex connected to shortest and longest side  
                ordering[1] = np.intersect1d(pairs[ii[0]],pairs[ii[1]]) # vertex connected to shortest and intermediate side  
                ordering[2] = np.intersect1d(pairs[ii[1]],pairs[ii[2]]) # vertex connected to intermediate and longest side  
                tx=tx[ordering]
                ty=ty[ordering]
                #print(np.sqrt(np.array([(tx[1]-tx[0])**2+(ty[1]-ty[0])**2,(tx[2]-tx[1])**2+(ty[2]-ty[1])**2,(tx[0]-tx[2])**2+(ty[0]-ty[2])**2])))
                r=np.sqrt(tl2[ii[2]]/tl2[ii[0]]) # ratio of longest to shortest side
                c=((tx[1]-tx[0])*(tx[2]-tx[0])+(ty[1]-ty[0])*(ty[2]-ty[0]))/np.sqrt( ((tx[1]-tx[0])**2+(ty[1]-ty[0])**2)*((tx[2]-tx[0])**2+(ty[2]-ty[0])**2)) # cos of angle of first vertex
                # orientation does not help here because many symmetric triangles, so I don't compute that
                #s=((tx[1]-tx[0])*(ty[2]-ty[0])-(tx[2]-tx[0])*(ty[1]-ty[0]))/np.sqrt( ((tx[1]-tx[0])**2+(ty[1]-ty[0])**2)*((tx[2]-tx[0])**2+(ty[2]-ty[0])**2)) # orientation whether vertices are traversed in a clockwise or counterclock-wise sense
                #print("r=",r,"c=",c)
                                
                tti.append(ti)
                ttx.append(tx)
                tty.append(ty)
                tr.append(r)
                tc.append(c)
                
    return np.array(tti),np.array(ttx),np.array(tty),np.array(tr),np.array(tc)


def findfiducials(spots,input_transform=None,separation=7.) :
    
    
    global metrology_table
    
    log = get_logger()
    log.info("findfiducials...")

    if input_transform is None :
        input_transform = resource_filename('desimeter',"data/default-fvc2fp.json")
    log.info("loading input tranform from {}".format(input_transform))
    input_tx = FVCFP_Polynomial.read_jsonfile(input_transform)

    
    if metrology_table is None :
        
        filename = resource_filename('desimeter',"data/fp-metrology.csv")
        if not os.path.isfile(filename) :
            log.error("cannot find {}".format(filename))
            raise IOError("cannot find {}".format(filename))
        log.info("reading metrology in {}".format(filename)) 
        metrology_table = Table.read(filename,format="csv")

        # keep only the pinholes
        metrology_table = metrology_table[:][metrology_table["DOTID"]>0]
        
        # use input transform to convert XFP,YFP to XPIX,YPIX
        xpix,ypix = input_tx.fp2fvc(metrology_table["XFP"],metrology_table["YFP"])
        metrology_table["XPIX"]=xpix
        metrology_table["YPIX"]=ypix
        
        
    # find fiducials candidates  
    # select spots with at least two close neighbors (in pixel units)
    xy   = np.array([spots["XPIX"],spots["YPIX"]]).T
    tree = KDTree(xy)
    measured_spots_distances,measured_spots_indices = tree.query(xy,k=4,distance_upper_bound=separation)
    number_of_neighbors = np.sum( measured_spots_distances<separation,axis=1)
    fiducials_candidates_indices = np.where(number_of_neighbors>=3)[0]  # including self, so at least 3 pinholes
    
    # match candidates to known fiducials
    # nearest neighbor
    
    spots_tree  = KDTree(np.array([spots["XPIX"][fiducials_candidates_indices],spots["YPIX"][fiducials_candidates_indices]]).T)
    xy   = np.array([metrology_table["XPIX"],metrology_table["YPIX"]]).T

    distances,indices = spots_tree.query(xy,k=1)
    log.debug("0- median distance = {:4.2f} pixels for {} candidate fiducials and {} known fiducials".format(np.median(distances),fiducials_candidates_indices.size,metrology_table["XPIX"].size))
    
    # fit offset
    dx = np.median(metrology_table["XPIX"]-spots["XPIX"][fiducials_candidates_indices][indices])
    dy = np.median(metrology_table["YPIX"]-spots["YPIX"][fiducials_candidates_indices][indices])
    
    # rematch
    distances,indices = spots_tree.query(np.array([metrology_table["XPIX"]-dx,metrology_table["YPIX"]-dy]).T,k=1)
    log.debug("1- median distance = {:4.2f} pixels for {} candidate fiducials and {} known fiducials".format(np.median(distances),fiducials_candidates_indices.size,metrology_table["XPIX"].size))
    rms = 1.48*np.median(distances)
    # print(indices)
    
    maxdistance = np.sqrt((3*rms)**2+2.**2) # keep fiducials within 3 sigma + systematic of 2 pixels
    log.debug("max distance = {} pixel".format(maxdistance))

    selection = np.where(distances<maxdistance)[0]
    
    fiducials_candidates_indices     = fiducials_candidates_indices[indices][selection]
    matching_known_fiducials_indices = selection
    
    log.debug("2- mean distance = {:4.2f} pixels for {} matched fiducials and {} known fiducials".format(np.mean(distances[distances<maxdistance]),fiducials_candidates_indices.size,metrology_table["XPIX"].size))
    
    
    
    #import matplotlib.pyplot as plt
    #plt.hist(distances[distances<maxdistance],bins=100)
    #plt.figure()
    #plt.plot(spots["XPIX"],spots["YPIX"],".")
    #plt.plot(spots["XPIX"][fiducials_candidates],spots["YPIX"][fiducials_candidates],"o")
    #plt.plot(metrology_table["XPIX"]-dx,metrology_table["YPIX"]-dy,"X",color="red")
    #plt.show()

    
    # now, identify pinholes

    nspots=spots["XPIX"].size
    if 'LOCATION' not in spots.dtype.names :
        spots.add_column(Column(np.zeros(nspots,dtype=int)),name='LOCATION')
    if 'PINHOLE_ID' not in spots.dtype.names :
        spots.add_column(Column(np.zeros(nspots,dtype=int)),name='PINHOLE_ID')
    
    
    for index1,index2 in zip ( fiducials_candidates_indices , matching_known_fiducials_indices ) :
        location = metrology_table["LOCATION"][index2]
        metrology_index1 = measured_spots_indices[index1][measured_spots_distances[index1]<separation]
        metrology_index2 = np.where(metrology_table["LOCATION"]==location)[0]
        print("LOCATION = ",location)

        x1 = spots["XPIX"][metrology_index1]
        y1 = spots["YPIX"][metrology_index1]

        x2 = metrology_table["XPIX"][metrology_index2]
        y2 = metrology_table["YPIX"][metrology_index2]

        # using http://articles.adsabs.harvard.edu/pdf/1986AJ.....91.1244G
        # compute all possible triangles in both dat sets
        # also computes a side length ratio 'r', a vertex cosine 'c', and an orientation 's' (see routine compute_triangles)
        tti1,ttx1,tty1,tr1,tc1 = compute_triangles(x1,y1)
        tti2,ttx2,tty2,tr2,tc2 = compute_triangles(x2,y2)
        
        # distance defined as difference of ratio and cosine
        # ignore orientation because many symmetric triangles 
        # could apply weights ...
        # the following is distances of all pairs of triangles
        dist2 = (tr1[:,None] - tr2.T)**2 + (tc1[:,None] - tc2.T)**2
        
        jj = np.argmin(dist2,axis=1) # this is the match

        dist2 = np.min(dist2,axis=1)

        goodpairs = np.where(dist2<1e-5)[0]

        matched = np.zeros(x1.size)
        for goodpair in goodpairs :
            print("a good pair")
            i1=tti1[goodpair]
            i2=tti2[jj[goodpair]]
            spots["LOCATION"][metrology_index1][i1] = location
            spots["PINHOLE_ID"][metrology_index1][i1] = metrology_table["PINHOLE_ID"][metrology_index2][i2]
            matched[i1] += 1
            if np.sum(matched>0)==x1.size : # all matched
                print("all pinholes matched for fiducial",location,"; pinholes=",spots["PINHOLE_ID"][metrology_index1])
                break
            
        
        """ 
        #for pair in dist2<1e-5
        #print(dist2)
        
        # get from this the translation
        dx = np.median(ttx1-ttx2[jj])
        dy = np.median(tty1-tty2[jj])

        
        import matplotlib.pyplot as plt
        plt.plot(x1,y1,"o",c="blue")
        plt.plot(x2+dx,y2+dy,"x",c="red")
        plt.show()
        
        
        import matplotlib.pyplot as plt
        
        for i,j in enumerate(jj) :
            plt.plot(x1,y1,"o",c="blue")
            plt.plot(x2,y2,"o",c="red")

            i1=i
            i2=j
            print("R: ",tr1[i1],tr2[i2])
            print("C: ",tc1[i1],tc2[i2])
            #print("S: ",ts1[i1],ts2[i2])
            print("dist2=",dist2[i,jj])
            plt.plot(ttx1[i1],tty1[i1],"-",c="blue",alpha=0.7)
            plt.plot(ttx2[i2],tty2[i2],"-",c="red",alpha=0.7)
            plt.plot(ttx1[i1][0],tty1[i1][0],"x",c="k")
            plt.plot(ttx2[i2][0],tty2[i2][0],"x",c="k")
            for k in range(3) :
                plt.plot([ ttx1[i1,k], ttx2[i2,k] ] , [ tty1[i1,k], tty2[i2,k] ] ,"-",c="gray")
            plt.show()
        
        
        

        for i1 in range(len(tri1)) :
            i2 = np,a
        print(tri1)
        sys.exit(12)
        
        
        xy_metro = np.array([metrology_table["XPIX"][metrology_index2]+dx,metrology_table["YPIX"][metrology_index2]+dy]).T
        metrology_tree = KDTree(xy_metro)
        xy_meas  = np.array([spots["XPIX"][metrology_index1],spots["YPIX"][metrology_index1]]).T
        distances,matched_indices = metrology_tree.query(xy_meas,k=1)

        if True : # median offset and then rematch
            dx = np.median( spots["XPIX"][metrology_index1] - metrology_table["XPIX"][metrology_index2][matched_indices] )
            dy = np.median( spots["YPIX"][metrology_index1] - metrology_table["YPIX"][metrology_index2][matched_indices] )
            xy_metro = np.array([metrology_table["XPIX"][metrology_index2]+dx,metrology_table["YPIX"][metrology_index2]+dy]).T
            metrology_tree = KDTree(xy_metro)
            distances,matched_indices = metrology_tree.query(xy_meas,k=1)

        spots["LOCATION"][metrology_index1] = metrology_table["LOCATION"][index2]
        spots["PINHOLE_ID"][metrology_index1] = metrology_table["PINHOLE_ID"][metrology_index2][matched_indices]
        """
    
    return spots
