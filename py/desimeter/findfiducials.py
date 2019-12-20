"""
Utility functions to find fiducials in a list of spots given a know pattern of pinholes
"""

import os
import numpy as np
from desiutil.log import get_logger
from astropy.table import Table,Column
from pkg_resources import resource_filename
from scipy.spatial import cKDTree as KDTree
from desimeter.fvc2fp import FVCFP_Polynomial


metrology_table = None

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
    if 'DOTID' not in spots.dtype.names :
        spots.add_column(Column(np.zeros(nspots,dtype=int)),name='DOTID')
    
    
    for index1,index2 in zip ( fiducials_candidates_indices , matching_known_fiducials_indices ) :
        metrology_index1 = measured_spots_indices[index1][measured_spots_distances[index1]<separation]
        metrology_index2 = np.where(metrology_table["LOCATION"]==metrology_table["LOCATION"][index2])[0]

        dx=spots["XPIX"][index1]-metrology_table["XPIX"][index2]
        dy=spots["YPIX"][index1]-metrology_table["YPIX"][index2]

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
        spots["DOTID"][metrology_index1] = metrology_table["DOTID"][metrology_index2][matched_indices]
        
    
    return spots
