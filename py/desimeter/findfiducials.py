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
from desimeter.match import match_same_system,match_arbitrary_translation_dilatation
from desimeter.simplecorr import SimpleCorr

metrology_pinholes_table = None
metrology_fiducials_table = None

def findfiducials(spots,input_transform=None,separation=8.) :
    
    
    global metrology_pinholes_table
    global metrology_fiducials_table
    log = get_logger()
    
    log.debug("load input tranformation we will use to go from FP to FVC pixels")
    if input_transform is None :
        input_transform = resource_filename('desimeter',"data/single-lens-fvc2fp.json")
    
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
    log.info("select spots with at least two close neighbors (in pixel units)")
    xy   = np.array([spots["XPIX"],spots["YPIX"]]).T
    tree = KDTree(xy)
    measured_spots_distances,measured_spots_indices = tree.query(xy,k=4,distance_upper_bound=separation)
    number_of_neighbors = np.sum( measured_spots_distances<separation,axis=1)
    fiducials_candidates_indices = np.where(number_of_neighbors>=3)[0]  # including self, so at least 3 pinholes
    
    # match candidates to fiducials from metrology

    log.info("first match {} fiducials candidates to metrology ({}) with iterative fit".format(fiducials_candidates_indices.size,len(metrology_fiducials_table)))
    x1 = spots["XPIX"][fiducials_candidates_indices]
    y1 = spots["YPIX"][fiducials_candidates_indices]
    x2 = metrology_fiducials_table["XPIX"] # do I need to do this?
    y2 = metrology_fiducials_table["YPIX"]
    
    nloop=20
    saved_median_distance=0
    for loop in range(nloop) :
        indices_2, distances = match_same_system(x1,y1,x2,y2)
        mdist = np.median(distances[indices_2>=0])
        if loop < nloop-1 :
            maxdistance = max(10,3.*1.4*mdist)
        else : # final iteration
            maxdistance = 10 # pixel
        selection = np.where((indices_2>=0)&(distances<maxdistance))[0]
        log.info("iter #{} median_dist={} max_dist={} matches={}".format(loop,mdist,maxdistance,selection.size))
        corr21 = SimpleCorr()
        corr21.fit(x2[indices_2[selection]],y2[indices_2[selection]],x1[selection],y1[selection])
        x2,y2 = corr21.apply(x2,y2)
        if np.abs(saved_median_distance-mdist)<0.0001 : break # no more improvement
        saved_median_distance = mdist
    
    # use same coord system match (note we now match the otherway around)
    indices_1, distances = match_same_system(x2,y2,x1,y1)
    maxdistance = 10. # FVC pixels
    selection = np.where((indices_1>=0)&(distances<maxdistance))[0]
    fiducials_candidates_indices     = fiducials_candidates_indices[indices_1[selection]]
    matching_known_fiducials_indices = selection
    
    log.debug("mean distance = {:4.2f} pixels for {} matched and {} known fiducials".format(np.mean(distances[distances<maxdistance]),fiducials_candidates_indices.size,metrology_fiducials_table["XPIX"].size))
    
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
        
        indices_2 , distances = match_arbitrary_translation_dilatation(x1,y1,x2,y2)
        
        metrology_pinhole_ids = metrology_pinholes_table["PINHOLE_ID"][pi2]
        pinhole_ids = np.zeros(x1.size,dtype=int)
        matched=(indices_2>=0)
        pinhole_ids[matched] = metrology_pinhole_ids[indices_2[matched]]

        spots["LOCATION"][pi1[matched]]   = location 
        spots["PINHOLE_ID"][pi1[matched]] = pinhole_ids[matched]
        
        if np.sum(pinhole_ids==0) > 0 :
            log.warning("only matched pinholes {} for {} detected at LOCATION {} xpix~{} ypix~{}".format(pinhole_ids[pinhole_ids>0],x1.size,location,int(np.mean(x1)),int(np.mean(y1))))
        
        # check duplicates
        if np.unique(pinhole_ids[pinhole_ids>0]).size != np.sum(pinhole_ids>0) :
            xfp=np.mean(metrology_pinholes_table[pi2]["X_FP"])
            yfp=np.mean(metrology_pinholes_table[pi2]["Y_FP"])
            log.warning("duplicate(s) pinhole ids in {} at LOCATION={} xpix~{} ypix~{} xfp~{} yfp~{}".format(pinhole_ids,location,int(np.mean(x1)),int(np.mean(y1)),int(xfp),int(yfp)))
            bc=np.bincount(pinhole_ids[pinhole_ids>0])
            duplicates = np.where(bc>1)[0]
            for duplicate in duplicates :
                log.warning("Unmatch ambiguous pinhole id = {}".format(duplicate))
                selection=(spots["LOCATION"]==location)&(spots["PINHOLE_ID"]==duplicate)
                spots["PINHOLE_ID"][selection]=0
        
    spots["PETAL_LOC"]=spots["LOCATION"]//1000
    spots["DEVICE_LOC"]=spots["LOCATION"]%1000
    
    n_matched_pinholes  = np.sum(spots["PINHOLE_ID"]>0)
    n_matched_fiducials = np.sum(spots["PINHOLE_ID"]==4)
    log.info("matched {} pinholes from {} fiducials".format(n_matched_pinholes,n_matched_fiducials))
    return spots
