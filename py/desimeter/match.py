"""
Utility functions to match in a list of spots given a list of expected fiber positions 
in the FP frame (after match of fiducials and fit of transform)
"""


import numpy as np
from scipy.spatial import cKDTree as KDTree

def match(x1,y1,x2,y2) :
    """
    match two catalogs, assuming in the same coordinate system (no transfo)
    """
    xy1=np.array([x1,y1]).T
    xy2=np.array([x2,y2]).T
    tree1 = KDTree(xy1)
    distances,indices1 = tree1.query(xy2,k=1)
    return indices1,distances

