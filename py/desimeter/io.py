import os
from pkg_resources import resource_filename
from astropy.table import Table
import yaml
from desimeter.log import get_logger

def desimeter_data_dir():
    '''
    Returns desimeter data dir
    '''
    if "DESIMETER_DATA" in os.environ :
        return os.environ["DESIMETER_DATA"]
    else :
        return resource_filename('desimeter', 'data')

def load_metrology():
    '''
    Returns metrology table
    '''
    filename = os.path.join(desimeter_data_dir(),'fp-metrology.csv')
    log=get_logger()
    log.debug("loading {}".format(filename))
    
    metrology = Table.read(filename)
    return metrology

def load_petal_alignement():
    filename = os.path.join(desimeter_data_dir(),'petal-alignments.yaml')
    with open(filename) as ifile :
        petal_alignment_dict = yaml.safe_load(ifile)
    return petal_alignment_dict

def fvc2fp_filename():
    return os.path.join(desimeter_data_dir(),'init-fvc2fp.json')
