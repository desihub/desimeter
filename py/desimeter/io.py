import os
from pkg_resources import resource_filename
from astropy.table import Table
import yaml
from desimeter.log import get_logger

def read_hexrot_deg(header) :
    '''Extract from header the hexapod rotation angle (works with both fitsio version 0.9 and 1.1)

    Args:
        header : fits header or dictionary

    Returns:
       hexapod rotation angle in degree
    '''
    log = get_logger()

    focus_params = header["FOCUS"]
    if isinstance(focus_params,str) :
        focus_params = focus_params.split(",")
    elif not isinstance(focus_params,tuple) :
        message="header['FOCUS']={} is not a str nor a tuple".format(focus_params)
        log.error(message)
        raise ValueError(message)
    if len(focus_params) != 6 :
        message="I expected 6 coefficient in header['FOCUS']={}".format(focus_params)
        log.error(message)
        raise ValueError(message)
    hexrot_deg = float(focus_params[5])/3600.
    log.debug("Use hexrot_deg={} from FOCUS={}".format(hexrot_deg,header["FOCUS"]))
    return hexrot_deg



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
    Returns metrology table.
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

def fvc_bias_filename():
    return os.path.join(desimeter_data_dir(),'bias.fits')
