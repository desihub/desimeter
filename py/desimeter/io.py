from pkg_resources import resource_filename
from astropy.table import Table

def load_metrology():
    '''
    Returns metrology table
    '''
    filename = resource_filename('desimeter', 'data/fp-metrology.csv')
    metrology = Table.read(filename)
    return metrology
