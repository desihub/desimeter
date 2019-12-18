#!/usr/bin/env python


import argparse
import sys,os
import numpy as np
import yaml
from desiutil.log import get_logger
from pkg_resources import resource_filename
from astropy.table import Table
from desicoord.pl2fp import apply_pl2fp

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Combine petal metrology and petal alignment data into a single CSV file""")

parser.add_argument('-o','--outfile', type = str, default = None, required = True, help = 'output csv file')


args  = parser.parse_args()
log = get_logger()

log.info("reading petal alignment data")
filename = resource_filename('desicoord',"data/petal-alignments.yaml")
log.info(" in {}".format(filename))

ifile=open(filename)
petal_alignment_dict = yaml.safe_load(ifile)
ifile.close()

log.info("reading petal metrology")
filename = resource_filename('desicoord',"data/UMT-DESI-5421-v1.csv")
log.info(" in {}".format(filename))
spots = Table.read(filename,format="csv")

spots = apply_pl2fp(spots,petal_alignment_dict)

spots["LOCATION"] = spots["Petal Loc ID"]*1000+spots["Device Loc ID"]
if 'Dot ID' in spots.dtype.names :
    spots.rename_column('Dot ID', 'DOTID')
    log.warning("rename_column('Dot ID', 'DOTID')")

spots.write(args.outfile,format='csv',overwrite=True)

log.info("wrote {}".format(args.outfile))



