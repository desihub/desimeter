#!/usr/bin/env python


import argparse
import sys,os
import numpy as np
import yaml
from desiutil.log import get_logger
from pkg_resources import resource_filename
from astropy.table import Table
from desimeter.pl2fp import apply_pl2fp

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Combine petal metrology and petal alignment data into a single CSV file""")

parser.add_argument('-o','--outfile', type = str, default = None, required = True, help = 'output csv file')


args  = parser.parse_args()
log = get_logger()

log.info("reading petal alignment data")
filename = resource_filename('desimeter',"data/petal-alignments.yaml")
log.info(" in {}".format(filename))

ifile=open(filename)
petal_alignment_dict = yaml.safe_load(ifile)
ifile.close()

log.info("reading petal metrology")
filename = resource_filename('desimeter',"data/UMT-DESI-5421-v1.csv")
log.info(" in {}".format(filename))
spots = Table.read(filename,format="csv")
r=np.sqrt(spots["X FCL"]**2+spots["Y FCL"]**2)
ii=np.where( r < 1.)[0]
for i in ii :
    log.warning("BAD DATA {} {} type={} xflc={} yflc={} xmnt={} ymnt={}".format(spots["Petal Loc ID"][i],spots["Device Loc ID"][i],spots["Device Type"][i],spots["X FCL"][i],spots["Y FCL"][i],spots["X MNT"][i],spots["Y MNT"][i]))
spots=spots[:][r>1.] # exclude bad data??

spots = apply_pl2fp(spots,petal_alignment_dict)

spots["LOCATION"] = spots["Petal Loc ID"]*1000+spots["Device Loc ID"]
if 'Dot ID' in spots.dtype.names :
    spots.rename_column('Dot ID', 'DOTID')
    log.warning("rename_column('Dot ID', 'DOTID')")

log.info("applying patch")
filename = resource_filename('desimeter',"data/fp-metrology-patch.csv")
log.info(" from {}".format(filename))
patch = Table.read(filename,format="csv")
for i in range(patch["XFP"].size) :
    selection=(spots["Petal Loc ID"]==patch["Petal Loc ID"][i])
    for k in ["Device Loc ID","DOTID"] :
        selection &= (spots[k]==patch[k][i])
    jj=np.where(selection)[0]
    if jj.size == 0 :
        log.info("Adding LOCATION={} DOTID={}".format(patch["LOCATION"][i],patch["DOTID"][i]))
        spots.add_row(patch[:][i])
    else :
        j=jj
        log.info("Replacing LOCATION={} DOTID={}".format(patch["LOCATION"][i],patch["DOTID"][i]))
        spots[:][j] = patch[:][i]

        
spots.write(args.outfile,format='csv',overwrite=True)

log.info("wrote {}".format(args.outfile))



