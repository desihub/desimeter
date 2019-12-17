#!/usr/bin/env python


import argparse
import sys,os
import time
import numpy as np
import fitsio

from astropy.table import Table
from desiutil.log import get_logger
from desicoord.detectspots import detectspots
from desicoord.findfiducials import findfiducials

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""FVC image processing""")
parser.add_argument('-i','--infile', type = str, default = None, required = True,
                    help = 'path to FVC image fits file or CSV file with spots positions')
parser.add_argument('-o','--outfile', type = str, default = None, required = True,
                    help = 'path to output CSV ASCII file')

args  = parser.parse_args()
log   = get_logger()

filename = args.infile
if filename.find(".fits")>0 :
    log.info("read FITS FVC image")
    image = fitsio.read(args.infile).astype(float)
    spots = detectspots(image)
elif filename.find(".csv")>0 :
    log.info("read CSV spots table")
    spots = Table.read(filename,format="csv")
else :
    log.info("sorry, I don't know what to do with input file {} because not .fits nor .csv")
    sys.exit(12)


spots = findfiducials(spots)



spots.write(args.outfile,format="csv",overwrite=True)
print("wrote {}".format(args.outfile))


