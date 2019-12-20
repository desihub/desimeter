#!/usr/bin/env python


import argparse
import sys,os
import time
import numpy as np
import fitsio

from astropy.table import Table
from desiutil.log import get_logger
from desimeter.detectspots import detectspots
from desimeter.findfiducials import findfiducials
from desimeter.fvc2fp import FVCFP_Polynomial

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""FVC image processing""")
parser.add_argument('-i','--infile', type = str, default = None, required = True,
                    help = 'path to FVC image fits file or CSV file with spots positions')
parser.add_argument('-o','--outfile', type = str, default = None, required = True,
                    help = 'path to output CSV ASCII file')
parser.add_argument('--extname', type = str, default = 'F0000', required = False,
                    help = 'input EXTNAME to use if more than one HDU')

args  = parser.parse_args()
log   = get_logger()

filename = args.infile
if filename.find(".fits")>0 :
    log.info("read FITS FVC image")
    with fitsio.FITS(args.infile) as fx:
        if len(fx) == 1:
            image = fx[0].read().astype(float)
        else:
            image = fx[args.extname].read().astype(float)

    spots = detectspots(image)
elif filename.find(".csv")>0 :
    log.info("read CSV spots table")
    spots = Table.read(filename,format="csv")
else :
    log.info("sorry, I don't know what to do with input file {} because not .fits nor .csv".format(filename))
    sys.exit(12)


spots = findfiducials(spots)

tx = FVCFP_Polynomial()
tx.fit(spots, update_spots=True)
# spots = fit_fvc2fp(spots)

spots.write(args.outfile,format="csv",overwrite=True)
print("wrote {}".format(args.outfile))


