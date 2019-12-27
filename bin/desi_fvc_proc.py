#!/usr/bin/env python


import argparse
import sys,os
import time
import numpy as np
import fitsio

from pkg_resources import resource_filename

from astropy.table import Table
from desiutil.log import get_logger
from desimeter.detectspots import detectspots
from desimeter.findfiducials import findfiducials
#from desimeter.transform.fvc2fp.poly2d import FVCFP_Polynomial
from desimeter.transform.fvc2fp.zb import FVCFP_ZhaoBurge

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""FVC image processing""")
parser.add_argument('-i','--infile', type = str, default = None, required = True,
                    help = 'path to FVC image fits file or CSV file with spots positions')
parser.add_argument('-o','--outfile', type = str, default = None, required = True,
                    help = 'path to output CSV ASCII file')
parser.add_argument('--extname', type = str, default = 'F0000', required = False,
                    help = 'input EXTNAME to use if more than one HDU')
parser.add_argument('--output-transform', type = str, default = None, required = False,
                    help = 'write transformation to this json file')
parser.add_argument('--input-transform', type = str, default = None, required = False,
                    help = 'use this json file as input for the match, defaut is data/default-fvc2fp.json')
parser.add_argument('--threshold', type = float, default = 500., required = False,
                    help = "threshold for spots detection")

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

    spots = detectspots(image,threshold=args.threshold,nsig=7)
elif filename.find(".csv")>0 :
    log.info("read CSV spots table")
    spots = Table.read(filename,format="csv")
else :
    log.info("sorry, I don't know what to do with input file {} because not .fits nor .csv".format(filename))
    sys.exit(12)


spots = findfiducials(spots,input_transform=args.input_transform)

#tx = FVCFP_Polynomial()
tx = FVCFP_ZhaoBurge()

tx.fit(spots, update_spots=True)
# spots = fit_fvc2fp(spots)

if args.output_transform is not None :
    if not args.output_transform.endswith(".json") :
        print("error, can only write json files, so please choose an output filename end ing with .json")
    else :
        tx.write_jsonfile(args.output_transform)
        print("wrote transform in {}".format(args.output_transform))
spots.write(args.outfile,format="csv",overwrite=True)
print("wrote {}".format(args.outfile))


