import os
import sys
import subprocess
from astropy.table import Table

def process_fvc(filename, overwrite=False):
    '''Calls desi_fvc_proc on the input file.

    INPUTS:  filename ... fits file (can also be csv, in which case desi_fvc_proc is skipped)
             overwrite ... if a processed table of the same basename already exists in /tmp,
                           reprocess and overwrite it with new table. (default behavior is
                           to use the already-existing table)

    OUTPUT:  astropy table
    '''
    if filename.find(".csv")>0 :
        # already a coordinates table
        return Table.read(filename)
    if filename.find(".fits")<0 :
        print("don't know what to do with",filename)
        sys.exit(12)
    outfilename = get_outfilename(filename)
    if os.path.isfile(outfilename) and not overwrite:
        print("using previously processed {}".format(outfilename))
    else :
        cmd="desi_fvc_proc -i {} -o {}".format(filename,outfilename)
        print("running '{}'".format(cmd))
        errcode=subprocess.call(cmd.split())
        if errcode != 0 :
            print("we got an error code = {}".format(errcode))
            sys.exit(errcode)
    return Table.read(outfilename)

def get_outfilename(path_to_fits):
    out_dir = '/tmp'
    csv_name= path_to_fits.replace(".fits.fz",".csv").replace(".fits",".csv")
    basename = os.path.basename(csv_name)
    out = os.path.join(out_dir, basename)
    return out
