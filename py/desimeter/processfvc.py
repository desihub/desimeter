import os
import subprocess
from astropy.table import Table

def process_fvc(filename) :
    if filename.find(".csv")>0 :
        # already a coordinates table
        return Table.read(filename)
    if filename.find(".fits")<0 :
        print("don't know what to do with",filename)
        sys.exit(12)
    outfilename=os.path.join("/tmp",os.path.basename(filename.replace(".fits",".csv")))
    if os.path.isfile(outfilename) :
        print("using previously processed {}".format(outfilename))
    else :
        cmd="desi_fvc_proc -i {} -o {}".format(filename,outfilename)
        print("running '{}'".format(cmd))
        errcode=subprocess.call(cmd.split())
        if errcode != 0 :
            print("we got an error code = {}".format(errcode))
            sys.exit(errcode)
    return Table.read(outfilename)
