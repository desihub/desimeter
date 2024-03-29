"""FVC image processing
"""
import os
import sys
import subprocess
import argparse
from astropy.table import Table
import numpy as np
import fitsio

from desimeter.log import get_logger
from desimeter.detectspots import detectspots
from desimeter.findfiducials import findfiducials
from desimeter.transform.fvc2fp import FVC2FP
from desimeter.match import match_same_system
from desimeter.transform.xy2qs import qs2xy
from desimeter.fieldmodel import FieldModel
from desimeter.io import load_metrology,fvc2fp_filename,fvc_bias_filename
from desimeter.spotmatch import spotmatch
from desimeter.turbulence import correct,correct_with_pol

def process_fvc(filename, overwrite=False, use_subprocess=True):
    """Process an FVC image.

    The default use_subprocess=True maintains backwards compatibility but
    new code should generally use use_subprocess=False to avoid the overhead
    of a subprocess call and for more predictable results.

    Parameters
    ----------
    filename: str
        Name of FITS file containing the FVC image to analyze, or else a CSV file of
        spot coordinates to read and return directly.
    overwrite: bool
        Read an ouput file from a previous call to this function when False.
        Otherwise, always generate the output to return.
    use_subprocess: bool
        Calls desi_fvc_proc in a subprocess when True.  Otherwise, process the FVC image
        in the current context. A subprocess will generally have a different python
        environment, especially if the current environment is managed by conda or
        a jupyter kernel.

    Returns
    -------
    astropy.table.Table
        A table of the found spot coordinates.
    """
    if filename.find(".csv")>0 :
        # already a coordinates table
        return Table.read(filename)
    if filename.find(".fits")<0 :
        print("don't know what to do with",filename)
        sys.exit(12)
    outfilename = get_outfilename(filename)
    if os.path.isfile(outfilename) and not overwrite:
        print("using previously processed {}".format(outfilename))
        return Table.read(outfilename)
    if use_subprocess:
        # Use a subprocess to invoke bin/desi_fvc_proc in a login environment.
        cmd="desi_fvc_proc -i {} -o {}".format(filename,outfilename)
        print("running '{}'".format(cmd))
        errcode=subprocess.call(cmd.split())
    else:
        # Use the current environment.
        parser = get_parser()
        args = parser.parse_args(['-i', filename, '-o', outfilename])
        log   = get_logger()
        errcode = fvc_proc(args, log)

    if errcode != 0 :
        print("we got an error code = {}".format(errcode))
        return None
    return Table.read(outfilename)

def get_outfilename(path_to_fits):
    out_dir = '/tmp'
    csv_name= path_to_fits.replace(".fits.fz",".csv").replace(".fits",".csv")
    basename = os.path.basename(csv_name)
    out = os.path.join(out_dir, basename)
    return out


def preproc(args):
    """Preprocess arguments.
    """
    if not args.outfile.endswith(".csv") :
        print("sorry output filename has to end with .csv")
        return 12

    if args.threshold is not None :
        print("sorry, please use option --min-counts-per-pixel instead of ambiguous and now deprecated --threshold.")
        return 1

    if args.no_bias :
        args.bias = None
    else :
        if args.bias :
            bias_filename = args.bias
        else :
            bias_filename = fvc_bias_filename()
        if os.path.isfile(bias_filename) :
            print("Will subtract bias. Reading {}".format(bias_filename))
            args.bias = fitsio.read(bias_filename).astype(float)
        else :
            args.bias = None

    if args.make_directory:
        directory, _ = os.path.split(args.outfile)
        if directory != '':
            if not os.path.exists(directory):
                os.makedirs(directory)
        if args.output_transform is not None:
            directory, _ = os.path.split(args.outfile)
            if directory != '':
                if not os.path.exists(directory):
                    os.makedirs(directory)

    return 0


def get_spots_list(args, log):

    spots_list=[]
    filename = args.infile
    extnames = None
    if filename.find(".fits")>0 :
        log.info("read FITS FVC image")
        with fitsio.FITS(args.infile) as fx:
            print(fx) # shockingly needed to get the attribute hdu_map
            extnames = []
            if args.sequence :
                for extname in fx.hdu_map.keys() :
                    if isinstance(extname, str) and extname.lower().find("f0")==0 :
                        extnames.append(extname)
            elif len(fx) == 1:
                extnames.append(0)
            elif args.extname.strip() != 'last':
                extnames.append(args.extname)
            else :
                tmp_extnames = [k for k in fx.hdu_map.keys() if isinstance(k, str)]
                tmp_extnames = [e.lower() for e in tmp_extnames]
                extname = None
                for i in range(len(tmp_extnames)):
                    extname0 = f'f{i:04d}'
                    if extname0 in tmp_extnames:
                        extname = extname0
                extnames.append(extname)

            for extname in extnames :
                log.info('reading image in extension {}'.format(extname))
                image = fx[extname].read().astype(float)
                if args.bias is not None :
                    image -= args.bias
                spots_list.append( detectspots(image,min_counts_per_pixel=args.min_counts_per_pixel,min_counts_per_spot=args.min_counts,nsig=7,psf_sigma=args.fvc_psf_sigma) )

    elif filename.find(".csv")>0 :
        log.info("read CSV spots table")
        spots_list.append( Table.read(filename,format="csv") )
    else :
        log.info("sorry, I don't know what to do with input file {} because not .fits nor .csv".format(filename))
        return None, None

    return spots_list, extnames


def get_expected_pos(args, log):
    if args.expected_positions is not None :
        log.info("reading expected positions in {}".format(args.expected_positions))
        expected_pos = Table.read(args.expected_positions) # auto-magically guess format

        if (not "X_FP" in expected_pos.keys()) and "X_FP_EXP" in expected_pos.keys() :
            log.warning("Rename X_FP_EXP,Y_FP_EXP -> X_FP,Y_FP")
            expected_pos.rename_column("X_FP_EXP","X_FP")
            expected_pos.rename_column("Y_FP_EXP","Y_FP")


        if not "X_FP" in expected_pos.keys() :
            if "EXP_Q_0" in  expected_pos.keys() :
                log.info("EXP_Q_0,EXP_S_0 -> X_FP,Y_FP")
                x,y = qs2xy(q=expected_pos["EXP_Q_0"],s=expected_pos["EXP_S_0"])
                expected_pos["X_FP"] = x
                expected_pos["Y_FP"] = y
                bad_expected_pos = (np.isnan(x)|np.isnan(y))
                if np.sum(bad_expected_pos) > 0 :
                    expected_pos["X_FP"][bad_expected_pos]=-99999.
                    expected_pos["Y_FP"][bad_expected_pos]=-99999.
            else :
                log.error("No EXP_Q_0 nor X_FP in expected positions file {}".format(args.expected_positions))
    else :
        log.info("since no input expected positions, use metrology to match the fibers to the positioner centers")
        expected_pos = load_metrology()

    if not "LOCATION" in expected_pos.keys() :
        # add useful location keyword
        expected_pos["LOCATION"] = np.array(expected_pos["PETAL_LOC"])*1000+np.array(expected_pos["DEVICE_LOC"])

    if "PINHOLE_ID" in expected_pos.dtype.names :
        # exclude pinhole because here we want to match fibers
        ii = np.where(expected_pos["PINHOLE_ID"]==0)[0]
        expected_pos = expected_pos[:][ii]

    return expected_pos


def fvc_proc(args, log):
    """Process an FVC image with options specified by args and output to log.
    """
    errcode = preproc(args)
    if errcode:
        return errcode

    spots_list, extnames = get_spots_list(args, log)
    if spots_list is None:
        return 13

    for seqid,spots in enumerate(spots_list) :

        if args.min_spots is not None :
            if len(spots) < args.min_spots :
                log.error("not enough spots, exiting")
                continue
        if args.max_spots is not None :
            if len(spots) > args.max_spots :
                log.error("too many spots, exiting")
                continue

        outfile = args.outfile
        if args.sequence:
            if extnames is not None:
                outfile = outfile.replace(".csv","-%s.csv" % extnames[seqid])
            else:
                outfile = outfile.replace(".csv","-F{:04d}.csv".format(seqid))

        if args.nomatch :
            # write spots
            spots.write(outfile,format="csv",overwrite=True)
            print("wrote {}".format(outfile))
            continue

        if args.use_spotmatch :
            spots = spotmatch(spots["XPIX"],spots["YPIX"],match_radius_pixels=args.spotmatch_match_radius_pixels)
            # drop pinholes, keep only center
            spots = spots[(spots["PINHOLE_ID"]==0)|(spots["PINHOLE_ID"]==99)]
            n_matched_fiducials = np.sum(spots['PINHOLE_ID'] == 99)
        else :
            spots = findfiducials(spots,input_transform=args.input_transform,pinhole_max_separation_mm=args.pinhole_max_separation_mm)
            n_matched_fiducials = np.sum(spots['PINHOLE_ID'] == 4)

        if n_matched_fiducials < 3:
            log.error('Fewer than three matched fiducials; exiting early.')
            return 13

        tx = FVC2FP.read_jsonfile(fvc2fp_filename())

        if args.use_spotmatch :
            metrology = load_metrology()
            # keep only the center of fiducials
            selection = (metrology["DEVICE_TYPE"]=="FIF")|(metrology["DEVICE_TYPE"]=="GIF")
            for loc in np.unique(metrology["LOCATION"][selection]) :
                # all the pinholes at that location
                ii = np.where(metrology["LOCATION"]==loc)[0]
                # replace first entry by mean of pinholes
                metrology["X_FP"][ii[0]] = np.mean(metrology["X_FP"][ii])
                metrology["Y_FP"][ii[0]] = np.mean(metrology["Y_FP"][ii])
                # set a dummy pinhole id = 10 just to make sure it's not interpreted as an existing pinhole
                metrology["PINHOLE_ID"][ii[0]] = 99
                # drop the others
                metrology.remove_rows(ii[1:])
        else :
            metrology = None

        tx.fit(spots, metrology=metrology, update_spots=True, zbfit=(args.zbfit), fixed_scale=args.fixed_scale, fixed_rotation=args.fixed_rotation)

        expected_pos = get_expected_pos(args, log)

        # select spots that are not already matched
        selection  = (spots["LOCATION"]==-1)

        if args.use_spotmatch :
            spots = spotmatch(spots["XPIX"],spots["YPIX"],expected_x_fp=expected_pos["X_FP"],expected_y_fp=expected_pos["Y_FP"],expected_location=expected_pos["LOCATION"],fvc2fp=tx,match_radius_pixels=args.spotmatch_match_radius_pixels)
            #spots = spotmatch(spots["XPIX"],spots["YPIX"],expected_x_fp=expected_pos["X_FP"],expected_y_fp=expected_pos["Y_FP"],expected_location=expected_pos["LOCATION"],fvc2fp=None,match_radius_pixels=args.spotmatch_match_radius_pixels)


            # add info
            nspots = len(spots["LOCATION"])
            for k in ["X_FP","Y_FP","X_FP_EXP","Y_FP_EXP","X_FP_METRO","Y_FP_METRO"] :
                if k not in spots.dtype.names :
                    spots[k] = np.zeros(nspots,dtype=float)

            spots["X_FP"],spots["Y_FP"] = tx.fvc2fp(spots["XPIX"],spots["YPIX"])

            is_matched = (spots["LOCATION"]>=0)
            loc2i = {loc:i for i,loc in enumerate(spots["LOCATION"])}

            ii=[]
            jj=[]
            for j,loc in enumerate(expected_pos["LOCATION"]) :
                if loc in loc2i :
                    ii.append(loc2i[loc])
                    jj.append(j)
            spots["X_FP_EXP"][ii]=expected_pos["X_FP"][jj]
            spots["Y_FP_EXP"][ii]=expected_pos["Y_FP"][jj]

            metrology = load_metrology()
            ii=[]
            jj=[]
            for j,loc in enumerate(metrology["LOCATION"]) :
                if loc in loc2i :
                    ii.append(loc2i[loc])
                    jj.append(j)
            spots["X_FP_METRO"][ii]=metrology["X_FP"][jj]
            spots["Y_FP_METRO"][ii]=metrology["Y_FP"][jj]


        else :

            # match
            indices_of_expected_pos,distances = match_same_system(spots["X_FP"][selection],spots["Y_FP"][selection],expected_pos["X_FP"],expected_pos["Y_FP"])
            is_matched = (distances<args.max_match_distance)&(indices_of_expected_pos>=0)
            ii=np.where(selection)[0]
            selection[ii]          &=  is_matched
            indices_of_expected_pos = indices_of_expected_pos[is_matched]
            distances               = distances[is_matched]

            # add columns after matching fibers
            for k1,k2 in zip(["X_FP","Y_FP"],["X_FP_EXP","Y_FP_EXP"]) :
                if k2 not in spots.keys() : spots[k2] = np.zeros(len(spots))
                spots[k2][selection]=expected_pos[k1][indices_of_expected_pos]
            for k in ["EXP_Q_0","EXP_S_0","PETAL_LOC","DEVICE_LOC","DEVICE_ID","DEVICE_TYPE","LOCATION"] :
                if k in expected_pos.keys() :
                    if k not in spots.keys() :
                        if k in ["DEVICE_ID","DEVICE_TYPE"] :
                            spots[k] = np.repeat("None           ",len(spots))
                        else :
                            spots[k] = np.zeros(len(spots))
                    spots[k][selection]=expected_pos[k][indices_of_expected_pos]


        if args.expected_positions is not None and ( args.turbulence_correction or args.turbulence_correction_with_pol ):

            if args.turbulence_correction_with_pol :
                log.info("Turbulence correction (local polynomial fit) ...")
            else :
                log.info("Turbulence correction (gaussian processes) ...")
            selection=(spots["LOCATION"]>=0)&(spots["X_FP_EXP"]!=0)&(spots["Y_FP_EXP"]!=0) # matched

            if args.turbulence_correction_with_pol :
                new_x,new_y = correct_with_pol(spots["X_FP"][selection],spots["Y_FP"][selection],spots["X_FP_EXP"][selection],spots["Y_FP_EXP"][selection])
            else :
                new_x,new_y = correct(spots["X_FP"][selection],spots["Y_FP"][selection],spots["X_FP_EXP"][selection],spots["Y_FP_EXP"][selection])
            rms_before = np.sqrt(np.mean((spots["X_FP"][selection]-spots["X_FP_EXP"][selection])**2+(spots["Y_FP"][selection]-spots["Y_FP_EXP"][selection])**2))
            rms_after  = np.sqrt(np.mean((new_x-spots["X_FP_EXP"][selection])**2+(new_y-spots["Y_FP_EXP"][selection])**2))
            log.info("rms(measured-expected)={:5.4f}mm rms(corrected-expected)={:5.4f}mm ".format(rms_before,rms_after))
            if rms_before < rms_after :
                log.warning("turbulence correction does not reduce the rms?")
                # apply it anyway because there might be a good reason for this
            spots["X_FP"][selection] = new_x
            spots["Y_FP"][selection] = new_y

        if "X_FP_METRO" in spots.dtype.names :
            # for spots with metrology X_FP_EXP=X_FP_METRO
            selection = (spots["X_FP_METRO"]!=0)
            spots["X_FP_EXP"][selection]=spots["X_FP_METRO"][selection]
            selection = (spots["Y_FP_METRO"]!=0)
            spots["Y_FP_EXP"][selection]=spots["Y_FP_METRO"][selection]

        # write transfo
        if args.output_transform is not None :
            if not args.output_transform.endswith(".json") :
                print("error, can only write json files, so please choose an output filename end ing with .json")
            else :
                tx.write_jsonfile(args.output_transform)
                print("wrote transform in {}".format(args.output_transform))

        if args.field_model is not None :
            log.info("Reading field model in {}".format(args.field_model))
            with open(args.field_model) as file :
                fm = FieldModel.fromjson(file.read())
            spots["RA"] = np.zeros(len(spots),dtype=float)
            spots["DEC"] = np.zeros(len(spots),dtype=float)
            ii=(spots["X_FP"]!=0)&(spots["Y_FP"]!=0)
            ra,dec = fm.fp2radec(spots["X_FP"][ii],spots["Y_FP"][ii])
            spots["RA"][ii]  = ra
            spots["DEC"][ii] = dec

        # write spots
        spots.write(outfile,format="csv",overwrite=True)
        print("wrote {}".format(outfile))

    return 0


def get_parser():
    """Define FVC processing options and defaults.
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                        description="""FVC image processing""")
    parser.add_argument('-i','--infile', type = str, default = None, required = True,
                        help = 'path to FVC image fits file or CSV file with spots positions')
    parser.add_argument('-o','--outfile', type = str, default = None, required = True,
                        help = 'path to output CSV ASCII file')
    parser.add_argument('--extname', type = str, default = 'F0000', required = False,
                        help = 'input EXTNAME to use if more than one HDU; last for last HDU')
    parser.add_argument('--sequence', action = 'store_true',
                        help = 'loop over FXXXX extensions, will write one output file per extension replacing .csv by -XXXX.csv')
    parser.add_argument('--output-transform', type = str, default = None, required = False,
                        help = 'write transformation to this json file')
    parser.add_argument('--input-transform', type = str, default = None, required = False,
                        help = 'use this json file as input for the match, defaut is data/default-fvc2fp.json')
    parser.add_argument('--pinhole-max-separation-mm', type = float, default = 1.5, required = False,
                        help = 'pinhole maximum separation in pixels')
    parser.add_argument('--max-match-distance', type = float, default = 7, required = False,
                        help = 'maximum match distance in mm')
    parser.add_argument('--min-counts-per-pixel', type = float, default = 500., required = False,
                        help = "threshold for spots detection in counts/pixel after convolution")
    parser.add_argument('--threshold', type = float, default = None, required = False,
                        help = "(deprecated)")
    parser.add_argument('--min-counts', type = float, default = 10000., required = False,
                        help = "(threshold for spots detection in counts (aperture flux around centroids)")
    parser.add_argument('--expected-positions', type = str, default = None, required = False,
                        help = 'path to fits file or CSV file with fibers expected position (table column EXP_Q and EXP_S')
    parser.add_argument('--hdu-for-expected-positions', type = str, default = "DATA", required = False,
                        help = 'specify HDU in expected position fits file')
    parser.add_argument('--field-model', type = str, default = None, required = False,
                        help = 'use this json file as field model to convert fiber coordinates to RA Dec')
    parser.add_argument('--fvc-psf-sigma', type = float, default = 1., required = False,
                        help = 'use this PSF sigma for the centroid measurement')
    parser.add_argument('--nomatch', action = 'store_true',
                        help = 'do not try to match data, just write FVC coordinates of detected spots')
    parser.add_argument('--min-spots', type = int , default = None, required = False,
                        help = 'minimum number of detected spots to run')
    parser.add_argument('--max-spots', type = int , default = None, required = False,
                        help = 'max number of detected spots to run')
    parser.add_argument('--no-zbfit', action = 'store_true',
                        help = '(deprecated)')
    parser.add_argument('--zbfit', action = 'store_true',
                        help = 'refit the high degree ZB polynomial coefficients')
    parser.add_argument('--fixed-scale', action = 'store_true',
                        help = 'do not refit the scale')
    parser.add_argument('--fixed-rotation', action = 'store_true',
                        help = 'do not refit the rotation')
    parser.add_argument('--bias', type = str, default = None, required = False,
                        help = 'use this master bias fits file instead of $DESIMETER_DATA/bias.fits if exist')
    parser.add_argument('--no-bias', action = 'store_true',
                        help = 'do not apply a bias correction (by default, use $DESIMETER_DATA/bias.fits if exist')
    parser.add_argument('--use-spotmatch', action = 'store_true',
                        help = 'use spotmatch instead of desimeter code for matching positioners and fiducials (requires the program match_positions installed and in the PATH)')
    parser.add_argument('--turbulence-correction', action = 'store_true',
                        help = 'turbulence correction assuming offsets between measured and expected coordinates are not spatially correlated. Only used in conjunction with option --expected-positions.')
    parser.add_argument('--turbulence-correction-with-pol', action = 'store_true',
                        help = 'turbulence correction using local polynomial fit instead of Gaussian processes.')
    parser.add_argument('--spotmatch-match-radius-pixels', type = float, default = 70, required = False, help = "spotmatch matching radius")
    parser.add_argument('--make-directory', action='store_true',
                        help='Make directory for output if needed.')
    return parser
