"""
Wrapper to the C program 'match_positions' from David Rabinowitz.
"""



import os
import numpy as np
import subprocess
import shutil
import tempfile

from astropy.table import Table

from desimeter.io import load_metrology,fvc2fp_filename
from desimeter.transform.fvc2fp import FVC2FP

def _compute_pixel_scale(fvc2fp) :
    """
    computes square root of determinant of jacobian to get the pixel scale in mm
    returns pixel scale
    """
    xc=3000.
    yc=3000.
    eps=0.1
    xpix=np.array([xc,xc+eps,xc])
    ypix=np.array([yc,yc,yc+eps])
    xfp,yfp = fvc2fp.fvc2fp(xpix,ypix)
    J=np.zeros((2,2))
    J[0,0]=(xfp[1]-xfp[0])/eps
    J[0,1]=(xfp[2]-xfp[0])/eps
    J[1,0]=(yfp[1]-yfp[0])/eps
    J[1,1]=(yfp[2]-yfp[0])/eps
    return np.sqrt(np.abs(np.linalg.det(J))) # mm per pixel


def _write_spotmatch_fiducial_config_file(filename) :


    """
    example:


    """
    metrology = load_metrology()

    selection  = (metrology["DEVICE_TYPE"]=="FIF")|(metrology["DEVICE_TYPE"]=="GIF")
    locations = np.unique(metrology["LOCATION"][selection])

    with open(filename,"w") as ofile :

        for location in locations :
            selection = (metrology["LOCATION"]==location)
            pinid     = metrology["PINHOLE_ID"][selection]
            xfp       = metrology["X_FP"][selection]
            yfp       = metrology["Y_FP"][selection]
            mx = np.mean(xfp)
            my = np.mean(yfp)
            dx = xfp-mx
            dy = yfp-my

            ofile.write("{} {:4.3f} {:4.3f} 000.0 0.001 {:05d} 8\n".format(location,0.,0.,0))
            for i in range(dx.size) :
                ofile.write("{} {:4.3f} {:4.3f} 000.0 0.001 {:05d} 2\n".format(location,dx[i],dy[i],pinid[i]))
    print("wrote",filename)

def _write_spotmatch_targets_file(x_fp,y_fp,location,filename,fvc2fp=None) :
    """
    example:

1000 179.640 5616.148  12.000   0.001 4
1001 336.988 5615.164  12.000   0.001 4
1002 259.003 5479.802  12.000   0.001 4
1003 493.786 5613.850  12.000   0.001 4
1004 415.894 5478.573  12.000   0.001 4
1005 650.448 5612.610  12.000   0.001 4
1006 572.602 5477.142  12.000   0.001 4
...
    """
    if fvc2fp is None:
        fvc2fp = FVC2FP.read(fvc2fp_filename())
        print("use default fvc2fp")

    x_fp = np.array(x_fp)
    y_fp = np.array(y_fp)
    location = np.array(location)
    flags = np.repeat(4,x_fp.size)

    # add fiducials centers
    metrology = load_metrology()
    selection  = (metrology["DEVICE_TYPE"]=="FIF")|(metrology["DEVICE_TYPE"]=="GIF")
    fid_locations = np.unique(metrology["LOCATION"][selection])
    for fid_location in fid_locations :
        selection = (metrology["LOCATION"]==fid_location)
        mx = np.mean(metrology["X_FP"][selection])
        my = np.mean(metrology["Y_FP"][selection])
        x_fp = np.append(x_fp,mx)
        y_fp = np.append(y_fp,my)
        location = np.append(location,fid_location)
        flags = np.append(flags,8)


    xpix,ypix = fvc2fp.fp2fvc(x_fp,y_fp)

    with open(filename,"w") as ofile :
        for i in range(xpix.size) :
            ofile.write("{} {:4.3f} {:4.3f} 12.000 0.001 {}\n".format(location[i],xpix[i],ypix[i],flags[i]))
    print("wrote",filename)



def _write_spotmatch_measured_pos_file(xpix,ypix,filename):
    """
    example:

3495.648 2544.892  16.170 1   3.101
4926.980 2742.090  12.956 2   2.388
4169.921 2742.754  12.985 3   2.271
4358.021 2747.092  12.928 4   2.346
4062.197 2792.177  12.908 5   2.293
3925.400 2874.306  13.002 6   2.067
4857.292 2876.005  12.943 7   2.323

    """

    assert(len(xpix)==len(ypix))

    with open(filename,"w") as ofile :
        for i in range(len(xpix)) :
            ofile.write("{:4.3f} {:4.3f} 12.000 {} 2.000\n".format(xpix[i],ypix[i],i+1))

    print("wrote",filename)

def _write_spotmatch_reference_pos_file(filename,fvc2fp=None) :
    """
    example:

 1541 4179.496 2557.798 13.670 515 1.897 pinhole
 1541 4177.541 2572.562 17.823 1027 1.090 pinhole
 1541 4191.178 2578.614 13.016 259 2.055 pinhole
 1541 4175.588 2587.326 12.784 771 1.851 pinhole
 1542 5090.489 2633.933 13.025 771 2.115 pinhole
 1542 5099.186 2646.061 12.641 1027 2.115 pinhole
    """

    metrology = load_metrology()
    if fvc2fp is None:
        fvc2fp = FVC2FP.read(fvc2fp_filename())

    selection = (metrology["DEVICE_TYPE"]=="FIF")|(metrology["DEVICE_TYPE"]=="GIF")
    xpix,ypix   = fvc2fp.fp2fvc(metrology["X_FP"][selection],metrology["Y_FP"][selection])
    locations  = metrology["LOCATION"][selection]
    pinhole_ids = metrology["PINHOLE_ID"][selection]

    num={1:515,2:1027,3:259,4:771} # the code wants that

    with open(filename,"w") as ofile :
        for i in range(xpix.size) :
            ofile.write(" {} {:4.3f} {:4.3f} 12.000 {} 2.000 pinhole\n".format(locations[i],xpix[i],ypix[i],num[pinhole_ids[i]]))

    print("wrote",filename)

def spotmatch(xpix,ypix,expected_x_fp=None,expected_y_fp=None,expected_location=None,verbose=0,match_radius_pixels=70) :
    """

    """

    if shutil.which("match_positions") is None :
        raise RuntimeError("match_positions is not in PATH. You need to install spotmatch first. It's in https://desi.lbl.gov/trac/browser/code/online/FVC/spotmatch.")


    image_rows=6000
    image_cols=6000
    target_x_dir=1
    target_y_dir=1
    target_x0=0
    target_y0=0
    fid_x_dir=-1
    fid_y_dir=1

    tmp_dir=tempfile.gettempdir()

    fvc2fp = FVC2FP.read(fvc2fp_filename())
    exp_pixel_scale = _compute_pixel_scale(fvc2fp)


    fiducial_config_filename=os.path.join(tmp_dir,"desimeter_spotmatch_fiducials.txt")
    _write_spotmatch_fiducial_config_file(fiducial_config_filename)

    if expected_x_fp is None :
        print("input expected = None, so we use the centers of positioners as the expected spots")
        metrology = load_metrology()
        spots = metrology[metrology["DEVICE_TYPE"]=="POS"] # select positioners
        expected_x_fp = spots["X_FP"]
        expected_y_fp = spots["Y_FP"]
        expected_location = spots["LOCATION"]


    targets_filename=os.path.join(tmp_dir,"desimeter_spotmatch_targets.txt")
    _write_spotmatch_targets_file(expected_x_fp,expected_y_fp,expected_location,targets_filename,fvc2fp)

    measured_pos_filename=os.path.join(tmp_dir,"desimeter_spotmatch_input_centroids.txt")
    _write_spotmatch_measured_pos_file(xpix,ypix,measured_pos_filename)

    reference_pos_filename=os.path.join(tmp_dir,"desimeter_spotmatch_pinhole_references.txt")
    _write_spotmatch_reference_pos_file(reference_pos_filename,fvc2fp=None)

    """
match_positions -verbose 1 -image_rows 6000 -image_cols 6000 -target_x_dir 1 -t
arget_y_dir 1 -target_x0 0.0 -target_y0 0.0 -fid_x_dir -1 -fid_y_dir 1 -exp_pixel_scale 0.073 -match_radius 50 -redu
ced_pos_file ./match_centers.tmp -fiducial_config_file ./match_fiducials.tmp -target_pos_file ./match_targets.tmp -m
easured_pos_file ./match_centroids.tmp -pos_save_file ./measured_pos.tmp -reference_pos_file ./pinhole_references_20
200410.dat


    """

    match_centers_filename=os.path.join(tmp_dir,"desimeter_spotmatch_output_centers.txt")
    saved_pos_filename=os.path.join(tmp_dir,"desimeter_spotmatch_saved_output.txt")

    cmd =  "match_positions"
    cmd += " -verbose {}".format(verbose)
    cmd += " -image_rows {}".format(image_rows)
    cmd += " -image_cols {}".format(image_cols)
    cmd += " -target_x_dir {}".format(target_x_dir)
    cmd += " -target_y_dir {}".format(target_y_dir)
    cmd += " -fid_x_dir {}".format(fid_x_dir)
    cmd += " -fid_y_dir {}".format(fid_y_dir)
    cmd += " -exp_pixel_scale {}".format(exp_pixel_scale)
    cmd += " -match_radius {}".format(match_radius_pixels)
    cmd += " -reduced_pos_file {}".format(match_centers_filename)
    cmd += " -fiducial_config_file {}".format(fiducial_config_filename)
    cmd += " -target_pos_file {}".format(targets_filename)
    cmd += " -measured_pos_file {}".format(measured_pos_filename)
    cmd += " -pos_save_file {}".format(saved_pos_filename)
    cmd += " -reference_pos_file {}".format(reference_pos_filename)


    print(cmd)

    subprocess.call(cmd.split(" "))



    location=[]
    xpix=[]
    ypix=[]
    mag=[]
    flag=[]
    err=[]
    device_type=[]
    print("reading",match_centers_filename)
    with open(match_centers_filename) as ifile :
        for line in ifile.readlines() :
            vals=line.strip().split()
            if len(vals)!=7 : continue
            location.append(int(vals[0]))
            xpix.append(float(vals[1]))
            ypix.append(float(vals[2]))
            mag.append(float(vals[3]))
            flag.append(float(vals[4]))
            err.append(float(vals[5]))
            device_type.append(vals[6])
    print("read {} entries in {}".format(len(xpix),match_centers_filename))

    res=Table()
    res["LOCATION"]=location
    res["XPIX"]=xpix
    res["YPIX"]=ypix
    res["FLAG"]=flag
    res["ERR"]=err
    nspots=len(location)
    res["SPOTMATCH_DEVICE_TYPE"]=device_type

    res["DEVICE_ID"]   = np.repeat("unknown",nspots)
    res["DEVICE_TYPE"] = np.repeat("unknown",nspots)

    # use metrology to get more info
    metrology = load_metrology()
    locmap = {loc:index for index,loc in enumerate(metrology["LOCATION"])}
    matched_index   = np.where(res["LOCATION"]>=0)[0]
    metrology_index = [locmap[loc] for loc in res["LOCATION"][matched_index]]
    keys_to_copy = []
    for k in ["DEVICE_TYPE","DEVICE_LOC","DEVICE_ID","PETAL_LOC","PETAL_ID","BUS_ID"] :
        if k in metrology.dtype.names :
            keys_to_copy.append(k)
            if not k in res.dtype.names :
                res[k]=np.repeat(-1,nspots)
            res[k][matched_index] = metrology[k][metrology_index]

    # special treatment for pinhole ids
    res["PINHOLE_ID"]  = np.repeat(0,nspots)
    fiducials=(res["SPOTMATCH_DEVICE_TYPE"]=="fiducial")

    # dummy pinhole numbers. 99 for center, 10,11,12,13 for pinholes
    res["PINHOLE_ID"][fiducials] = 99
    for loc in np.unique(res["LOCATION"][fiducials]) :
        ii = (res["LOCATION"]==loc)&(res["SPOTMATCH_DEVICE_TYPE"]=="pinhole")
        res["PINHOLE_ID"][ii] = 10+np.arange(np.sum(ii))

    return res
