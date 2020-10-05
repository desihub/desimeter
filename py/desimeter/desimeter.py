import os

from desimeter.log import get_logger
from desimeter.detectspots import detectspots
from desimeter.findfiducials import findfiducials
from desimeter.transform.fvc2fp import FVC2FP
from desimeter.io import fvc2fp_filename, load_metrology
from desimeter.match import match_same_system
from desimeter.circles import fit_circle
from desimeter.transform.xy2qs import xy2uv, uv2xy
from desimeter.io import desimeter_data_dir
from desimeter.transform.zhaoburge import fit_scale_rotation_offset
from desimeter.transform.gfa2fp import fit_gfa2fp
from desimeter.fieldmodel import FieldModel        
from desimeter.time import mjd2lst

from astropy.table import Table

import fitsio
import numpy as np

class Desimeter(object):
    def __init__(self, desimeter_dir=None, data_dir=None,
                 gfa_dir=None,
                 proc_data_dir=None):
        '''
        desimeter_dir: directory containing metrology and other Desimeter state files.
        data_dir: directory containing data from the DESI spectrographs & FVC.
        proc_data_dir: directory containing intermediate processed results.
        '''
        if desimeter_dir is None:
            # will check for $DESIMETER_DATA, else use desimeter/data in desimeter repo.
            desimeter_dir = desimeter_data_dir()
        self.desimeter_dir = desimeter_dir

        if data_dir is None:
            # check $DESI_DATA_DIR, else assume "data" in current dir?
            data_dir = os.environ.get('DESI_DATA_DIR', 'data')
        self.data_dir = data_dir

        if gfa_dir is None:
            # check $DESI_GFA_DIR, else 'data'
            gfa_dir = os.environ.get('DESI_GFA_DIR', 'data')
        self.gfa_dir = gfa_dir

        fn = self.find_file('fvc2fp')
        self.fvc2fp = FVC2FP.read_jsonfile(fn)
        fn = self.find_file('metrology')
        self.metro = Table.read(fn)

        if not "LOCATION" in self.metro.keys() :
            # add useful location keyword
            self.metro["LOCATION"] = np.array(self.metro["PETAL_LOC"])*1000+np.array(self.metro["DEVICE_LOC"])

        # Set up GFA-pixel-to-FP-mm transforms
        self.gfa_trans = fit_gfa2fp(self.metro)

        if proc_data_dir is None:
            proc_data_dir = '.'
        self.proc_dir = proc_data_dir

    def find_file(self, filetype, night=None, expnum=None, frame=None,
                  tag=None,
                  desimeter_dir=None):
        from glob import glob

        if desimeter_dir is None:
            desimeter_dir = self.desimeter_dir

        if filetype in ['fvc2fp', 'metrology']:
            fn = { 'fvc2fp': 'init-fvc2fp.json',
                   'metrology': 'fp-metrology.csv',
            }[filetype]
            return os.path.join(desimeter_dir, fn)

        if filetype in ['fvc', 'guide', 'coordinates']:
            suff = dict(coordinates='').get(filetype, '.fz')
            if night is None:
                pat = os.path.join(self.data_dir, '*', '%08d' % expnum,
                                   '%s-%08d.fits%s' % (filetype, expnum, suff))
                print('Checking', pat)
                fns = glob(pat)
                if len(fns) == 0:
                    return None
                return fns[0]
            fn = os.path.join(self.data_dir, night, '%s-%08d.fits.fz' % (filetype, expnum))
            return fn
        if filetype == 'fvc-spots':
            fn = os.path.join(self.proc_dir,
                              'fvc-spots-%08d-F%04d.fits' % (expnum, frame))
            return fn
        if filetype == 'fvc-circles':
            fn = os.path.join(self.proc_dir,
                              'fvc-circles-%08d-%s.fits' % (expnum, tag))
            return fn
        if filetype == 'fieldmodel':
            fn = os.path.join(self.proc_dir,
                              'fm-%08d-F%04d.json' % (expnum, frame))
            return fn
        if filetype == 'guide-catalog':
            pat = os.path.join(self.gfa_dir, '*', '%08d' % expnum,
                               'guide-%08d_catalog-%05d.fits' % (expnum, frame))
            print('Checking', pat)
            fns = glob(pat)
            if len(fns) == 0:
                return None
            return fns[0]
            
        raise RuntimeError('Unknown file type "%s"' % filetype)

    def write_desimeter(self, dirname,
                        metrology=True,
                        fvc2fp=True):
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        outfn = self.find_file('fvc2fp', desimeter_dir=dirname)
        # fake overwriting an existing file
        if os.path.exists(outfn):
            os.remove(outfn)
        self.fvc2fp.write(outfn)
        print('Wrote', outfn)
        outfn = self.find_file('metrology', desimeter_dir=dirname)
        if os.path.exists(outfn):
            os.remove(outfn)
        self.metro.write(outfn)
        print('Wrote', outfn)

    def measure_spots(self, expnum, frame):
        infn = self.find_file('fvc', expnum=expnum)
        F = fitsio.read(infn, ext='F%04i' % frame)

        spots = detectspots(F, nsig=7, psf_sigma=1.)
        spots = findfiducials(spots,
                              input_transform_func=self.fvc2fp,
                              metrology=self.metro,
                              pinhole_max_separation_mm=1.5)
        spots = self.refit_spots(spots)
        return spots

    def refit_spots(self, spots, zbfit=True):
        self.fvc2fp.fit(spots, metrology=self.metro,
                        update_spots=True, zbfit=zbfit,
                        fixed_scale=False, fixed_rotation=False)

        # select spots that are not already matched
        selection  = (spots["LOCATION"]==-1)

        # match
        indices_of_expected_pos,distances = match_same_system(spots["X_FP"][selection],spots["Y_FP"][selection],self.metro["X_FP"],self.metro["Y_FP"])

        max_match_distance = 7
        is_matched = (distances<max_match_distance)&(indices_of_expected_pos>=0)
        ii=np.where(selection)[0]
        selection[ii]          &=  is_matched
        indices_of_expected_pos = indices_of_expected_pos[is_matched]
        distances               = distances[is_matched]

        # add columns after matching fibers
        for k1,k2 in zip(["X_FP","Y_FP"],["X_FP_EXP","Y_FP_EXP"]) :
            if k2 not in spots.keys() : spots[k2] = np.zeros(len(spots))
            spots[k2][selection]=self.metro[k1][indices_of_expected_pos]
        for k in ["EXP_Q_0","EXP_S_0","PETAL_LOC","DEVICE_LOC","LOCATION"] :
            if k in self.metro.keys() :
                if k not in spots.keys() : spots[k] = np.zeros(len(spots))
                spots[k][selection]=self.metro[k][indices_of_expected_pos]

        # for spots with metrology X_FP_EXP=X_FP_METRO
        selection = (spots["X_FP_METRO"]!=0)
        spots["X_FP_EXP"][selection]=spots["X_FP_METRO"][selection]
        selection = (spots["Y_FP_METRO"]!=0)
        spots["Y_FP_EXP"][selection]=spots["Y_FP_METRO"][selection]

        # Lots of those with poor matches are bad spots -- eg cosmic rays or something
        #dist = np.hypot(spots['X_FP_EXP'] - spots['X_FP'], spots['Y_FP_EXP'] - spots['Y_FP'])
        #bad = np.flatnonzero(dist > 50)
        return spots

    def measure_circles(self, allspots, do_plot=False, nmaxplot=10):
        x={}
        y={}
        xexp={}
        yexp={}
        first=True
        for spots in allspots:
            t = spots
            selection=(spots["LOCATION"]>0)
            location_and_pinhole=(np.array(t["LOCATION"])*10+np.array(t["PINHOLE_ID"])).astype(int)
            if first :
                for loc in location_and_pinhole[selection] :
                    x[loc] = []
                    y[loc] = []
                    xexp[loc] = float(t["X_FP_EXP"][location_and_pinhole==loc][0])
                    yexp[loc] = float(t["Y_FP_EXP"][location_and_pinhole==loc][0])
                    #print(loc,xexp[loc],yexp[loc])
                first=False

            for loc in location_and_pinhole[selection] :
                ii = np.where(location_and_pinhole==loc)[0]
                if ii.size > 1 :
                    print("several matched for LOCATION ",loc)
                    continue
                i=ii[0]
                if not loc in x.keys() :
                    x[loc] = []
                    y[loc] = []
                    xexp[loc] = float(t["X_FP_EXP"][location_and_pinhole==loc][0])
                    yexp[loc] = float(t["Y_FP_EXP"][location_and_pinhole==loc][0])
                x[loc].append(float(t["X_FP"][i]))
                y[loc].append(float(t["Y_FP"][i]))

        location_and_pinhole=np.array(list(x.keys()),dtype=int)
        location=location_and_pinhole//10
        pinhole=location_and_pinhole%10
        print("number of positioners:",np.sum(pinhole==0))
        print("number of fiducials:",np.sum(pinhole==1))
        print("number of pinholes:",np.sum(pinhole>=1))
        ndots=len(location_and_pinhole)
        theta=np.linspace(0,2*np.pi,50)
        xfp_metro=np.zeros(ndots)
        yfp_metro=np.zeros(ndots)
        xfp_meas=np.zeros(ndots)
        yfp_meas=np.zeros(ndots)
        count=0
        for iloc,loc in enumerate(x.keys()) :
            if len(x[loc])<6 : continue
            x[loc]=np.array(x[loc])
            y[loc]=np.array(y[loc])
            ii=np.where(x[loc]!=0)[0]
            x[loc]=x[loc][ii]
            y[loc]=y[loc][ii]
            if pinhole[iloc] == 0 and np.std(x[loc])<1. :
                # this is a non-moving positioner, I don't use this
                continue
            count += 1

            xc=np.median(x[loc])
            yc=np.median(y[loc])

            if pinhole[iloc] == 0 : # it's a positioner
                # here is the fit
                try:
                    #- Transform to curved focal surface which is closer to a real circle
                    x_cfs, y_cfs = xy2uv(x[loc], y[loc])
                    #- Do the fit
                    xc_cfs,yc_cfs,r = fit_circle(x_cfs, y_cfs)
                    #- Convert center back into CS5 x,y
                    xc, yc = uv2xy(xc_cfs, yc_cfs)
                except ValueError:
                    print("fit circle failed for loc={} x={} y={}".format(loc,xc,yc))
                    continue

                if iloc%100==0 :
                    print("{}/{} loc={} x={} y={} r={}".format(iloc,len(x),loc,xc,yc,r))
                if r<0.1 : continue

                if do_plot and count<nmaxplot :
                    plt.figure("circles")
                    plt.plot(x[loc],y[loc],"o")
                    plt.plot(xexp[loc],yexp[loc],"x")
                    theta=np.linspace(0,2*np.pi,50)
                    plt.plot(xc+r*np.cos(theta),yc+r*np.sin(theta),"-",color="green")
                    plt.plot(xc,yc,"+",color="green")
            xfp_metro[iloc]=xexp[loc]
            yfp_metro[iloc]=yexp[loc]
            xfp_meas[iloc]=xc
            yfp_meas[iloc]=yc
        dx=xfp_meas-xfp_metro
        dy=yfp_meas-yfp_metro
        dr=np.sqrt(dx**2+dy**2)
        print("median offset = %4.1f um" % (np.median(dr[dr!=0])*1000.));
        ii=np.where((xfp_metro!=0)&(dr<3.))[0]

        # make a table out of that
        t2=Table([location[ii],pinhole[ii],xfp_metro[ii],yfp_metro[ii],xfp_meas[ii],yfp_meas[ii]],names=["LOCATION","PINHOLE_ID","X_FP_METRO","Y_FP_METRO","X_FP","Y_FP"],dtype=[int,int,float,float,float,float])
        return t2

    def refit_zb(self, x1, y1, x2, y2, zbpolids):
        '''
        x1,y1: measured {X,Y}_FP
        x2,y2: metrology {X,Y}_FP
        zbpolids: list of integers of Zhao-Burge polynomial numbers
        '''
        transfo = self.fvc2fp

        # apply transfo back to FVC pixels
        xpix,ypix = transfo.fp2fvc(x1,y1)

        # set polynomial ids
        transfo.zbpolids = zbpolids
        # and redo the fit, now globally
        rxpix, rypix = transfo._reduce_xyfvc(xpix,ypix)
        rxfp, ryfp = transfo._reduce_xyfp(x2,y2)
        scale, rotation, offset_x, offset_y, zbpolids, zbcoeffs = fit_scale_rotation_offset(
            rxpix, rypix, rxfp, ryfp,
            fitzb=True, zbpolids=transfo.zbpolids, zbcoeffs=transfo.zbcoeffs)
        transfo.scale = scale
        transfo.rotation = rotation
        transfo.offset_x = offset_x
        transfo.offset_y = offset_y
        transfo.zbpolids = zbpolids
        transfo.zbcoeffs = zbcoeffs
        # and apply it now
        x1b,y1b = transfo.fvc2fp(xpix,ypix)
        dist=np.sqrt((x1b-x2)**2+(y1b-y2)**2)
        ok=(dist<0.08)

        scale, rotation, offset_x, offset_y, zbpolids, zbcoeffs = fit_scale_rotation_offset(rxpix[ok], rypix[ok], rxfp[ok], ryfp[ok], fitzb=True, zbpolids=transfo.zbpolids, zbcoeffs=transfo.zbcoeffs)
        transfo.scale = scale
        transfo.rotation = rotation
        transfo.offset_x = offset_x
        transfo.offset_y = offset_y
        transfo.zbpolids = zbpolids
        transfo.zbcoeffs = zbcoeffs
        # and apply it now
        x1b,y1b = transfo.fvc2fp(xpix,ypix)

        return x1b, y1b

    def fit_guide_star_coords(self, expnum, frame):
        log = get_logger()
        hdrfn = self.find_file('guide', expnum=expnum)
        header = fitsio.read_header(hdrfn)
        #if not "TARGTRA" in header :
        #    log.warning("no TARGTRA in header of HDU 0, try HDU 1")
        #    header = fitsio.read_header(args.fits_header,1)
        if not "TARGTRA" in header :
            log.error("no TARGTRA in header of file {}".format(hdrfn))
            return None

        fm = FieldModel()
        fm.ra  = header["TARGTRA"]
        fm.dec = header["TARGTDEC"]
        fm.expid = header["EXPID"]
        fm.hexrot_deg = float(header["FOCUS"][5])/3600.
        fm.adc1 = header["ADC1PHI"]
        fm.adc2 = header["ADC2PHI"]

        catfn = self.find_file('guide-catalog', expnum=expnum, frame=frame)
        catalog = fm.read_guide_stars_catalog(catfn)

        # MJD and LST are needed for sky transform
        fm.mjd    = np.mean(catalog["mjd_obs"])
        fm.lst    = mjd2lst(fm.mjd)
        fm.fit_tancorr(catalog, gfa_transform=self.gfa_trans)
        return fm
