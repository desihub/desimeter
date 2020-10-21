import os
import numpy as np
import pylab as plt
import astropy.table
from desimeter.desimeter import Desimeter
import fitsio

def outline_petals(petals=None, **kwargs):
    '''
    kwargs: matplotlib plot() kwargs.
    '''
    th = np.linspace(0,np.pi/5,100)
    if petals is None:
        petals = list(range(10))
    for i in petals:
        oth = np.deg2rad(i*36)
        r1 = 420
        r0 = 40
        plt.plot(r1*np.cos(oth+th), r1*np.sin(oth+th), '-', **kwargs)
        plt.plot(r0*np.cos(oth+th), r0*np.sin(oth+th), '-', **kwargs)
        plt.plot([r0*np.cos(oth+th[0]), r1*np.cos(oth+th[0])],
                 [r0*np.sin(oth+th[0]), r1*np.sin(oth+th[0])], '-', **kwargs)
        plt.plot([r0*np.cos(oth+th[-1]), r1*np.cos(oth+th[-1])],
                 [r0*np.sin(oth+th[-1]), r1*np.sin(oth+th[-1])], '-', **kwargs)

def plot_fiducial_offsets(table, expnum=-1, frame=-1):
    # Takes an "fvc-spots" table
    from desimeter.twrapper import Twrapper
    T = Twrapper(table)
    # Find fiducials (3-4 pinholes per fiducial)
    F = T[T.pinhole_id > 0]
    # Compute a "device id" that is unique across the focal plane
    F.devid = F.device_loc + F.petal_loc * 1000
    devs = np.unique(F.devid)

    # Average pinholes per device.
    F.dev_x    = np.zeros(len(F), np.float32)
    F.dev_y    = np.zeros(len(F), np.float32)
    F.dev_stdx = np.zeros(len(F), np.float32)
    F.dev_stdy = np.zeros(len(F), np.float32)
    F.dev_dx   = np.zeros(len(F), np.float32)
    F.dev_dy   = np.zeros(len(F), np.float32)
    idevs = []
    for d in devs:
        I = np.flatnonzero(F.devid == d)
        F.dev_x[I] = np.mean(F.x_fp[I])
        F.dev_y[I] = np.mean(F.y_fp[I])
        dx = F.x_fp[I] - F.x_fp_metro[I]
        dy = F.y_fp[I] - F.y_fp_metro[I]
        F.dev_dx[I] = np.mean(dx)
        F.dev_dy[I] = np.mean(dy)
        F.dev_stdx[I] = np.std(dx)
        F.dev_stdy[I] = np.std(dy)
        idevs.append(I[0])
        #plt.clf()
        #plt.plot(F.x_fp[I], F.y_fp[I], 'b.')
        #plt.plot(F.x_fp_metro[I], F.y_fp_metro[I], 'r.')
        #plt.show()
    D = F[np.array(idevs)]

    D.dev_radi = np.zeros(len(D), np.float32)
    D.dev_tang = np.zeros(len(D), np.float32)
    for i in range(len(D)):
        if D.device_loc[i] in [541, 542]:
            # decompose into radial and tangential terms.
            th = np.arctan2(D.dev_y[i], D.dev_x[i])
            rx =  np.cos(th)
            ry =  np.sin(th)
            tx = -np.sin(th)
            ty =  np.cos(th)
            radi = D.dev_dx[i] * rx + D.dev_dy[i] * ry
            tang = D.dev_dx[i] * tx + D.dev_dy[i] * ty
            print('GIF % 4i' % D.location[i], 'dx,dy',
                  '%+.3f, %+.3f' % (D.dev_dx[i], D.dev_dy[i]),
                  'Radial %+.3f, Tangential %+.3f' % (radi, tang))
            D.dev_radi[i] = radi
            D.dev_tang[i] = tang
    D.is_gif = np.logical_or(D.device_loc == 541, D.device_loc == 542)
    D.is_guide = np.isin(D.petal_loc, [0, 2, 3, 5, 7, 8])

    print('Average GUIDE GIF radial: %+.3f.  Average tangential: %+.3f' %
          (np.mean(D.dev_radi[D.is_gif * D.is_guide]),
           np.mean(D.dev_tang[D.is_gif * D.is_guide])))
    for p in np.unique(D.petal_loc):
        I = np.flatnonzero(D.petal_loc == p)
        plt.plot(D.dev_x[I], D.dev_y[I], 'o', mec='none', ms=25, alpha=0.1)
        #plt.text(np.mean(D.dev_x[I]), np.mean(D.dev_y[I]), 'Petal loc %i' % p)
        th = np.arctan2(np.mean(D.dev_y[I]), np.mean(D.dev_x[I]))
        pp = int(np.round((th / (2.*np.pi / 10.)) - 0.5))
        pth = (pp + 0.5) * (2.*np.pi/10.)
        R = 300.
        plt.text(np.cos(pth)*R, np.sin(pth)*R, 'Petal loc %i' % p, ha='center')
    #for x,y,d in zip(D.dev_x, D.dev_y, D.devid):
    #    plt.text(x, y, '%i' % d)
    qargs = dict(pivot='middle', angles='xy', scale_units='xy',
                 scale=0.0005)
    Q = plt.quiver(D.dev_x, D.dev_y, D.dev_dx, D.dev_dy, **qargs)
    plt.quiver(D.dev_x[D.is_gif], D.dev_y[D.is_gif], D.dev_dx[D.is_gif], D.dev_dy[D.is_gif], color='b', **qargs)
    # add quiver scale marker!
    sx = 20
    plt.quiverkey(Q, -400, 400, sx/1000., '%i um' % sx, coordinates='data')
    # Draw lines around the petals.
    outline_petals(color='k', alpha=0.1)
    plt.axis('equal')
    rms2d = 1000. * np.sqrt(np.mean(D.dev_dx**2 + D.dev_dy**2))
    #plt.title('Offsets of fiducials vs whole-focal-plane fit: %.1f um RMS2d, expnum %i frame %i' % (rms2d, expnum, frame))
    plt.xlabel('Focal plane x (mm)')
    plt.ylabel('Focal plane y (mm)')
    plt.title('Fiducial resids: %.1f um RMS2d, expnum %i frame %i' % (rms2d, expnum, frame))
    return D

def main():
    #plt.figure(figsize=(10,10))
    plt.figure(figsize=(6,6))

    dm = Desimeter(proc_data_dir='proc')

    fidfn = 'fiducials.fits'
    if not os.path.exists(fidfn):
        fids = []
        # All exposures on 2020-03-14,15
        for keep,expnum in [(False, 52644)] + [(True, e) for e in list(range(55353, 55692+1))]:
        #for expnum in [52644] + list(range(55353, 55357+1)):
            fvcfn = dm.find_file('fvc', expnum=expnum)
            if fvcfn is None:
                continue
            print('Expnum', expnum)
            print('FVC', fvcfn)
            frame = 0
            fn = dm.find_file('fvc-spots', expnum=expnum, frame=frame)
            if os.path.exists(fn):
                spots = astropy.table.Table.read(fn)
            else:
                spots = dm.measure_spots(expnum, frame)
                spots.write(fn, overwrite=True)
            spots = dm.refit_spots(spots)
            #print('Fit fvc2fp transform:', dm.fvc2fp.tojson())
            plt.clf()
            D = plot_fiducial_offsets(spots, expnum=expnum, frame=frame)
            fn = 'fvc-fid-offsets-%08i-F%04i.png' % (expnum, frame)
            plt.savefig(fn)
            print('Wrote', fn)
            plt.clf()
            if not keep:
                continue
            # Format table of per-fiducial results
            D.expnum = np.zeros(len(D), np.int32) + expnum
            D.frame  = np.zeros(len(D), np.int16) + frame
            hdr = fitsio.read_header(fvcfn)
            for key in ['MJD-OBS', 'TARGTRA', 'TARGTDEC', 'TARGTAZ', 'TARGTEZ', 'AIRMASS',
                        'ADC1PHI', 'ADC2PHI', 'TRUSTEMP']:
                D.set(key.replace('-','_').lower(), np.zeros(len(D), np.float32) + hdr.get(key, -99))
            fids.append(D._table)
        fids = astropy.table.vstack(fids)
        fids.write('fiducials.fits', overwrite=True)
    else:
        fids = astropy.table.Table.read(fidfn)

    from desimeter.twrapper import Twrapper
    fids = Twrapper(fids)
    print(len(fids), 'fiducial measurements,', len(np.unique(fids.devid)), 'unique fids')
    d = np.hypot(fids.dev_dx, fids.dev_dy)
    I = np.flatnonzero(d < 0.1)
    fids = fids[I]
    print('Cut to', len(fids), 'based on residuals')

    plt.clf()
    plt.hist(d * 1000., range=(0, 60), bins=30)
    rms = np.sqrt(np.mean((d*1000.)**2))
    plt.xlabel('Offset of fiducials (um)')
    plt.title('Before: Fiducials, 2020-03-14 & 15: RMS2D %.1f um' % rms)
    plt.savefig('fid-offset-hist-before.png')

    # Plots wrt some observing properties.
    X = fids[fids.devid == 1542]
    X.dx = X.dev_dx * 1000.
    X.dy = X.dev_dy * 1000.
    tt = 'Fiducial 1542, 2020-03-(14,15)'
    X.adc_dphi = np.fmod(360 + X.adc2phi-X.adc1phi, 360)

    for k in ('airmass', 'expnum', 'adc_dphi', 'trustemp'):
        plt.clf()
        plt.plot(X.get(k), X.dx, 'b.', label='dx')
        plt.plot(X.get(k), X.dy, 'r.', label='dy')
        plt.title(tt)
        plt.legend()
        plt.xlabel(k)
        plt.ylabel('Fiducial offset (um)')
        plt.savefig('fid-%s.png' % k)

    # We're going to re-read the metrology table, because all the positions have been updated
    # during dm.measure_spots.

    fn = dm.find_file('metrology')
    from astropy.table import Table
    dm.metro = Table.read(fn)

    mcopy = dm.metro.copy()
    M = Twrapper(mcopy)
    #M.orig_x_fp = M.x_fp.copy()
    #M.orig_y_fp = M.y_fp.copy()

    del mcopy['DEVICE_ID']
    del mcopy['BUS_ID']
    del mcopy['CAN_ID']
    mcopy.write('fp-metro-pre.csv', overwrite=True)
    
    applied_dx = []
    applied_dy = []

    devids = np.unique(fids.devid)
    for d in devids:
        I = np.flatnonzero(fids.devid == d)
        # drop largest and smallest dx,dy (cheap 'sigma-clipping')
        dx = 1000. * fids.dev_dx[I]
        dy = 1000. * fids.dev_dy[I]
        print(len(I), 'measurements for dev', d, 'with mean dx,dy  %+5.1f, %+5.1f and std %4.1f, %4.1f' %
              (np.mean(dx), np.mean(dy), np.std(dx), np.std(dy)))
        Kx = np.argsort(dx)
        Ky = np.argsort(dy)
        dx = dx[Kx[1:-1]]
        dy = dy[Ky[1:-1]]
        print('After dropping min & max:', len(dx), 'meas, mean dx,dy %+5.1f, %+5.1f and std %4.1f, %4.1f' %
              (np.mean(dx), np.mean(dy), np.std(dx), np.std(dy)))
        dx = np.mean(dx) / 1000.
        dy = np.mean(dy) / 1000.

        IM = np.flatnonzero(M.location == d)
        print(len(IM), 'metrology entries for location', d, M.device_type[IM[0]])
        #print('Metro x_fp,y_fp', M.x_fp[IM], M.y_fp[IM])
        #print('Fid x_fp, y_fp:', fids.dev_x[I], fids.dev_y[I])

        applied_dx.append(dx)
        applied_dy.append(dy)

        M.x_fp[IM] += dx
        M.y_fp[IM] += dy
        M.provenance[IM] = 'fiducial-systematics.py'

    mcopy.write('fp-metro-post.csv', overwrite=True)

    # I attempted to convert the shift + rotation of the GIFs into shift + rotation of the GFAs
    # -- not complete!!
    # petals = devids // 1000
    # applied_dx = np.array(applied_dx)
    # applied_dy = np.array(applied_dy)
    # for p in np.unique(petals):
    #     I = np.flatnonzero(petals == p)
    #     print('Petal', p, ': average dx,dy %.1f, %.1f um' %
    #           (1000. * np.mean(applied_dx[I]), 1000. * np.mean(applied_dy[I])))
    #
    #     ig1 = np.flatnonzero(devids == 1000*p + 541)
    #     ig2 = np.flatnonzero(devids == 1000*p + 542)
    #     if not(len(ig1) == 1 and len(ig2) == 1):
    #         print('Petal', p, ': found', len(ig1), 'and', len(ig2), 'offsets for GIFs')
    #         continue
    #     dx1 = applied_dx[ig1]
    #     dy1 = applied_dy[ig1]
    #     dx2 = applied_dx[ig2]
    #     dy2 = applied_dy[ig2]
    #     G1 = np.flatnonzero((M.device_loc == 541) * (M.petal_loc == p))
    #     G2 = np.flatnonzero((M.device_loc == 542) * (M.petal_loc == p))
    #     if not(len(G1) == 4 and len(G2) == 4):
    #         print('Petal', p, ': Found', len(G1), 'and', len(G2), 'GIFs -- not fixing GFAs')
    #         continue
    #     GIF1 = M[G1]
    #     GIF2 = M[G2]
    #
    #     I = np.flatnonzero((M.petal_loc == p) * (M.device_type == 'GFA'))
    #     if len(I) != 4:
    #         print('Petal', p, ': Found', len(I), 'GFA points -- not fixing GFAs')
    #         continue
    #     GFA = M[I]
    #
    #     # Vectors from GIF1 to GIF2, before and after moving them.
    #     gdx2 = np.mean(GIF2.x_fp) - np.mean(GIF1.x_fp)
    #     gdy2 = np.mean(GIF2.y_fp) - np.mean(GIF1.y_fp)
    #     gdx = np.mean(GIF2.orig_x_fp) - np.mean(GIF1.orig_x_fp)
    #     gdy = np.mean(GIF2.orig_y_fp) - np.mean(GIF1.orig_y_fp)
    #
    #     print('gdx2,gdy2', gdx2,gdy2)
    #     print('gdx ,gdy ', gdx ,gdy )
    #
    #     rot = (gdx * gdy2 - gdx2 * gdy) / (np.hypot(gdx, gdy) * np.hypot(gdx2, gdy2))
    #     print('rot:', rot)
    #     th = np.arcsin(rot)
    #     print('angle:', np.rad2deg(th))
    #
    #     # check:
    #     R = np.array([[np.cos(th), np.sin(th)], [-np.sin(th), np.cos(th)]])
    #     rx,ry = np.dot(R, np.array([gdx, gdy]))
    #     print('Rotated gdx,gdy:', rx,ry)
    #     print('Rotated other way:', np.dot(R.T, np.array([gdx, gdy])))
    #
    #     avdx = (dx1 + dx2) / 2.
    #     avdy = (dy1 + dy2) / 2.
    #     print('Average GIF dx,dy: %.1f, %.1f um' % (1000.*avdx, 1000.*avdy))
    #
    #     cx = (np.mean(GIF1.orig_x_fp) + np.mean(GIF2.orig_x_fp)) / 2.
    #     cy = (np.mean(GIF1.orig_y_fp) + np.mean(GIF2.orig_y_fp)) / 2.
    #     print('max dist from GIF midpoint:', np.max(np.hypot(GFA.x_fp - cx, GFA.y_fp - cy)))
    #
    #     def transform(x,y):
    #         dxy = np.array([x - cx, y - cy])
    #         print('tx: dxy', dxy)
    #         R = np.array([[np.cos(th), np.sin(th)], [-np.sin(th), np.cos(th)]])
    #         print('R:', R)
    #         rxy = np.dot(R, dxy)
    #         print('tx: rxy', rxy)
    #         newx = cx + rxy[0] + avdx
    #         newy = cy + rxy[1] + avdy
    #         return newx,newy
    #
    #     tx,ty = transform(GIF1.orig_x_fp, GIF1.orig_y_fp)
    #     print('Transformed GIF1: dists', np.hypot(GIF1.x_fp - tx, GIF1.y_fp - ty))
    #     tx,ty = transform(GIF2.orig_x_fp, GIF2.orig_y_fp)
    #     print('Transformed GIF2: dists', np.hypot(GIF2.x_fp - tx, GIF2.y_fp - ty))
    # M.delete_column('orig_x_fp')
    # M.delete_column('orig_y_fp')

    outdir = 'dm-fid-sys'
    dm.write_desimeter(outdir)

    import sys
    sys.exit(0)

    # Check!
    newdm = Desimeter(desimeter_dir=outdir, proc_data_dir='proc-fid-sys')
    fids = []
    #for expnum in [52644] + list(range(55353, 55357+1)):
    for keep,expnum in [(False, 52644)] + [(True, e) for e in list(range(55353, 55692+1))]:
        fvcfn = newdm.find_file('fvc', expnum=expnum)
        if fvcfn is None:
            continue
        print('Expnum', expnum)
        print('FVC', fvcfn)
        frame = 0
        fn = newdm.find_file('fvc-spots', expnum=expnum, frame=frame)
        if os.path.exists(fn):
            spots = astropy.table.Table.read(fn)
        else:
            spots = newdm.measure_spots(expnum, frame)
            spots.write(fn, overwrite=True)
        # At first I thought we could get away without Z-B, or a constant Z-B, but that very
        # much does not work due to, I assume, how different offsets & scales interact with Z-B.
        # Zero out the Z-B corrections (zbfit=False does not touch them)
        #newdm.fvc2fp.zbcoeffs = newdm.fvc2fp.zbpolids = None
        #spots = newdm.refit_spots(spots, zbfit=False)
        spots = newdm.refit_spots(spots)
        D = plot_fiducial_offsets(spots, expnum=expnum, frame=frame)
        fn = 'fvc-fid-sys-%08i-F%04i.png' % (expnum, frame)
        plt.savefig(fn)
        print('Wrote', fn)
        plt.clf()
        if not keep:
            continue
        fids.append(D._table)
    fids = astropy.table.vstack(fids)
    fids = Twrapper(fids)
    print('After:', len(fids), 'fiducial measurements,', len(np.unique(fids.devid)), 'unique fids')
    d = np.hypot(fids.dev_dx, fids.dev_dy)
    I = np.flatnonzero(d < 0.1)
    fids = fids[I]
    print('Cut to', len(fids), 'based on residuals')
    plt.clf()
    plt.hist(d * 1000., range=(0, 60), bins=30)
    rms = np.sqrt(np.mean((d*1000.)**2))
    plt.xlabel('Offset of fiducials (um)')
    plt.title('After: Fiducials, 2020-03-14 & 15: RMS2D %.1f um' % rms)
    plt.savefig('fid-offset-hist-after.png')

if __name__ == '__main__':
    main()
