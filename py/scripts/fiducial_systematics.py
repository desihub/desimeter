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

    for p in np.unique(D.petal_loc):
        I = np.flatnonzero(D.petal_loc == p)
        plt.plot(D.dev_x[I], D.dev_y[I], 'o', mec='none', ms=25, alpha=0.1)
    #for x,y,d in zip(D.dev_x, D.dev_y, D.devid):
    #    plt.text(x, y, '%i' % d)
    qargs = dict(pivot='middle', angles='xy', scale_units='xy',
                scale=0.0005)
    plt.quiver(D.dev_x, D.dev_y, D.dev_dx, D.dev_dy, **qargs)
    # Draw lines around the petals.
    outline_petals(color='k', alpha=0.1)
    plt.axis('equal')
    rms2d = 1000. * np.sqrt(np.mean(D.dev_dx**2 + D.dev_dy**2))
    plt.title('Offsets of fiducials vs whole-focal-plane fit: %.1f um RMS2d, expnum %i frame %i' % (rms2d, expnum, frame))
    return D

if __name__ == '__main__':
    dm = Desimeter(proc_data_dir='proc')

    fids = []
    plt.figure(figsize=(10,10))
    for expnum in range(55353, 55692+1):
    #for expnum in range(55353, 55357+1):
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

        D = plot_fiducial_offsets(spots, expnum=expnum, frame=frame)
        # Format table of per-fiducial results
        D.expnum = np.zeros(len(D), np.int32) + expnum
        D.frame  = np.zeros(len(D), np.int16) + frame
        hdr = fitsio.read_header(fvcfn)
        for key in ['MJD-OBS', 'TARGTRA', 'TARGTDEC', 'TARGTAZ', 'TARGTEZ', 'AIRMASS', 'ADC1PHI', 'ADC2PHI']:
            D.set(key.replace('-','_').lower(), np.zeros(len(D), np.float32) + hdr.get(key, -99))
        fids.append(D._table)
        
        fn = 'fvc-fid-offsets-%08i-F%04i.png' % (expnum, frame)
        plt.savefig(fn)
        print('Wrote', fn)
        plt.clf()

    fids = astropy.table.vstack(fids)
    fids.write('fiducials.fits', overwrite=True)
    
