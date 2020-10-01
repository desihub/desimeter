import os
import numpy as np
import pylab as plt
import astropy.table
from desimeter.desimeter import Desimeter

def main():
    dm = Desimeter(data_dir='arccal', proc_data_dir='proc')

    expnum = 52644

    allspots = []
    for frame in range(6):
        fn = dm.find_file('fvc-spots', expnum=expnum, frame=frame)
        if os.path.exists(fn):
            spots = astropy.table.Table.read(fn)
        else:
            spots = dm.measure_spots(expnum, frame)
            spots.write(fn, overwrite=True)

        plot_fiducial_offsets(spots)
        plt.savefig('fvc-fid-offsets-%08i-F%04i.png' % (expnum, frame))

        allspots.append(spots)

    cfn = dm.find_file('fvc-circles', expnum=expnum, tag='F0_5')
    if os.path.exists(cfn):
        circles = astropy.table.Table.read(cfn)
    else:
        circles = dm.measure_circles(allspots)
        circles.write(cfn, overwrite=True)

    xfp_meas  = circles['X_FP']
    yfp_meas  = circles['Y_FP']
    xfp_metro = circles['X_FP_METRO']
    yfp_metro = circles['Y_FP_METRO']
    pinhole   = circles['PINHOLE_ID']

    def report_result(x, y, xmetro,ymetro, pinhole,
                      title, plotfn):
        dx = x - xmetro
        dy = y - ymetro
        nmad2d=1.20*np.median(np.sqrt((dx)**2+(dx)**2))
        rms2d=np.sqrt(np.mean((dx)**2+(dy)**2))
        dr = np.sqrt(dx**2+dy**2)
        ii=np.where((xmetro!=0)&(dr<3.))[0]
        plt.figure("quiver",figsize=(10,10))
        plt.clf()

        # smaller "scale" makes the arrows longer.
        qargs = dict(pivot='middle', angles='xy', scale_units='xy',
                     scale=0.005)
        plt.quiver(x[ii],y[ii],dx[ii],dy[ii], **qargs)

        i_fid = ii[pinhole[ii] > 0]
        plt.quiver(x[i_fid],y[i_fid],dx[i_fid],dy[i_fid], color='r', **qargs)
        plt.plot(x[i_fid],y[i_fid],".",c="red")
        txt = ('%s: MAD %.1f um, RMS %.1f um' %
               (title, nmad2d*1000., rms2d*1000.))
        plt.title(txt)
        plt.savefig(plotfn)
        print(txt)

    x1,y1 = xfp_meas,  yfp_meas
    x2,y2 = xfp_metro, yfp_metro

    report_result(x1,y1, x2,y2, pinhole, 'Before refitting Z-B',
                  'quiver1.png')

    zbpolids = np.array([0,1,2,3,4,5,6,9,20,27,28,29,30])
    x1b, y1b = dm.refit_zb(xfp_meas, yfp_meas, xfp_metro, yfp_metro, zbpolids)

    report_result(x1b,y1b, x2,y2, pinhole, 'After  refitting Z-B',
                  'quiver2.png')

    def apply_rigid(sx, sy, rot, x1, y1):
        c,s = np.cos(np.deg2rad(rot)), np.sin(np.deg2rad(rot))
        R = np.array([[c, s], [-s, c]])
        rxy = np.dot(R, np.vstack((x1,y1)))
        xx = rxy[0,:] + sx
        yy = rxy[1,:] + sy
        return xx, yy

    def rms_rigid(X, x1,y1, x2,y2):
        sx, sy, rot = X
        xx,yy = apply_rigid(sx, sy, rot, x1,y1)
        return np.sqrt(np.mean((xx - x2)**2 + (yy - y2)**2))

    # Try moving each petal independently...
    x1c = np.zeros_like(x1b)
    y1c = np.zeros_like(y1b)
    from scipy.optimize import minimize
    loc = circles['LOCATION']
    for petal in range(10):
        I = np.flatnonzero(loc // 1000 == petal)
        R = minimize(rms_rigid, np.zeros(3),
                     args=(x1b[I],y1b[I], xfp_metro[I],yfp_metro[I]))
        sx,sy,rot = R.x
        xr,yr = apply_rigid(sx, sy, rot, x1b[I], y1b[I])
        x1c[I] = xr
        y1c[I] = yr

    report_result(x1c,y1c, x2,y2, pinhole, 'After shift/rot',
                  'quiver3.png')

    x1d, y1d = dm.refit_zb(x1c,y1c, xfp_metro, yfp_metro, zbpolids)

    report_result(x1d,y1d, x2,y2, pinhole, 'After second Z-B',
                  'quiver4.png')

    dm.write_desimeter('arccal-dm')

def plot_zb_terms():
    from desimeter.transform import zhaoburge as zb
    xx,yy = np.meshgrid(np.linspace(-1, 1, 20),
                        np.linspace(-1, 1, 20))
    rr = xx**2 + yy**2
    I = np.flatnonzero(rr < 1)
    xx = xx.flat[I]
    yy = yy.flat[I]
    for i in range(zb.NCOEFS):
        zbx,zby,label = zb.getZhaoBurgeTerm(i, xx, yy)
        plt.clf()
        plt.quiver(xx, yy, zbx, zby, pivot='middle')
        plt.axis('equal')
        plt.title('Zhao-Burge term %i: %s' % (i, label))
        plt.savefig('zb-%02i.png' % i)

if __name__ == '__main__':
    main()
