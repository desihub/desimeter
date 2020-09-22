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
    dx = xfp_meas - xfp_metro
    dy = yfp_meas - yfp_metro
    dr = np.sqrt(dx**2+dy**2)
    print("median offset = %4.1f um" % (np.median(dr[dr!=0])*1000.));
    ii=np.where((xfp_metro!=0)&(dr<3.))[0]
    plt.figure("quiver",figsize=(10,10))
    plt.quiver(xfp_meas[ii],yfp_meas[ii],dx[ii],dy[ii])
    plt.plot(xfp_meas[ii][pinhole[ii]>0],yfp_meas[ii][pinhole[ii]>0],".",c="red")
    plt.title('Before refitting Z-B')
    plt.savefig('quiver1.png')

    zbpolids = np.array([0,1,2,3,4,5,6,9,20,27,28,29,30])
    x1b, y1b = dm.refit_zb(xfp_meas, yfp_meas, xfp_metro, yfp_metro, zbpolids)

    x1,y1 = xfp_meas,  yfp_meas
    x2,y2 = xfp_metro, yfp_metro
    nmad2d=1.20*np.median(np.sqrt((x1-x2)**2+(y1-y2)**2))
    rms2d=np.sqrt(np.mean((x1-x2)**2+(y1-y2)**2))
    print('Before refitting Z-B: MAD2D %.1f um, RMS2D %.1f um' %
          (nmad2d*1000., rms2d*1000.))

    dxb = x1b - xfp_metro
    dyb = y1b - yfp_metro
    plt.clf()
    plt.quiver(x1b[ii], y1b[ii], dxb[ii], dyb[ii])
    plt.plot(x1b[ii][pinhole[ii]>0], y1b[ii][pinhole[ii]>0], ".", c="red")
    plt.title('After refitting Z-B')
    plt.savefig('quiver2.png')

    nmad2d_b = 1.20*np.median(np.sqrt((x1b-x2)**2+(y1b-y2)**2))
    rms2d_b = np.sqrt(np.mean((x1b-x2)**2+(y1b-y2)**2))
    print('After  refitting Z-B: MAD2D %.1f um, RMS2D %.1f um' %
          (nmad2d_b*1000., rms2d_b*1000.))

if __name__ == '__main__':
    main()
