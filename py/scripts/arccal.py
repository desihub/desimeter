import os
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

if __name__ == '__main__':
    main()
