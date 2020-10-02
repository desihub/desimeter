import os
from desimeter.desimeter import Desimeter

def main():
    dm = Desimeter(desimeter_dir='dm-fid-sys', proc_data_dir='proc-dither')

    expnum = 55670
    frame = 1

    fn = dm.find_file('fieldmodel', expnum=expnum, frame=frame)
    if not os.path.exists(fn):
        fm = dm.fit_guide_star_coords(expnum, frame)
        with open(fn, 'w') as f:
            f.write(fm.tojson())
    else:
        with open(fn) as f:
            fm = FieldModel.fromjson(f.read())


if __name__ == '__main__':
    main()
