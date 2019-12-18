#!/usr/bin/env python

# first version of code copied from https://desi.lbl.gov/trac/browser/code/focalplane/plate_control/trunk/pecs/poscalibrationfits.py#L48

import argparse
import sys
import numpy as np
import psycopg2
import yaml
import datetime

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""FVC image processing""")

parser.add_argument('-o','--outfile', type = str, default = None, required = True, help = 'output yaml file')
parser.add_argument('--host', type = str, default = "db.replicator.dev-cattle.stable.spin.nersc.org", required = False, help = 'db host, use desi-db at KPNO')
parser.add_argument('--port', type = str, default = "60042", required = False, help = 'db port, use 5442 at KPNO')
parser.add_argument('--database', type = str, default = "desi_dev", required = False, help = 'db port')


args  = parser.parse_args()

           
conn = psycopg2.connect(host=args.host, port=args.port, database=args.database, user="desi_reader", password="reader")

# find the latest snapshot version, may want and older one
cx = conn.cursor()
cx.execute('SELECT version from snapshots')
res = cx.fetchall()
snap_ver = np.max(res)

# read petal alignment for all petals
query="SELECT constants.elements.name, constants.elements.constants FROM constants.groups, constants.elements,constants.group_to_snapshot, constants.snapshots,constants.element_to_group WHERE constants.group_to_snapshot.snapshot_name = 'DESI' AND constants.group_to_snapshot.snapshot_version = '{snap_ver}' AND constants.groups.name = 'focal_plane_metrology' AND constants.groups.version = constants.group_to_snapshot.group_version AND constants.snapshots.version = constants.group_to_snapshot.snapshot_version AND constants.group_to_snapshot.group_name = constants.element_to_group.group_name AND constants.group_to_snapshot.group_version = constants.element_to_group.group_version AND constants.element_to_group.element_name = constants.elements.name AND constants.element_to_group.element_version = constants.elements.version ORDER BY constants.elements.name".format(snap_ver=snap_ver)
cx.execute(query)
res = cx.fetchall()
alignments = {int(petal_loc): {'Tx': data['petal_offset_x'],
                               'Ty': data['petal_offset_y'],
                               'Tz': data['petal_offset_z'],
                               'alpha': data['petal_rot_1'],
                               'beta': data['petal_rot_2'],
                               'gamma': data['petal_rot_3']} for ( petal_loc, data ) in res }# zip(res.index, res['constants'])}

if args.outfile.find(".yaml")>0 :
    ofile = open(args.outfile, 'w')
    ofile.write("# generate by {} on {}, snapshot version = {}\n".format(sys.argv[0],datetime.date.today().strftime("%Y-%m-%d"),snap_ver))
    yaml.dump(alignments, ofile)
    ofile.close()
else :
    print("outfile {} does not have a yaml extension so either change the outfile extension as .yaml or write another output style in code".format(args.outfile))
    sys.exit(12)
print("wrote {}".format(args.outfile))


      
