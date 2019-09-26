import fitsio
import argparse
import numpy as np


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--dla-cat", type=str, default=None, required=True, help='path to dla catalog file')
parser.add_argument("--master", type=str, default=None, required=True, help='path to master file')
parser.add_argument("--out", type=str, default=None, required=True, help='output drq')
args = parser.parse_args()

### Data
print("Reading {}".format(args.dla_cat))
h = fitsio.FITS(args.dla_cat)
master = fitsio.read(args.master, ext=1)
data = {}
data['Z'] = h[1]['Z_DLA'][:]
data['THING_ID'] = h[1]['MOCKID'][:]
data['PLATE'] = h[1]['MOCKID'][:]
data['MJD'] = h[1]['MOCKID'][:]
data['FIBERID'] = h[1]['MOCKID'][:]
qso_id = master['THING_ID']
data['RA'] = []
data['DEC'] = []
for i in data['THING_ID']:
    msk = (i == qso_id)
    data['RA'].append(master['RA'][msk])
    data['DEC'].append(master['DEC'][msk])
data['RA'] = np.array(data['RA'])
data['DEC'] = np.array(data['DEC'])
h.close()

### Save data
print(data['RA'].size)
print("Writting out {}".format(args.out))
out = fitsio.FITS(args.out,'rw',clobber=True)
cols = [ v for k,v in data.items() ]
names = [ k for k in data.keys() ]
out.write(cols,names=names)
out.close()
