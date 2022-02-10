import fitsio
import argparse
import numpy as np

"""
This codes convert the DLA catalog produced by the SaclayMocks to a 
DLA catalog readable by picca (to compute DLA auto-correlation, or
to remove DLA from delta files for instance
"""

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", type=str, default=None, required=True, help='path to input dla catalog file')
parser.add_argument("-o", type=str, default=None, required=True, help='path to output dla catalog file')
# parser.add_argument("--master", type=str, default=None, required=True, help='path to master file')
parser.add_argument("--nhi-min", type=float, default=None, required=False, help='remove DLA with log(n_HI) < NHI_MIN')
parser.add_argument("--nhi-max", type=float, default=None, required=False, help='remove DLA with log(n_HI) > NHI_MIN')
args = parser.parse_args()

### Data
print("Reading {}".format(args.i))
data = fitsio.read(args.i, ext='DLACAT')
if args.nhi_min is not None:
    msk1 = data['N_HI_DLA'] > args.nhi_min
    print("Removing {} DLAs with log(n_HI) < {} out of {} total DLAs".format((~msk1).sum(), args.nhi_min, data.size))
else:
    msk1 = np.bool_(np.ones(data.size))

if args.nhi_max is not None:
    msk2 = data['N_HI_DLA'] < args.nhi_max
    print("Removing {} DLAs with log(n_HI) > {} out of {} total DLAs".format((~msk2).sum(), args.nhi_max, data.size))
else:
    msk2 = np.bool_(np.ones(data.size))

msk = msk1 & msk2
outfits = fitsio.FITS(args.o, 'rw', clobber=True)
names = ['RA', 'DEC', 'THING_ID', 'Z', 'PLATE', 'MJD', 'FIBERID', 'NHI', 'ZQSO']
table = []
table.append(data['RA'][msk])
table.append(data['DEC'][msk])
table.append(data['MOCKID'][msk])
table.append(data['Z_DLA_RSD'][msk])
table.append(data['MOCKID'][msk])
table.append(data['MOCKID'][msk])
table.append(data['MOCKID'][msk])
table.append(data['N_HI_DLA'][msk])
table.append(data['Z_QSO_RSD'][msk])

print("Writting out {}".format(args.o))
outfits.write(table, names=names, extname='DLACAT')
outfits.close()
