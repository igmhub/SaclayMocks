import fitsio
import argparse
import time
from LyaMocks import util
import scipy as sp


parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str)
parser.add_argument("-o", type=str)
parser.add_argument("--to-do", type=str, help='QSO or DLA', default='QSO')
parser.add_argument("--downsampling-z-cut", default=0, type=float)
parser.add_argument("--downsampling-z-cut-max", default=10, type=float)
parser.add_argument("--downsampling-nb", default=None, type=float)
parser.add_argument("--rsd", type=str, default='True')
parser.add_argument("--rand", type=str, default='False')
args = parser.parse_args()

print("Starging conversion of master.fits")
print("Input is {}".format(args.i))
print("Output is {}".format(args.o))
t0 = time.time()
rsd_cond = util.str2bool(args.rsd)
rand_cond = util.str2bool(args.rand)
downsampling_z_cut = args.downsampling_z_cut
downsampling_z_cut_max = args.downsampling_z_cut_max
downsampling_nb = args.downsampling_nb
master = fitsio.FITS(args.i)
ra = master[1]['RA'][:].astype('float64')
dec = master[1]['DEC'][:].astype('float64')
if not rand_cond:
    if rsd_cond:
        zqso = master[1]['Z_{}_RSD'.format(args.to_do)][:]
    else:
        zqso = master[1]['Z_{}_NO_RSD'.format(args.to_do)][:]
else:
    zqso = master[1]['Z'][:]

thid = master[1]['MOCKID'][:]
master.close()

print("Applying downsampling...")
print("zmin={}".format(downsampling_z_cut))
print("zmax={}".format(downsampling_z_cut_max))
if zqso[(zqso>downsampling_z_cut)&(zqso<downsampling_z_cut_max)].size<downsampling_nb:
    print('WARNING:: Trying to downsample, when nb cat = {} and nb downsampling = {}'.format(zqso[(zqso>downsampling_z_cut)&(zqso<downsampling_z_cut_max)].size,downsampling_nb) )
    select = (zqso>downsampling_z_cut)&(zqso<downsampling_z_cut_max)
else:
    select_fraction = downsampling_nb/((zqso>downsampling_z_cut)&(zqso<downsampling_z_cut_max)).sum()
    select = sp.random.choice(sp.arange(ra.size),size=int(ra.size*select_fraction),replace=False)

ra = ra[select]
dec = dec[select]
zqso = zqso[select]
thid = thid[select]

### Save
out = fitsio.FITS(args.o,'rw',clobber=True)
cols = [ra,dec,thid,thid,thid,thid,zqso]
names = ['RA','DEC','THING_ID','PLATE','MJD','FIBERID','Z']
out.write(cols, names=names)
out.close()
print("Took {} s".format(time.time() -t0))
