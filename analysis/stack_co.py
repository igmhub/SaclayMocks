import numpy as np
import fitsio
import argparse
import time


parser = argparse.ArgumentParser()
parser.add_argument("--in-files", type=str, nargs='*')
parser.add_argument("--out", type=str)
parser.add_argument("--nrp", type=int, default=50)
parser.add_argument("--nrt", type=int, default=50)
args = parser.parse_args()

t0 = time.time()
nfiles = len(args.in_files)
z = np.zeros((nfiles, args.nrt*args.nrp))
nb = np.zeros((nfiles, args.nrt*args.nrp))
da = np.zeros((nfiles, args.nrt*args.nrp))
err = np.zeros((nfiles, args.nrt*args.nrp))
co = np.zeros((nfiles, args.nrt*args.nrp, args.nrt*args.nrp))
fits = fitsio.FITS(args.in_files[0])
rp = fits[1].read()['RP']
rt = fits[1].read()['RT']
dm = fits[1].read()['DM']
head = fits[1].read_header()
fits.close()

for i,f in enumerate(args.in_files):
    print("Reading {}".format(f))
    fits = fitsio.FITS(f)
    z[i] = fits[1].read()['Z']
    nb[i] = fits[1].read()['NB']
    da[i] = fits[1].read()['DA']
    err[i] = np.diag(fits[1].read()['CO'])
    co[i] = fits[1].read()['CO']
    fits.close()

we = 1 / err
we_sum = np.sum(we, axis=0)
z = np.sum(z*we, axis=0) / we_sum
nb = np.sum(nb*we, axis=0) / we_sum
da = np.sum(da*we, axis=0) / we_sum
co = np.mean(co, axis=0) / np.sqrt(nfiles)

print("Writting output in {}".format(args.out))
outfits = fitsio.FITS(args.out, 'rw', clobber=True)
table = [rp, rt, z, da, co, dm, nb]
names = ['RP', 'RT', 'Z', 'DA', 'CO', 'DM', 'NB']
outfits.write(table, names=names, extname='COR', header=head)
outfits.close()
print("Done. {} s".format(time.time() - t0))
