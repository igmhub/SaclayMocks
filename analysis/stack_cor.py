import numpy as np
import fitsio
import argparse
import time


parser = argparse.ArgumentParser()
parser.add_argument("--in-files", type=str, nargs='*')
parser.add_argument("--out", type=str)
# parser.add_argument("--nrp", type=int, default=50)
# parser.add_argument("--nrt", type=int, default=50)
args = parser.parse_args()

t0 = time.time()
nfiles = len(args.in_files)
nrp = fitsio.read_header(args.in_files[0], ext='COR')['NP']
nrt = fitsio.read_header(args.in_files[0], ext='COR')['NT']
z = np.zeros((nfiles, nrt*nrp))
nb = np.zeros((nfiles, nrt*nrp))
da = np.zeros((nfiles, nrt*nrp))
err = np.zeros((nfiles, nrt*nrp))
co = np.zeros((nfiles, nrt*nrp, nrt*nrp))
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
    co_tmp = fits[1].read()['CO']
    if 0 in co_tmp: print("Warning! There are bins where CO is null")
    co[i] = co_tmp
    fits.close()

# we = np.zeros_like(err)
# msk = err != 0
# we[msk] = 1 / err[msk]
# we_sum = np.sum(we, axis=0)
# z = np.sum(z*we, axis=0) / we_sum
# nb = np.sum(nb*we, axis=0) / we_sum
# da = np.sum(da*we, axis=0) / we_sum
# cov = 1 / np.sum(1/co,axis=0)

# z = np.mean(z, axis=0)
# nb = np.mean(nb, axis=0)
# da = np.mean(da, axis=0)
# cov = np.mean(co, axis=0) / nfiles

we = np.zeros_like(err)
msk = err != 0
we[msk] = 1 / err[msk]
we_sum = np.sum(we, axis=0)
z = np.sum(z*we, axis=0) / we_sum
nb = np.sum(nb*we, axis=0) / we_sum
da = np.sum(da*we, axis=0) / we_sum
# cov = np.mean(co, axis=0) / nfiles
cov = np.sum(co*np.sum(we,axis=1).reshape(-1,1,1), axis=0) / np.sum(we)
cov /= nfiles


print("Writting output in {}".format(args.out))
outfits = fitsio.FITS(args.out, 'rw', clobber=True)
table = [rp, rt, z, da, cov, dm, nb]
names = ['RP', 'RT', 'Z', 'DA', 'CO', 'DM', 'NB']
outfits.write(table, names=names, extname='COR', header=head)
outfits.close()
print("Done. {} s".format(time.time() - t0))
