import numpy as np
import numpy.ma as ma
import fitsio
from SaclayMocks import powerspectrum, util, transmissions, constant
import argparse
import time


parser = argparse.ArgumentParser()
parser.add_argument("--in-dir", type=str, help="v4.7/v4.7.01/")
parser.add_argument("--out-dir", type=str, help="v4.7/v4.7.01/P1d/")
parser.add_argument("--n-files", type=float, default=None)
args = parser.parse_args()

k_bins = np.linspace(0, 2, 40)
z_bins = [2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4]
dz = 0.1
pixel_size = 0.2
# wav_max = constant.lya
wav_max = 1200
wav_min = constant.lylimit
# wav_min = 1040
print("dz is {}".format(dz))
print("considering pixels in {} < lambda < {}".format(wav_min, wav_max))

print("Reading transmission files ...")
t0 = time.time()
trans = transmissions.ReadTransmission(args.in_dir, nfiles=args.n_files, read_dla=False)
print("Done. {} s.".format(time.time() - t0))
z_trans = trans.wavelength / constant.lya - 1
wav_rf = trans.wavelength / (trans.metadata['Z'].reshape(-1,1) + 1)  # (nspec, npix)
msk_wav = ((wav_rf < wav_max) & (wav_rf > wav_min))

names = ['k', 'p1d', 'p1derr']
len_forest = []

for zi in z_bins:
    print("Computing P1D for {} < z < {} ...".format(zi-dz, zi+dz))
    t0 = time.time()
    msk_1d = (z_trans >= zi - dz) & (z_trans < zi + dz)
    msk_2d = np.bool_(msk_1d * np.ones_like(trans.transmission))
    msk = msk_wav & msk_2d
    mean_trans = ma.mean(trans.transmission[msk])
    mean_z = ma.mean((z_trans*np.ones_like(trans.transmission))[msk])
    print("<F> = {} at z_eff = {}".format(mean_trans, mean_z))
    p1d = powerspectrum.ComputeP1D(pixel_size)
    ll = 0
    cc = 0
    for i in range(len(trans.transmission)):
        w_rf = trans.wavelength / (trans.metadata['Z'][i] + 1)
        msk = (w_rf < wav_max) & (w_rf > wav_min) & msk_1d
        tt = trans.transmission[i][msk]
        if len(tt[tt.mask==False]) < 1: continue
        else:
            flux = tt[tt.mask==False]
        delta = flux / mean_trans - 1
        p1d.add_spectrum(delta)
        ll += len(delta)
        cc += 1
    k, pk, pkerr = p1d.P1D(k_bins)
    table = [k, pk, pkerr]
    outname = args.out_dir+"/p1d_z_{}_{}.fits".format(np.round(zi-dz,3), np.round(zi+dz,3))
    outfits = fitsio.FITS(outname, 'rw', clobber=True)
    outfits.write(table, names=names, extname='P1D')
    outfits[1].write_key('zmin', zi-dz)
    outfits[1].write_key('zmax', zi+dz)
    outfits.close()
    print("P1D file {} written. {} s".format(outname, time.time()-t0))
    len_forest.append(ll/cc)
print(len_forest)
np.save(args.out_dir+"/len_forest.npy", len_forest)
