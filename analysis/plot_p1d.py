import numpy as np
import matplotlib.pyplot as plt
import fitsio
from SaclayMocks import util, constant
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--in-dir", type=str, help="v4.7/v4.7.01/P1d")
args = parser.parse_args()

# indir="/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/v4.7.01/eboss-raw/P1d"
# indir="/global/cscratch1/sd/tetourne/DesiMocks/debug/check_p1d/z_3.0/mock_0/output/P1d"
indir = args.in_dir

z_bins = [2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4]
dz = 0.1

filenames = []
for zi in z_bins:
    filename = indir+"/p1d_z_{}_{}.fits".format(np.round(zi-dz,1), np.round(zi+dz,1))
    filenames.append(filename)

### Plot
f, ax = plt.subplots()
for ff, zi, cc in zip(filenames, z_bins, constant.colors):
    # # Data
    # try:
    #     k_data, pk_data, pkerr_data = util.read_P1D(zi, mpc=True)
    #     ax.errorbar(k_data, pk_data, yerr=pkerr_data, fmt='+', color='green')
    # except: print("No data for this bin")
    # # Model
    k_mod, pk_mod = util.read_P1D_model(zi, mpc=True)
    ax.plot(k_mod, pk_mod, color=cc, label='z = {}'.format(np.round(zi,1)))
    # # Mocks
    fits = fitsio.FITS(ff)
    k_mock = fits[1].read()['k']
    pk_mock = fits[1].read()['p1d']
    pkerr_mock = fits[1].read()['p1derr']
    fits.close()
    # x, ya, yb = util.array_interp(k_mod, pk_mod, k_mock, pk_mock)
    # ax.plot(x, (yb - ya)/ya, color=cc, label='z = {}'.format(np.round(zi,1)))
    # ax.plot(x, (yb)/ya, color=cc, label='z = {}'.format(np.round(zi,1)))
    ax.errorbar(k_mock, pk_mock, yerr=pkerr_mock, fmt='o', color=cc)
    # ax.errorbar(k_mock, pk_mock, yerr=pkerr_mock, label='z = {}'.format(np.round(zi,1)))

# plt.xlabel(r'$k [h^{-1}\mathrm{Mpc}]$')
# plt.ylabel(r'$P^{\mathrm{1D}} [h \mathrm{Mpc}^{-1}]$')
ax.set_xlabel(r'$k [h^{-1}\mathrm{Mpc}]$')
ax.set_ylabel(r'$P^{\mathrm{1D}} [h \mathrm{Mpc}^{-1}]$')
ax.grid()
ax.legend()
plt.tight_layout()
plt.show()
