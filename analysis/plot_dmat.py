import fitsio
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from SaclayMocks import util
from iminuit import Minuit
import h5py


# SMALL_SIZE = 11
# MEDIUM_SIZE = 14
# BIGGER_SIZE = 14
# plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
# plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
# plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
# plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
# plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
# plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize

# plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
# plt.rc('figure', figsize=(11,8))

cf_raw = fitsio.FITS("/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-raw/cf_z_0_10-exp.fits")
cf_00 = fitsio.FITS("/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.0/cf_z_0_10-exp.fits")
fit_raw = h5py.File("/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-raw/cf_z_0_10-exp.h5")
fit_00 = h5py.File("/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.0/cf_z_0_10-exp.h5")

da_raw = cf_raw[1].read()['DA']
co_raw = cf_raw[1].read()['CO']
dmat_raw = cf_raw[1].read()['DM']
da_00 = cf_00[1].read()['DA']
co_00 = cf_00[1].read()['CO']
dmat_00 = cf_00[1].read()['DM']
cf_fit_raw = fit_raw['cf_z_0_10-exp/fit'][...]
cf_fit_00 = fit_00['LYA(LYA)xLYA(LYA)/fit'][...]

mu_bins = np.array([0, 0.5, 0.8, 0.95, 1])
for i in range(len(mu_bins)-1):
    plt.figure()
    util.add_wedge(da_raw, co_raw, mumin=mu_bins[i], mumax=mu_bins[i+1], color='blue', fmt='+', label='raw mocks')
    util.add_wedge(cf_fit_raw, co_raw, mumin=mu_bins[i], mumax=mu_bins[i+1], color='blue', linestyle='--')
    util.add_wedge(da_00, co_00, mumin=mu_bins[i], mumax=mu_bins[i+1], color='red', fmt='+', label='mock eboss-0.0')
    util.add_wedge(cf_fit_00, co_00, mumin=mu_bins[i], mumax=mu_bins[i+1], color='red', linestyle='--')
    # util.add_wedge(np.dot(dmat_00,da_raw), co_raw, mumin=mu_bins[i], mumax=mu_bins[i+1], color='green', marker='.', errorbar=False)
    plt.xlabel(r'$r\;[h^{-1}\mathrm{Mpc}]$')
    plt.ylabel(r'$r^2\xi(r)\;[(h^{-1}\mathrm{Mpc})^2]$')
    plt.grid()
    plt.legend()
    plt.title(r'${} < \mu < {}$'.format(mu_bins[i], mu_bins[i+1]))
    plt.tight_layout()

plt.show()
