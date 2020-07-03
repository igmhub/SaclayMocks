import os, sys
import numpy as np
import matplotlib.pyplot as plt
import h5py
import fitsio
import picca.wedgize
import glob


SMALL_SIZE = 11
MEDIUM_SIZE = 14
BIGGER_SIZE = 14
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize

plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.rc('figure', figsize=(11,8))

# Reference correlation
ref_file = "/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-raw/xcf_z_0_10-exp.fits"

# List of correlations
patern_file = "/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/v4.7.*/eboss-raw/xcf_z_0_10-exp.fits"
files = np.sort(glob.glob(patern_file))


mu_bins = np.array([0, 0.5, 0.8, 0.95, 1])

r_pow = 2
rtmin = 0
rtmax = 200
rpmin = -200
rpmax = 200
nrt = 50
nrp = 100

colors = ['b', 'g', 'orange', 'red']
# colors = ['b', 'g', 'red']

ylabel = r"$r^{2}\xi(r) \, [(h^{-1} \mathrm{Mpc})^2]$"
if nrp == 100:
    mu_label = r"${} < |\mu| < {}$"
else:
    mu_label = r"${} < \mu < {}$"

# Read and Plot
f, axs = plt.subplots(nrows=2, ncols=2)
for i in range(len(mu_bins)-1):
    print("\n Bin: "+mu_label.format(mu_bins[i], mu_bins[i+1]))
    print("Reading {}".format(ref_file))
    data = fitsio.FITS(ref_file)
    da = data[1].read()['DA']
    co = data[1].read()['CO']
    data.close()
    w = picca.wedgize.wedge(mumin=mu_bins[i],mumax=mu_bins[i+1], rtmax=rtmax, rpmax=rpmax, rtmin=rtmin, rpmin=rpmin, nrt=nrt, nrp=nrp,absoluteMu=True)
    wedge_data = w.wedge(da,co)
    coef = wedge_data[0]**r_pow
    axs[i//2, i%2].plot(wedge_data[0],coef*wedge_data[1], color=colors[i]) #, label=mu_label.format(mu_bins[i], mu_bins[i+1]))

    for f in files:
        if "old" in f: continue
        # if 'v4.7.10' in f: continue
        print("Reading {}".format(f))
        data = fitsio.FITS(f)
        da = data[1].read()['DA']
        co = data[1].read()['CO']
        data.close()
        wedge_data = w.wedge(da,co)
        coef = wedge_data[0]**r_pow
        axs[i//2, i%2].plot(wedge_data[0],coef*wedge_data[1], color=colors[i], alpha=0.15)

    axs[i//2, i%2].grid()
    # axs[i//2, i%2].legend()
    if i%2 == 0:
        axs[i//2, i%2].set_ylabel(ylabel)
    if i//2 == 1:
        axs[i//2, i%2].set_xlabel(r"$r \, [h^{-1} \mathrm{Mpc}]$")
    axs[i//2, i%2].set_title(mu_label.format(mu_bins[i], mu_bins[i+1]))
    plt.tight_layout()
plt.show()

