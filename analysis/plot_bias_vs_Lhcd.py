import glob
import fitsio
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from SaclayMocks import util
from iminuit import Minuit


SMALL_SIZE = 15
MEDIUM_SIZE = 18
BIGGER_SIZE = 20
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize

plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.rc('figure', figsize=(16,7))


# files = np.sort(glob.glob("/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2/cf_z_0_10-exp_L*.h5"))
files = np.sort(glob.glob("/global/project/projectdirs/eboss/lya_forest/dr16/redo_4_zbins/Fits/cf/Kaiser_sky_met_hcd/z_0_10/result_L*.h5"))
# files = files[[0,2,4,5,6,7,8,9]]
files = files[:-2]
print(files)

f1, ax1 = plt.subplots(ncols=2)
f2, ax2 = plt.subplots(ncols=2)
# Plot
for f in files:
    pars = util.extract_h5file(f)
    if pars[2]['L0_hcd'] == 10:
        color = 'tab:red'
    else:
        color = 'tab:blue'
    ax1[0].errorbar(pars[2]['L0_hcd'], np.abs(pars[2]['beff_LYA']), yerr=pars[3]['beff_LYA'], fmt='o', color=color)
    ax1[1].errorbar(pars[2]['L0_hcd'], np.abs(pars[2]['beta_LYA']), yerr=pars[3]['beta_LYA'], fmt='o', color=color)
    ax2[0].errorbar(pars[2]['L0_hcd'], np.abs(pars[2]['bias_hcd']), yerr=pars[3]['bias_hcd'], fmt='o', color=color)
    ax2[1].errorbar(pars[2]['L0_hcd'], np.abs(pars[2]['beta_hcd']), yerr=pars[3]['beta_hcd'], fmt='o', color=color)

ax1[0].set_xlabel(r'$L_{\mathrm{HCD}}$')
ax1[1].set_xlabel(r'$L_{\mathrm{HCD}}$')
ax2[0].set_xlabel(r'$L_{\mathrm{HCD}}$')
ax2[1].set_xlabel(r'$L_{\mathrm{HCD}}$')

ax1[0].grid()
ax1[1].grid()
ax2[0].grid()
ax2[1].grid()

ax1[0].set_ylabel(r'$|b_{\mathrm{eff}, \mathrm{Ly}\alpha}|$')
ax1[1].set_ylabel(r'$\beta_{\mathrm{Ly}\alpha}$')
ax2[0].set_ylabel(r'$|b_{\mathrm{HCD}}|$')
ax2[1].set_ylabel(r'$\beta_{\mathrm{HCD}}$')

f1.tight_layout()
f2.tight_layout()

plt.show()
