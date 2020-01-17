import glob
import h5py
import fitsio
import numpy as np
import matplotlib.pyplot as plt
import picca.wedgize
from SaclayMocks import util


'''
This code plots the evolution of bias and beta when we change rmin in the fitter
with all the other parameters fixed to the value found with the standard fit
'''

SMALL_SIZE = 18
MEDIUM_SIZE = 20
BIGGER_SIZE = 22
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize

plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.rc('figure', figsize=(9,7))

# Parameters
indir = "/global/cscratch1/sd/tetourne/Analysis/redo_dr16/Fits/cf/Bias_beta_vs_rmin/z_0_10"

files = glob.glob(indir+"/*.h5")
rmin = []
bias_eta = []
bias_eta_err = []
beta = []
beta_err = []
bias = []
beff = []
for f in files:
    pars = util.extract_h5file(f)
    i = f.rfind('rmin')
    rmin.append(np.int32(f[i+4
                           :-3]))
    bias_eta.append(pars[2]['bias_eta_LYA'])
    bias_eta_err.append(pars[3]['bias_eta_LYA'])
    beta.append(pars[2]['beta_LYA'])
    beta_err.append(pars[3]['beta_LYA'])
    bias.append(pars[2]['bias_eta_LYA']*pars[2]['growth_rate']/pars[2]['beta_LYA'])
    beff.append(pars[2]['bias_eta_LYA']*pars[2]['growth_rate']/pars[2]['beta_LYA']*np.sqrt(1+2/3*pars[2]['beta_LYA']+1/5*pars[2]['beta_LYA']**2))

f1, ax1 = plt.subplots()
ax1.errorbar(rmin, bias_eta, yerr=bias_eta_err, fmt='o')
ax1.set_xlabel(r'$r_{min}$ [Mpc/h]')
ax1.set_ylabel(r'$b_{\eta}$')
ax1.grid()
plt.tight_layout()

f2, ax2 = plt.subplots()
ax2.errorbar(rmin, beta, yerr=beta_err, fmt='o')
ax2.set_xlabel(r'$r_{min}$ [Mpc/h]')
ax2.set_ylabel(r'$\beta$')
ax2.grid()
plt.tight_layout()

f3, ax3 = plt.subplots()
ax3.plot(rmin, bias, 'o')
ax3.set_xlabel(r'$r_{min}$ [Mpc/h]')
ax3.set_ylabel(r'$b$')
ax3.grid()
plt.tight_layout()

f4, ax4 = plt.subplots()
ax4.plot(rmin, beff, 'o')
ax4.set_xlabel(r'$r_{min}$ [Mpc/h]')
ax4.set_ylabel(r'$b_{eff}$')
ax4.grid()
plt.tight_layout()

plt.show()
