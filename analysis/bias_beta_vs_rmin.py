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
# indir = "/global/cscratch1/sd/tetourne/Analysis/redo_dr16/Fits/cf/Bias_beta_vs_rmin/z_0_10"
# indir = "/global/cscratch1/sd/tetourne/Out/mock_0.2/Fits"
indir = "/global/cscratch1/sd/tetourne/Out/v4.7.22/from_transmissions/Fit_pred"
# title = 'dr16'
title = 'prediction'

files = np.sort(glob.glob(indir+"/*0.h5"))
print(files)
rmin_list = []
bias_eta_list = []
bias_eta_err_list = []
beta_list = []
beta_err_list = []
bias_list = []
bias_err_list = []
beff_list = []
beff_err_list = []

for f in files:
    pars = util.extract_h5file(f)
    i = f.rfind('rmin')
    rmin_list.append(np.int32(f[i+4:i+6]))
    bias_eta = pars[2]['bias_eta_LYA']
    beta = pars[2]['beta_LYA']
    f = pars[2]['growth_rate']
    bias_eta_err = pars[3]['bias_eta_LYA']
    beta_err = pars[3]['beta_LYA']
    cov = pars[4]['cov[beta_LYA, bias_eta_LYA]']
    bias = bias_eta * f / beta
    # bias_err = np.sqrt((f*bias_eta_err/beta)**2 + (bias_eta*f*beta_err/beta**2)**2 - 2*cov/bias_eta/beta)
    bias_err = np.sqrt((f*bias_eta_err/beta)**2 + (bias_eta*f*beta_err/beta**2)**2 - 2*cov*bias_eta*f**2/beta**3)
    beff = bias_eta*f * np.sqrt(1+2/3*beta+1/5*beta**2) / beta
    # db_dbiaseta = f*np.sqrt(1+2/3*beta+1/5*beta**2)/beta
    # db_dbeta = (bias_eta*f/beta)*np.sqrt(1+2/3*beta+1/5*beta**2)*((1/3+1/5*beta)/(1+2/3*beta+1/5*beta**2) - 1/beta)
    # beff_err = np.sqrt((db_dbiaseta*bias_eta_err)**2 + (db_dbeta*beta_err)**2 + 2*db_dbiaseta*db_dbeta*cov)
    beff_err = util.beff_err(bias_eta, bias_eta_err, beta, beta_err, cov, f)
    bias_eta_list.append(bias_eta)
    bias_eta_err_list.append(bias_eta_err)
    beta_list.append(beta)
    beta_err_list.append(beta_err)
    bias_list.append(bias)
    bias_err_list.append(bias_err)
    beff_list.append(beff)
    beff_err_list.append(beff_err)

print("\nbias_eta:")
print(bias_eta_list)
print(bias_eta_err_list)
print("\beta:")
print(beta_list)
print(beta_err_list)
print("\nbias:")
print(bias_list)
print(bias_err_list)
print("\beff:")
print(beff_list)
print(beff_err_list)

f1, ax1 = plt.subplots()
ax1.errorbar(rmin_list, bias_eta_list, yerr=bias_eta_err_list, fmt='o')
ax1.set_xlabel(r'$r_{min}$ [Mpc/h]')
ax1.set_ylabel(r'$b_{\eta}$')
ax1.grid()
ax1.set_title(title)
plt.tight_layout()

f2, ax2 = plt.subplots()
ax2.errorbar(rmin_list, beta_list, yerr=beta_err_list, fmt='o')
ax2.set_xlabel(r'$r_{min}$ [Mpc/h]')
ax2.set_ylabel(r'$\beta$')
ax2.grid()
ax2.set_title(title)
plt.tight_layout()

f3, ax3 = plt.subplots()
ax3.errorbar(rmin_list, bias_list, yerr=bias_err_list, fmt='o')
ax3.set_xlabel(r'$r_{min}$ [Mpc/h]')
ax3.set_ylabel(r'$b$')
ax3.grid()
ax3.set_title(title)
plt.tight_layout()

f4, ax4 = plt.subplots()
ax4.errorbar(rmin_list, beff_list, yerr=beff_err_list, fmt='o')
ax4.set_xlabel(r'$r_{min}$ [Mpc/h]')
ax4.set_ylabel(r'$b_{eff}$')
ax4.grid()
ax4.set_title(title)
plt.tight_layout()

plt.show()
