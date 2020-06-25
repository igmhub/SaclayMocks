import numpy as np
import matplotlib.pyplot as plt
from SaclayMocks import util


indir = "/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_10/eboss-raw/"
# files = [indir+"/co_qso_z_0_2.35-exp.h5",
#          indir+"/co_qso_z_2.35_2.65-exp.h5",
#          indir+"/co_qso_z_2.65_3.05-exp.h5",
#          indir+"/co_qso_z_3.05_10-exp.h5"]
files = [indir+"/co_qso_z_0_2.35-exp_fixed_f.h5",
         indir+"/co_qso_z_2.35_2.65-exp_fixed_f.h5",
         indir+"/co_qso_z_2.65_3.05-exp_fixed_f.h5",
         indir+"/co_qso_z_3.05_10-exp_fixed_f.h5"]

fmt = '.'
label = 'raw mocks'
color = 'orangered'
color_model = 'royalblue'
use_legends = False
plot_model = True

redshift = np.zeros(len(files))
bias = np.zeros(len(files))
beta = np.zeros(len(files))
beta_err = np.zeros(len(files))
growth_rate = np.zeros(len(files))
growth_rate_err = np.zeros(len(files))
cov = np.zeros(len(files))

for i, f in enumerate(files):
    pars = util.extract_h5file(f)
    redshift[i] = pars[2]['zeff']
    beta[i] = pars[2]['beta_QSO']
    beta_err[i] = pars[3]['beta_QSO']
    growth_rate[i] = pars[2]['growth_rate']
    growth_rate_err[i] = pars[3]['growth_rate']
    if 'cov[beta_QSO, growth_rate]' in pars[4]:
        cov[i] = pars[4]['cov[beta_QSO, growth_rate]']

bias = growth_rate / beta
bias_err = util.bias_err(growth_rate, growth_rate_err, beta, beta_err, cov)

# Model
zz = np.linspace(1.8, 3.6, 10000)
bias_model = util.bias_qso(zz)
beta_model = util.growthRateStructure(zz) / bias_model

# Plot bias
fig1, ax1 = plt.subplots()
ax1.errorbar(redshift, bias, yerr=bias_err, fmt=fmt, label=label, color=color)
if plot_model:
    plt.plot(zz, bias_model, color=color_model)
if use_legends:
    ax1.legend()
ax1.grid()
ax1.set_xlabel('z')
ylabel = r'$b_{\mathrm{QSO}}$'
ax1.set_ylabel(ylabel)
plt.tight_layout()

# Plot beta
fig2, ax2 = plt.subplots()
ax2.errorbar(redshift, beta, yerr=beta_err, fmt=fmt, label=label, color=color)
if plot_model:
    plt.plot(zz, beta_model, color=color_model)
if use_legends:
    ax2.legend()
ax2.grid()
ax2.set_xlabel('z')
ylabel = r'$\beta_{\mathrm{QSO}}$'
ax2.set_ylabel(ylabel)
plt.tight_layout()

plt.show()
