import numpy as np
import fitsio
import matplotlib.pyplot as plt
import glob
from SaclayMocks import util


# indir="/global/cscratch1/sd/tetourne/DesiMocks/delta_g/mock_0/output/eboss-raw"
indir="/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.0"
filenames="cf_z_0_10-exp_rmin*.h5"
cf_file = indir+"/"+"cf_z_0_10-exp.fits"
fits = fitsio.FITS(cf_file)
r = np.sqrt(fits[1].read()['RP']**2 + fits[1].read()['RT']**2)
fits.close()
rmax = 180
free_pars = 4
par1 = 'chi2'
par2 = 'beff_LYA'
par3 = 'beta_LYA'

files = np.sort(glob.glob(indir+"/"+filenames))
f1, ax1 = plt.subplots()
f2, ax2 = plt.subplots()
f3, ax3 = plt.subplots()

for f in files:
    if 'smoothing' in f: continue
    print("Reading {}".format(f))
    pars = util.extract_h5file(f)
    idx1 = f.find('rmin')
    rmin = int(f[idx1+4:idx1+6])
    dof = ((r>rmin)&(r<rmax)).sum()
    print("dof is {}".format(dof))
    ax1.plot(rmin, pars[2][par1]/(dof-free_pars), marker='o', color='blue')
    ax2.errorbar(rmin, pars[2][par2], yerr=pars[3][par2], marker='o', color='blue')
    ax3.errorbar(rmin, pars[2][par3], yerr=pars[3][par3], marker='o', color='blue')

ax1.set_ylabel(r'$\chi^2_{red}$')
ax2.set_ylabel(r'$b_{\mathrm{eff},\mathrm{Ly}\alpha}$')
ax3.set_ylabel(r'$\beta_{\mathrm{Ly}\alpha}$')
for ax in [ax1, ax2, ax3]:
    ax.set_xlabel(r'$r_{\mathrm{min}}$')
    ax.grid()
for f in [f1, f2, f3]:
    f.tight_layout()

plt.show()
