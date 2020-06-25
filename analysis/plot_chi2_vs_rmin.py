import numpy as np
import fitsio
import matplotlib.pyplot as plt
import glob
from SaclayMocks import util


# indir="/global/cscratch1/sd/tetourne/DesiMocks/delta_g/mock_0/output/eboss-raw"
indir="/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.0"
filenames="cf_z_0_10-exp_rmin*.h5"
para = "beta_LYA"


files = np.sort(glob.glob(indir+"/"+filenames))
f1, ax1 = plt.subplots()
f2, ax2 = plt.subplots()

for f in files:
    if 'smoothing' in f: continue
    pars = util.extract_h5file(f)
    idx1 = f.find('rmin')
    rmin = int(f[idx1+4:idx1+6])
    val = pars[2][para]
    if para == "chi2" :
        dof = util.dof(rmin, 180) - 4
        print(rmin, val, dof)
        ax1.plot(rmin, val, color='blue', marker='o')
        ax1.plot(rmin, dof, color='green', marker='o')
        ax2.plot(rmin, val/dof, color='blue', marker='o')
    else :
        err_val = pars[3][para]
        print(rmin, val)
        ax1.errorbar(rmin, val, yerr=err_val, marker='+', color='blue')

ax1.plot([], [], color='blue', label=para)
if para == "chi2":
    ax1.plot([], [], color='green', label='d.o.f.')
    ax2.plot([], [], color='blue', label='chi2/d.o.f.')
ax1.set_xlabel(r'$r_{\mathrm{min}}$')
ax2.set_xlabel(r'$r_{\mathrm{min}}$')
ax1.grid()
ax2.grid()
ax1.legend()
ax2.legend()

plt.show()
