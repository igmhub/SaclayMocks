import numpy as np
import fitsio
import matplotlib.pyplot as plt
import glob
from SaclayMocks import util


# indir="/global/cscratch1/sd/tetourne/DesiMocks/delta_g/mock_0/output/eboss-raw"
# indir="/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2"
indir="/global/project/projectdirs/eboss/lya_forest/dr16/redo_4_zbins/Fits/cf/Kaiser_sky_met_hcd/z_0_10"
# filenames="cf_z_0_10-exp_L0*.h5"
# ref_file = 'cf_z_0_10-exp_Rogers_free_L0.h5'
filenames = 'result_L0*.h5'
ref_file = 'result_L0_free_betaLYA_free.h5'
title = 'DR16 data'
par = 'chi2'

files = np.sort(glob.glob(indir+"/"+filenames))
f1, ax1 = plt.subplots()

for f in files:
    if 'free' in f: continue
    if 'prior' in f: continue
    print("Reading {}".format(f))
    pars = util.extract_h5file(f)
    idx1 = f.find('L0')
    idx2 = f.rfind('.')
    L_hcd = float(f[idx1+2:idx2])
    ax1.plot(L_hcd, pars[2][par], marker='o', color='blue')

pars = util.extract_h5file(indir+"/"+ref_file)
L_hcd = pars[2]['L0_hcd']
L_hcd_err = pars[3]['L0_hcd']
ax1.errorbar(L_hcd, pars[2][par], xerr=L_hcd_err, marker='o', color='red')
ax1.set_ylabel(r'$\chi^2_{red}$')
ax1.set_xlabel(r'$L_{\mathrm{HCD}}$')
ax1.grid()
f1.tight_layout()

plt.show()
