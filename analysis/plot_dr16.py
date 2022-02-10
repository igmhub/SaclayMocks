import matplotlib.pyplot as plt
import numpy as np
import h5py
import fitsio
from SaclayMocks import util


cor_dir = "/global/cscratch1/sd/tetourne/Analysis/redo_dr16/Correlations"
fit_dir = "/global/cscratch1/sd/tetourne/Analysis/redo_dr16/Fits/cf/Kaiser_sky_met_hcd"

### Read CF files
data1 = "/cf_z_0_2.35-exp.fits"
data2 = "/cf_z_2.35_2.65-exp.fits"
data3 = "/cf_z_2.65_3.05-exp.fits"
data4 = "/cf_z_3.05_10-exp.fits"

### Read fit files
fit1 = "/z_0_2.35/result.h5"
fit2 = "/z_2.35_2.65/result.h5"
fit3 = "/z_2.65_3.05/result.h5"
fit4 = "/z_3.05_10/result.h5"


data = [data1, data2, data3, data4]
fit = [fit1, fit2, fit3, fit4]
zbins = [r'$0<z<2.35$', r'$2.35<z<2.65$', r'$2.65<z<3.05$', r'$3.05<z<10$']
mubins = [0, 0.5, 0.8, 0.95, 1]
colors = ['mediumblue', 'green', 'darkorange', 'red']

### Plot 
for i in range(len(data)):
    f = fitsio.FITS(cor_dir+data[i])
    da = f[1].read()['DA']
    co = f[1].read()['CO']
    rp = f[1].read()['RP']
    rt = f[1].read()['RT']
    f.close()
    rr = np.sqrt(rp**2 + rt**2)
    f = h5py.File(fit_dir+fit[i])
    idx1 = fit[i].rfind("/")+1
    idx2 = fit[i].rfind("-")
    # fit_name = fit[i][idx1:idx2]
    fit_name = 'cf_z_0_10'
    fi = f[fit_name]['fit'][...]
    f.close()
    plt.figure()
    for j in range(len(mubins)-1):
        util.add_wedge(da, co, errorbar=True, mumin=mubins[j], mumax=mubins[j+1], color=colors[j], fmt='+')
        util.add_wedge(fi, co, errorbar=False, mumin=mubins[j], mumax=mubins[j+1], color=colors[j], label=r'${} < \mu < {}$'.format(mubins[j], mubins[j+1]))
    plt.grid()
    plt.legend()
    plt.xlabel('r [Mpc/h]')
    plt.ylabel(r'$r^2 \xi(r)$')
    plt.title(zbins[i])
    plt.tight_layout()

plt.show()
