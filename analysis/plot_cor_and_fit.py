import os, sys
import numpy as np
import matplotlib.pyplot as plt
import h5py
import fitsio
import picca.wedgize

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

# indir_cor = "/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2"
indir_cor = "/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-raw"
# indir_cor = "/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_10/eboss-raw"
# indir_cor = "/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/v4.7.30/eboss-raw"
# indir_cor = "/global/project/projectdirs/eboss/lya_forest/dr16/redo_4_zbins/Correlations/"

indir_fit = indir_cor
# indir_fit = "/global/project/projectdirs/eboss/lya_forest/dr16/redo_4_zbins/Fits/cf/Kaiser_sky_met_hcd/"

patern_cf = "xcf_z_{zmin}_{zmax}-exp.fits"
# patern_cf = "co_qso_z_{zmin}_{zmax}-exp.fits"

patern_fit = patern_cf.replace('.fits', '.h5')
# patern_fit = "z_{zmin}_{zmax}/result.h5"


ylabel = r"$r^{2}\xi(r) \, [(h^{-1} \mathrm{Mpc})^2]$"

mu_bins = np.array([0, 0.5, 0.8, 0.95, 1])
# mu_bins = np.array([0,1])
# mu_bins = np.array([0, 0.2, 0.5, 1])
colors = ['b', 'g', 'orange', 'red']
# colors = ['b', 'g', 'red']

z_bins = np.array(['0', '2.35', '2.65', '3.05', '10'])

pred = False
patern_pred = patern_cf.replace('-exp', '-exp_pred')
r_pow = 2

nrp = fitsio.read_header(indir_cor+"/"+patern_cf.format(zmin=z_bins[0],zmax=z_bins[1]), ext='COR')['NP']
nrt = fitsio.read_header(indir_cor+"/"+patern_cf.format(zmin=z_bins[0],zmax=z_bins[1]), ext='COR')['NT']
rtmin = 0
rtmax = 200
rpmax = 200
if nrp == 100:
    mu_label = r"${} < |\mu| < {}$"
    rpmin = -200
else:
    mu_label = r"${} < \mu < {}$"
    rpmin = 0

f, axs = plt.subplots(nrows=2, ncols=2)

for i in range(len(z_bins)-1):
    filename_fits = indir_cor+"/"+patern_cf.format(zmin=z_bins[i], zmax=z_bins[i+1])
    print("Reading {}".format(filename_fits))
    data = fitsio.FITS(filename_fits)
    da = data[1].read()['DA']
    co = data[1].read()['CO']
    data.close()
    if pred:
        filename_fits = indir_cor+"/"+patern_pred.format(zmin=z_bins[i], zmax=z_bins[i+1])
        print("Reading {}".format(filename_fits))
        data = fitsio.FITS(filename_fits)
        da_pred = data[1].read()['DA']
        co_pred = data[1].read()['CO']
        data.close()
    filename_h5 = indir_fit+"/"+patern_fit.format(zmin=z_bins[i], zmax=z_bins[i+1])
    print("Reading {}".format(filename_h5))
    ff = h5py.File(filename_h5)
    try :
        fit = ff["LYA(LYA)xLYA(LYA)/fit"][...]
    except:
        print("Can't read LYA(LYA)xLYA(LYA)/fit")
        try:
            fit = ff["LYA(LYA)xQSO/fit"][...]
        except:
            print("Can't read LYA(LYA)xQSO/fit")
            try:
                fit = ff["QSOxQSO/fit"][...]
            except:
                print("Can't read QSOxQSO/fit")
                try:
                    fit = ff["HCDxHCD/fit"][...]
                except:
                    print("Can't read HCDxHCD/fit")
                    try:
                        fit = ff["cf_z_0_10/fit"][...]
                    except:
                        print("Can't read cf_z_0_10/fit")
                        try:
                            idx1 = int(filename_h5.find(".h5"))
                            idx2 = int(filename_h5.rfind("/"))+1
                            fit = ff[filename_h5[idx2:idx1]+"/fit"][...]
                        except:
                            print("Can't read "+filename_h5[idx2:idx1]+"/fit")
                            try :
                                idx1 = int(filename_fits.find(".fits"))
                                idx2 = int(filename_fits.rfind("/"))+1
                                fit = ff[filename_fits[idx2:idx1]+"/fit"][...]
                            except:
                                print("Can't read "+filename_fits[idx2:idx1]+"/fit")
                                sys.exit()
    ff.close()
    for j in range(len(mu_bins)-1):
        w = picca.wedgize.wedge(mumin=mu_bins[j],mumax=mu_bins[j+1], rtmax=rtmax, rpmax=rpmax, rtmin=rtmin, rpmin=rpmin, nrt=nrt, nrp=nrp,absoluteMu=True)
        wedge_data = w.wedge(da,co)
        coef = wedge_data[0]**r_pow
        wedge_fit = w.wedge(fit,co)
        axs[i//2, i%2].errorbar(wedge_data[0],coef*wedge_data[1],yerr=coef*np.sqrt(np.diag(wedge_data[2])),fmt='+', label=mu_label.format(mu_bins[j], mu_bins[j+1]), color=colors[j])
        axs[i//2, i%2].plot(wedge_fit[0],coef*wedge_fit[1], color=colors[j])
        if pred:
            wedge_pred = w.wedge(da_pred, co_pred)
            axs[i//2, i%2].plot(wedge_pred[0],coef*wedge_pred[1], linestyle='--', color=colors[j])
    axs[i//2, i%2].grid()
    if i == 0:
        axs[i//2, i%2].legend()
    if i%2 == 0:
        axs[i//2, i%2].set_ylabel(ylabel)
    if i//2 == 1:
        axs[i//2, i%2].set_xlabel(r"$r \, [h^{-1} \mathrm{Mpc}]$")
    # axs[i//2, i%2].text(0.85, 0.95, r'{} < z < {}'.format(z_bins[i], z_bins[i+1]))
    axs[i//2, i%2].set_title(r'{} < z < {}'.format(z_bins[i], z_bins[i+1]))
    plt.tight_layout()
plt.show()
