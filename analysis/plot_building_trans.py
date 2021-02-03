import fitsio
import numpy as np
import matplotlib.pyplot as plt
from SaclayMocks import util, constant


SMALL_SIZE = 13
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


i = -5  # -5 have nice dla
zcat = fitsio.FITS("/global/project/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/v4.7.07/eboss-0.2/zcat.fits")
flux =  zcat[1].read()['FLUX_G'] + zcat[1].read()['FLUX_R'] + zcat[1].read()['FLUX_Z']
idx_zcat = np.argsort(flux)
hp_pix = util.radec2pix(ra=zcat[1].read()['RA'][idx_zcat[i]], dec=zcat[1].read()['DEC'][idx_zcat[i]], nside=16)
th_id = zcat[1].read()['TARGETID'][idx_zcat[i]]
z = zcat[1].read()[idx_zcat[i]]['Z']
spec = fitsio.FITS("/global/project/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/v4.7.07/eboss-0.2/spectra-16/{}/{}/spectra-16-{}.fits".format(hp_pix//100, hp_pix, hp_pix))
idx_spec = np.where(spec[1]['TARGETID'].read()==th_id)[0][0]

trans = fitsio.FITS("/global/project/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_07/chunk_{}/spectra_merged_with_delta_g/spectra_merged-{}-{}.fits.gz".format(str(th_id)[0], hp_pix, str(th_id)[1:4]))
idx_trans = np.where(trans[1].read()['THING_ID']==th_id)[0][0]
wav = trans['LAMBDA'].read()
c_of_z = util.InterpFitsTable('../etc/params.fits', 'z', 'c')
zz = wav / constant.lya - 1
cc = c_of_z.interp(zz)
delta_l = trans['DELTA_L'].read()[idx_trans]
delta_s = trans['DELTA_S'].read()[idx_trans]
eta_par = trans['ETA_PAR'].read()[idx_trans]
flux = trans['FLUX'].read()[idx_trans]


# Plot spectra
f0, ax0 = plt.subplots()
ax0.plot((spec['B_WAVELENGTH'].read())[:-600], spec['B_FLUX'].read()[idx_spec][:-600], color='b')
ax0.plot((spec['R_WAVELENGTH'].read())[600:], spec['R_FLUX'].read()[idx_spec][600:], color='r')
# ax0.axvline(constant.lya, linestyle='--', color='grey')
# ax0.axvline(1031.48, linestyle='--', color='grey')
ax0.grid()
ax0.set_xlabel(r'$\lambda_{\mathrm{obs}}\,[\AA{}]$')
ax0.set_ylabel(r"$flux \;[ 10^{-17} \; erg \cdot{} s^{-1} \cdot{} cm^{-2} \cdot{} \AA{}^{-1}]$")

# Plot transmission
f, axs = plt.subplots(nrows=4, sharex=True, figsize=(11,14))
# axs[0].set_xlim(constant.lyb-10, constant.lya+5)
axs[0].set_xlim(1120, 1200)
axs[-1].set_xlabel(r'$\lambda_{\mathrm{RF}}\,[\AA{}]$')

axs[0].plot(wav/(1+z), delta_l, alpha=0.5, color='b', label=r'$\delta_l$')
axs[0].plot(wav/(1+z), cc*eta_par, alpha=0.5, color='g', label=r'$c(z)\eta_{\parallel}$')
axs[0].plot(wav/(1+z), delta_l + cc*eta_par, alpha=0.5,color='r', label=r'$\delta_l + c(z)\eta_{\parallel}$')
axs[0].set_ylim(-4,4)
# axs[0].set_ylabel(r"$\delta_l$")

# axs[1].plot(wav/(1+z), delta_l+cc*eta_par)
# axs[1].set_ylim(-5,5)
# axs[1].set_ylabel(r"$\delta_l + c(z)\eta_{\parallel}$")
axs[1].plot(wav/(1+z), delta_l+delta_s, alpha=0.5, label=r'$\delta_l + \delta_s$')
axs[1].plot(wav/(1+z), delta_l+delta_s+cc*eta_par, alpha=0.5, label=r'$\delta_g$')
axs[1].set_ylim(-18,18)
# axs[1].set_ylabel(r"$\delta_l + c(z)\eta_{\parallel}$")


# axs[2].plot(wav/(1+z), delta_l+cc*eta_par + delta_s)
# axs[2].set_ylim(-18,18)
# axs[2].set_ylabel(r"$\delta_g$")

axs[2].plot(wav/(1+z), flux, label=r'$F$')
axs[2].set_ylim(-0.1,1.1)
# axs[2].set_ylabel(r"$F$")

axs[3].plot(spec['B_WAVELENGTH'].read()/(1+z), spec['B_FLUX'].read()[idx_spec], label=r'$flux$')
axs[3].set_ylim(-5,100)
# axs[3].set_ylabel(r"$flux$")

for ax in axs:
    ax.grid()
    ax.legend()


f.tight_layout()
f.savefig("building_trans.pdf")
print("{} saved.".format("building_trans.pdf"))
# plt.show()
