import fitsio
import numpy as np
import matplotlib.pyplot as plt
from SaclayMocks import util


cflo00 = fitsio.FITS("/global/cfs/cdirs/desi/science/lya/picca_on_mocks/london/v9.0/global/eboss-0.0/stack/cf_z_0_10-exp.fits")
xcflo00 = fitsio.FITS("/global/cfs/cdirs/desi/science/lya/picca_on_mocks/london/v9.0/global/eboss-0.0/stack/xcf_z_0_10-exp.fits")
cflo02 = fitsio.FITS("/global/cfs/cdirs/desi/science/lya/picca_on_mocks/london/v9.0/global/eboss-0.2/stack/cf_z_0_10-exp.fits")
xcflo02 = fitsio.FITS("/global/cfs/cdirs/desi/science/lya/picca_on_mocks/london/v9.0/global/eboss-0.2/stack/xcf_z_0_10-exp.fits")

cfsa00 = fitsio.FITS("/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-11_20/eboss-0.0/cf_z_0_10-exp.fits")
xcfsa00 = fitsio.FITS("/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-11_20/eboss-0.0/xcf_z_0_10-exp.fits")
cfsa02 = fitsio.FITS("/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-11_20/eboss-0.2/cf_z_0_10-exp.fits")
xcfsa02 = fitsio.FITS("/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-11_20/eboss-0.2/xcf_z_0_10-exp.fits")

cfda = fitsio.FITS("/global/cscratch1/sd/tetourne/Analysis/dr16_paper/Correlations/cf_z_0_10-exp.fits.gz")
xcfda = fitsio.FITS("/global/cscratch1/sd/tetourne/Analysis/dr16_paper/Correlations/xcf_z_0_10-exp.fits.gz")

f1 = xcflo02
f2 = xcfsa02
nrp = 100
nrt = 50

da1 = util.convert1DTo2D(f1[1].read()['DA'], nbX=nrp, nbY=nrt)
co1 = util.convert1DTo2D(np.diag(f1[1].read()['CO']), nbX=nrp, nbY=nrt)
da2 = util.convert1DTo2D(f2[1].read()['DA'], nbX=nrp, nbY=nrt)
co2 = util.convert1DTo2D(np.diag(f2[1].read()['CO']), nbX=nrp, nbY=nrt)

# dadata = util.convert1DTo2D(cfda[1].read()['DA'], nbX=nrp, nbY=nrt)
# codata = util.convert1DTo2D(np.diag(cfda[1].read()['CO']), nbX=nrp, nbY=nrt)
dadata = util.convert1DTo2D(xcfda[1].read()['DA'], nbX=nrp, nbY=nrt)
codata = util.convert1DTo2D(np.diag(xcfda[1].read()['CO']), nbX=nrp, nbY=nrt)

rp = util.convert1DTo2D(f2[1].read()['RP'], nbX=nrp, nbY=nrt)
rt = util.convert1DTo2D(f2[1].read()['RT'], nbX=nrp, nbY=nrt)
rr = np.sqrt(rp**2 + rt**2)

# y = (da1 - da2) / np.sqrt(co2)
# y = (da1 - dadata) / np.sqrt(codata)
# y = (da1 - da2)*rr**2
y = (da1 - dadata)*rr**2

# extent = (0, 200, 0, 200)
extent = (0, 200, -200, 200)

# cbar = plt.imshow(y, aspect='auto', origin='lower', cmap='bwr', extent=extent, vmin=-3, vmax=3)
cbar = plt.imshow(y, origin='lower', cmap='bwr', extent=extent, vmin=-6, vmax=6)
plt.colorbar(cbar)

plt.xlabel(r'$r_{\perp}$ [Mpc/h]')
plt.ylabel(r'$r_{\parallel}$ [Mpc/h]')
# plt.title(r'$(XCF_{London} - XCF_{Saclay}) / \sigma_{Saclay}$ - eboss-0.2')
plt.title(r'$(XCF_{London} - XCF_{DR16})\times r^2$')
# plt.title(r'$(XCF_{\mathrm{London-0.2}} - XCF_{\mathrm{DR16}}) / \sigma_{\mathrm{DR16}}$')

plt.tight_layout()
plt.show()
