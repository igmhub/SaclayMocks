import numpy as np
import fitsio
import matplotlib.pyplot as plt


# indir = "/global/projecta/projectdirs/desi/mocks/lya_forest/develop/london/v9.0/v9.0.42"
indir = "/global/cscratch1/sd/tetourne/DesiMocks/v4.7/mock_1/output/"

master = fitsio.FITS(indir+"/master.fits")[1].read()
dla = fitsio.FITS(indir+"/master_DLA.fits")[1].read()
zmin = 1.8
cut_nhi = 20.3
print(indir)
print("n_dla / n_qso (z>{}) : {} / {} = {}".format(zmin, dla.size, master[master['Z_QSO_RSD']>zmin].size, np.float32(dla.size)/master[master['Z_QSO_RSD']>zmin].size))
print("n_dla / n_qso (z>{}, NHI>{}) : {} / {} = {}".format(zmin, cut_nhi, dla[dla['N_HI_DLA']>cut_nhi].size, master[master['Z_QSO_RSD']>zmin].size, np.float32(dla[dla['N_HI_DLA']>cut_nhi].size)/master[master['Z_QSO_RSD']>zmin].size))

# N_HI distrib
f, ax = plt.subplots()
ax.hist(dla['N_HI_DLA'], bins=100, density=True)
ax.set_xlabel('log(N_HI)')
ax.set_yscale('log')

# dn/dz
f, ax = plt.subplots()
ax.hist(dla['Z_DLA_RSD'], bins=100, histtype='step', label='17.2 < log(NHI) < 22.5')
ax.hist(dla['Z_DLA_RSD'][dla['N_HI_DLA']>20.3], bins=100, histtype='step', label='20.3 < log(NHI) < 22.5')
ax.hist(dla['Z_DLA_RSD'][dla['N_HI_DLA']<20.3], bins=100, histtype='step', label='17.2 < log(NHI) < 20.3')
ax.legend()
ax.set_xlabel('z')

# dndz ratio with z of DLA
redshift_bins = np.linspace(1.5, 4, 101)
redshift = (redshift_bins[1:] + redshift_bins[:-1]) / 2
dla_hist = np.histogram(dla['Z_DLA_RSD'][dla['N_HI_DLA']>cut_nhi], bins=redshift_bins)[0]
qso_hist = np.histogram(master['Z_QSO_RSD'], bins=redshift_bins)[0]
y = np.zeros_like(redshift)
msk = qso_hist != 0
y[msk] = np.float32(dla_hist[msk]) / qso_hist[msk]
f, ax = plt.subplots()
ax.plot(redshift, y)
ax.set_xlabel('z')
ax.set_ylabel('nz_dla / nz_qso')
ax.set_title('DLA cut_nhi at log(N_HI) > {}'.format(cut_nhi))
ax.grid()

# dndz ratio with z of host QSO
redshift_bins = np.linspace(1.5, 4, 101)
redshift = (redshift_bins[1:] + redshift_bins[:-1]) / 2
dla_hist = np.histogram(dla['Z_QSO_RSD'][dla['N_HI_DLA']>cut_nhi], bins=redshift_bins)[0]
qso_hist = np.histogram(master['Z_QSO_RSD'], bins=redshift_bins)[0]
y = np.zeros_like(redshift)
msk = qso_hist != 0
y[msk] = np.float32(dla_hist[msk]) / qso_hist[msk]
f, ax = plt.subplots()
ax.plot(redshift, y)
ax.set_xlabel('z')
ax.set_ylabel('nz_qso_host / nz_qso')
ax.set_title('DLA cut_nhi at log(N_HI) > {}'.format(cut_nhi))
ax.grid()

plt.show()
