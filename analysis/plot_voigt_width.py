import scipy as sp
import numpy as np
import glob
import fitsio
from SaclayMocks import util
import matplotlib.pyplot as plt


# pyplot config
SMALL_SIZE = 16
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

lamvoi,tau,Fvoi,logn,nindex=np.loadtxt('../etc/voigttable.out', unpack=True)
nindex=nindex.astype('int')
ln=np.zeros(60)
for i in range(60):
    ln[i]=logn[nindex==i][0]
dl1=np.zeros(60)
dl2=np.zeros(60)
for i in range(60):
    m1=(nindex==i)&(Fvoi<0.8)
    dl1[i]=max(lamvoi[m1])-min(lamvoi[m1])
    m2=(nindex==i)&(Fvoi<0.5)
    dl2[i]=max(lamvoi[m2])-min(lamvoi[m2])
mmm=(ln>=17.3)&(ln<=22.4)
#mmm=(ln>=20.)&(ln<=21.)
plt.plot(ln[mmm],dl1[mmm],lw=3,color='red',label='Fvoigt<0.8') 
plt.plot(ln[mmm],dl2[mmm],lw=3,color='blue',label='Fvoigt<0.1')
plt.plot([20.3,20.3],[0.,90.],color='black')
plt.text(19.51,70., 'Not Masked',fontsize=16)
plt.text(20.35,70., 'Masked',fontsize=16)
plt.legend(loc='upper left',fontsize=14)
plt.grid()
plt.xlabel('logNHI',fontsize=14)
plt.ylabel(r'$\Delta\lambda$ (Angstrom)  $\approx r_\parallel$(Mpc)',fontsize=14)

for i in [17.6, 18.2, 18.8, 19.4, 20]:
    idx = np.argsort(np.abs(ln-i))[0]
    print("log n_HI = {} give delta_lambda = {} (Fvoigt < 0.8)".format(i, dl1[idx]))
    print("log n_HI = {} give delta_lambda = {} (Fvoigt < 0.5)".format(i, dl2[idx]))

# Read DLAs
dla = fitsio.FITS("/global/project/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output_old/master_DLA.fits")
val, nhi_bins = np.histogram(dla[1].read()['N_HI_DLA'], bins=70)
bin_center = (nhi_bins[1:] + nhi_bins[:-1]) / 2
f_nhi = sp.interpolate.interp1d(bin_center, val)
weights = f_nhi(ln[mmm])
mean_dl1 = np.average(dl1[mmm], weights=weights)
mean_dl2 = np.average(dl2[mmm], weights=weights)

print("delta_lambda mean is {} (Fvoigt < 0.8)".format(mean_dl1))
print("delta_lambda mean is {} (Fvoigt < 0.5)".format(mean_dl2))

plt.show()
