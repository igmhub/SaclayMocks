import glob
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from SaclayMocks import util, powerspectrum, constant
from iminuit import Minuit


SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 16
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize

plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.rc('figure', figsize=(17,7))

# files = np.sort(glob.glob("/global/homes/t/tetourne/Git/picca/py/picca/fitter2/models/fvoigt_models/Fvoigt_v4.7.22_nhi_1*.txt"))
# files = files[1:]
files = []
files = np.append(files, "/global/homes/t/tetourne/Git/picca/py/picca/fitter2/models/fvoigt_models/Fvoigt_v4.7.22_highcut_nhi20.3.txt")
files = np.append(files, "/global/homes/t/tetourne/Git/picca/py/picca/fitter2/models/fvoigt_models/Fvoigt_v4.7.22.txt")
# files = np.append(files, "/global/homes/t/tetourne/Git/picca/py/picca/fitter2/models/fvoigt_models/Fvoigt_v4.7.22_lowcut_nhi20.3.txt")
# files = np.append(files, "/global/homes/t/tetourne/Git/picca/py/picca/fitter2/models/fvoigt_models/Fvoigt_v4.7.22.txt")
# labels = [r'$\log n_{\mathrm{HI}} = 17.6$', r'$\log n_{\mathrm{HI}} = 18.2$', r'$\log n_{\mathrm{HI}} = 18.8$', r'$\log n_{\mathrm{HI}} = 19.4$', r'$\log n_{\mathrm{HI}} = 20$', r'$\log n_{\mathrm{HI}} < 20.3$', r'$\log n_{\mathrm{HI}} > 20.3$', r'$17.2 < \log n_{\mathrm{HI}} < 22.5$']
labels = [r'$\log n_{\mathrm{HI}} = 17.6$', r'$\log n_{\mathrm{HI}} = 18.2$', r'$\log n_{\mathrm{HI}} = 18.8$', r'$\log n_{\mathrm{HI}} = 19.4$', r'$\log n_{\mathrm{HI}} = 20$', r'$17.2 < \log n_{\mathrm{HI}} < 20.3$', r'$17.2 < \log n_{\mathrm{HI}} < 22.5$']
labels = np.array(labels)[-2:]
# l_hcd = [0.22, 0.43, 0.75, 1.45, 2.80, 0.74]
# l_hcd = [0.50, 0.79, 1.46, 2.77, 5.36, 2.]
l_hcd = [0.28, 0.6, 1.2, 2.4, 4.8, 3., 6.5]
l_hcd = np.array(l_hcd)[-2:]

fig, ax = plt.subplots(ncols=2)
for i,f in enumerate(files):
    data = np.loadtxt(f).T
    k, fk = data[0], data[1]
    r, fr = powerspectrum.xi_from_pk_1D(k, fk)
    idx = np.argsort(np.abs(fr-0.5))[0]
    print(r[idx])
    ax[0].plot(k, fk, label=labels[i], color=constant.colors[i])
    ax[1].plot(r, fr/fr[0], label=labels[i], color=constant.colors[i])
    fk_rogers = np.exp(-l_hcd[i]*k)
    # fk_rogers = np.exp(-10*k)
    r, fr_rogers = powerspectrum.xi_from_pk_1D(k, fk_rogers)
    ax[0].plot(k, fk_rogers, color=constant.colors[i], linestyle='--')
    ax[1].plot(r, fr_rogers/fr_rogers[0], color=constant.colors[i], linestyle='--')

ax[1].axhline(0.5, linestyle='--', color='grey')
ax[0].set_xlim(0, 0.4)
# ax[1].set_xlim(0, 10)
ax[0].set_xlabel(r'$k \; [h\;\mathrm{Mpc}^{-1}]$')
ax[1].set_xlabel(r'$r \; [h^{-1}\mathrm{Mpc}]$')
ax[0].set_ylabel(r'$F_{\mathrm{HCD}}$')
ax[1].set_ylabel(r'$F_{\mathrm{HCD}}$')
ax[0].grid()
ax[1].grid()
ax[0].legend()
ax[1].legend()
fig.tight_layout()
plt.show()
