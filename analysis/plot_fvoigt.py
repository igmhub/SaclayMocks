import glob
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from SaclayMocks import util, powerspectrum, constant
from iminuit import Minuit


SMALL_SIZE = 15
MEDIUM_SIZE = 18
BIGGER_SIZE = 16
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize

plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.rc('figure', figsize=(12,8))

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
l_hcd = [0.28, 0.6, 1.2, 2.5, 5, 2.8, 6.5]
l_hcd = np.array(l_hcd)[-2:]

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
for i,f in enumerate(files):
    data = np.loadtxt(f).T
    k, fk = data[0], data[1]
    r, fr = powerspectrum.xi_from_pk_1D(k, fk)
    idx = np.argsort(np.abs(fr-0.5))[0]
    print(r[idx])
    ax1.plot(k, fk, label=labels[i], color=constant.colors[i])
    ax2.plot(r, fr/fr[0], label=labels[i], color=constant.colors[i])
    fk_rogers = np.exp(-l_hcd[i]*k)
    # fk_rogers = np.exp(-10*k)
    r, fr_rogers = powerspectrum.xi_from_pk_1D(k, fk_rogers)
    ax1.plot(k, fk_rogers, color=constant.colors[i], linestyle='--')
    ax2.plot(r, fr_rogers/fr_rogers[0], color=constant.colors[i], linestyle='--')

fk_rogers = np.exp(-10*k)
r, fr_rogers = powerspectrum.xi_from_pk_1D(k, fk_rogers)
ax1.plot(k, fk_rogers, color=constant.colors[i+1], linestyle='--', label=r'$L_{\mathrm{HCD}}=10$')
ax2.plot(r, fr_rogers/fr_rogers[0], color=constant.colors[i+1], linestyle='--', label=r'$L_{\mathrm{HCD}}=10$')

ax2.axhline(0.5, linestyle='--', color='grey')
ax1.set_xlim(0, 0.8)
# ax2.set_xlim(0, 10)
ax1.set_xlabel(r'$k \; [h\;\mathrm{Mpc}^{-1}]$')
ax2.set_xlabel(r'$r \; [h^{-1}\mathrm{Mpc}]$')
ax1.set_ylabel(r'$F_{\mathrm{HCD}}$')
ax2.set_ylabel(r'$F_{\mathrm{HCD}}$')
ax1.grid()
ax2.grid()
ax1.legend()
ax2.legend()
fig1.tight_layout()
fig2.tight_layout()
plt.show()
