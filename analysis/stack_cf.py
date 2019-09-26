import glob
import fitsio
import picca.wedgize
import h5py
import matplotlib.pyplot as plt
import numpy as np


plot_fit = True
mumin=-1
mumax=1
rtmin=0
rtmax=200
rpmin=-200
rpmax=200
nrt=50
nrp=100
mu0, mu1, mu2, mu3 = 0, 0.5, 0.8, 1
r_pow = 2
title = 'stack of 10 Lya-QSO cross-correlations'

files = glob.glob("/global/projecta/projectdirs/desi/mocks/lya_forest/picca/saclay/v4.4/v4.4.*/eboss-0.0/xcf_z_0_10-exp.fits")

data_wedge1_stack = []
data_wedge2_stack = []
data_wedge3_stack = []
f1_stack = []
f2_stack = []
f3_stack = []
fig, ax = plt.subplots(figsize=(12,8))
for f in files:
    # Read data
    data = fitsio.FITS(f)
    da = data[1].read()['DA']
    co = data[1].read()['CO']
    data.close()

    # Read fit
    if plot_fit:
        fit_file = f.replace('.fits', '.h5')
        ff = h5py.File(fit_file)

    w1 = picca.wedgize.wedge(mumin=mu0,mumax=mu1, rtmax=rtmax, rpmax=rpmax, rtmin=rtmin, rpmin=rpmin, nrt=nrt, nrp=nrp)
    w2 = picca.wedgize.wedge(mumin=mu1,mumax=mu2, rtmax=rtmax, rpmax=rpmax, rtmin=rtmin, rpmin=rpmin, nrt=nrt, nrp=nrp)
    w3 = picca.wedgize.wedge(mumin=mu2,mumax=mu3, rtmax=rtmax, rpmax=rpmax, rtmin=rtmin, rpmin=rpmin, nrt=nrt, nrp=nrp)
    data_wedge1 = w1.wedge(da,co)
    coef1 = data_wedge1[0]**r_pow
    data_wedge2 = w2.wedge(da,co)
    coef2 = data_wedge2[0]**r_pow
    data_wedge3 = w3.wedge(da,co)
    coef3 = data_wedge3[0]**r_pow
    if plot_fit:
        try:
            r1,f1,_ = w1.wedge(ff["LYA(LYA)xLYA(LYA)/fit"][...],co)
            r2,f2,_ = w2.wedge(ff["LYA(LYA)xLYA(LYA)/fit"][...],co)
            r3,f3,_ = w3.wedge(ff["LYA(LYA)xLYA(LYA)/fit"][...],co)
        except KeyError:
            print("Can't find LYA(LYA)xLYA(LYA)/fit")
            try:
                i = int(fit_file.find(".h5"))
                j = int(fit_file.rfind("/"))+1
                r1,f1,_ = w1.wedge(ff[fit_file[j:i]+"/fit"][...],co)
                r2,f2,_ = w2.wedge(ff[fit_file[j:i]+"/fit"][...],co)
                r3,f3,_ = w3.wedge(ff[fit_file[j:i]+"/fit"][...],co)
            except:
                print("Can't find {}".format(fit_file[j:i]+"/fit"))
    
    ax.errorbar(data_wedge1[0],coef1*data_wedge1[1],yerr=coef1*np.sqrt(np.diag(data_wedge1[2])), color='b', alpha=0.4)
    ax.errorbar(data_wedge2[0],coef2*data_wedge2[1],yerr=coef2*np.sqrt(np.diag(data_wedge2[2])), color='g', alpha=0.4)
    ax.errorbar(data_wedge3[0],coef3*data_wedge3[1],yerr=coef3*np.sqrt(np.diag(data_wedge3[2])), color='r', alpha=0.4)
    # ax.plot(r1, f1*r1**r_pow, color='b')
    # ax.plot(r2, f2*r2**r_pow, color='g')
    # ax.plot(r3, f3*r3**r_pow, color='r')

    # Compute the stack
    data_wedge1_stack.append([data_wedge1[0], data_wedge1[1], np.diag(data_wedge1[2])])
    data_wedge2_stack.append([data_wedge2[0], data_wedge2[1], np.diag(data_wedge2[2])])
    data_wedge3_stack.append([data_wedge3[0], data_wedge3[1], np.diag(data_wedge3[2])])
    if plot_fit:
        f1_stack.append(f1)
        f2_stack.append(f2)
        f3_stack.append(f3)

data_wedge1_stack = np.sum(data_wedge1_stack, axis=0) / len(files)
data_wedge2_stack = np.sum(data_wedge2_stack, axis=0) / len(files)
data_wedge3_stack = np.sum(data_wedge3_stack, axis=0) / len(files)
f1_stack = np.sum(f1_stack, axis=0) / len(files)
f2_stack = np.sum(f2_stack, axis=0) / len(files)
f3_stack = np.sum(f3_stack, axis=0) / len(files)

# Plot stack
ax.errorbar(data_wedge1_stack[0],coef1*data_wedge1_stack[1],yerr=coef1*np.sqrt(data_wedge1_stack[2]),fmt='--', label=r"${}<\mu<{}$".format(mu0, mu1), color='b')
ax.errorbar(data_wedge2_stack[0],coef2*data_wedge2_stack[1],yerr=coef2*np.sqrt(data_wedge2_stack[2]),fmt='--', label=r"${}<\mu<{}$".format(mu1, mu2), color='g')
ax.errorbar(data_wedge3_stack[0],coef3*data_wedge3_stack[1],yerr=coef3*np.sqrt(data_wedge3_stack[2]),fmt='--', label=r"${}<\mu<{}$".format(mu2, mu3), color='r')
ax.plot(r1, f1_stack*r1**r_pow, color='b')
ax.plot(r2, f2_stack*r2**r_pow, color='g')
ax.plot(r3, f3_stack*r3**r_pow, color='r')

ax.grid()
ax.set_title(title, fontsize=20)
ax.set_xlabel(r"$r \, [\mathrm{Mpc \, h^{-1}}]$",fontsize=20)
if r_pow == 2:
    ax.set_ylabel(r"$r^{2}\xi(r) \, [\mathrm{Mpc \, h^{-1}}]$",fontsize=20)
if r_pow == 1:
    ax.set_ylabel(r"$r\xi(r) \, [\mathrm{Mpc \, h^{-1}}]$",fontsize=20)
if r_pow == 0:
    ax.set_ylabel(r"$\xi(r) \, [\mathrm{Mpc \, h^{-1}}]$",fontsize=20)

ax.legend()
plt.show()
