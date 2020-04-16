import fitsio
import picca.wedgize
import argparse
import matplotlib.pyplot as plt
import numpy as np
import h5py


# Plot options
TINY_SIZE = 14
SMALL_SIZE = 18
MEDIUM_SIZE = 20
BIGGER_SIZE = 22
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=TINY_SIZE)    # legend fontsize

plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.rc('figure', figsize=(9,7))

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, nargs="*", help="list of files")
parser.add_argument("--labels", type=str, nargs='*', help='list of labels', default=None)
parser.add_argument("--rt-min", type=float, default=0.)
parser.add_argument("--rt-max", type=float, default=200.)
parser.add_argument("--nrt", type=int, default=50)
parser.add_argument("--rp-min", type=float, default=0.)
parser.add_argument("--rp-max", type=float, default=200.)
parser.add_argument("--nrp", type=int, default=50)
parser.add_argument("--mu-min", type=float, default=-1.)
parser.add_argument("--mu-max", type=float, default=1.)
parser.add_argument("--r-pow", type=int, default=2)
parser.add_argument("--error-bar", action='store_true')
parser.add_argument("--sigma", action='store_true')
parser.add_argument("--title", default="")

args = parser.parse_args()

files = args.i
labels = args.labels
rpmin = args.rp_min
rpmax = args.rp_max
nrp = args.nrp
rtmin = args.rt_min
rtmax = args.rt_max
nrt = args.nrt
mumin = args.mu_min
mumax = args.mu_max
r_pow = args.r_pow

colors = ['darkblue', 'red']
linestyles = ['-.', '--']

if labels is None:
    labels = np.arange(len(files))
fmt = ['.', '.', '.', 'x', '+', 'o', '.', 'x', 'o', '-.', ':']

mu0, mu1, mu2, mu3, mu4 = 0, 0.5, 0.8, 0.95, 1
w = picca.wedgize.wedge(mumin=mumin,mumax=mumax, rtmax=rtmax, rpmax=rpmax, rtmin=rtmin, rpmin=rpmin, nrt=nrt, nrp=nrp,absoluteMu=True)
w1 = picca.wedgize.wedge(mumin=mu0,mumax=mu1, rtmax=rtmax, rpmax=rpmax, rtmin=rtmin, rpmin=rpmin, nrt=nrt, nrp=nrp,absoluteMu=True)
w2 = picca.wedgize.wedge(mumin=mu1,mumax=mu2, rtmax=rtmax, rpmax=rpmax, rtmin=rtmin, rpmin=rpmin, nrt=nrt, nrp=nrp,absoluteMu=True)
w3 = picca.wedgize.wedge(mumin=mu2,mumax=mu3, rtmax=rtmax, rpmax=rpmax, rtmin=rtmin, rpmin=rpmin, nrt=nrt, nrp=nrp,absoluteMu=True)
w4 = picca.wedgize.wedge(mumin=mu3,mumax=mu4, rtmax=rtmax, rpmax=rpmax, rtmin=rtmin, rpmin=rpmin, nrt=nrt, nrp=nrp,absoluteMu=True)

# Read file 1
data1 = fitsio.FITS(files[0])
da1 = data1[1]['DA'][:]
co1 = data1[1]['CO'][:]

# Read file 2 and +
da2 = []
co2 = []
for f in files[1:]:
    if '.fits' in f:
        data = fitsio.FITS(f)
        da2.append(data[1]['DA'][:])
        co2.append(data[1]['CO'][:])
    
    if '.h5' in f:
        ff = h5py.File(f, 'r')
        if "LYA(LYA)xLYA(LYA)" in ff.keys():
            da2.append(ff["LYA(LYA)xLYA(LYA)/fit"][...])
        elif "cf_z_0_10" in ff.keys():
            da2.append(ff["cf_z_0_10/fit"][...])
        else:
            idx1 = int(f.rfind("/"))+1
            idx2 = int(f.find(".h5"))
            da2.append(ff[f[idx2:idx2]+"/fit"][...])
        ff.close()
        co2.append(np.copy(co1))

# Compute wedges for file 1
data_wedge01 = w.wedge(da1,co1)
coef01 = data_wedge01[0]**r_pow
data_wedge11 = w1.wedge(da1,co1)
coef11 = data_wedge11[0]**r_pow
data_wedge21 = w2.wedge(da1,co1)
coef21 = data_wedge21[0]**r_pow
data_wedge31 = w3.wedge(da1,co1)
coef31 = data_wedge31[0]**r_pow
data_wedge41 = w4.wedge(da1,co1)
coef41 = data_wedge41[0]**r_pow

# Compute wedges for file 2 and +
data_wedge02 = []
data_wedge12 = []
data_wedge22 = []
data_wedge32 = []
data_wedge42 = []
coef02 = []
coef12 = []
coef22 = []
coef32 = []
coef42 = []
for da, co in zip(da2, co2):
    data_wedge02.append(w.wedge(da,co))
    data_wedge12.append(w1.wedge(da,co))
    data_wedge22.append(w2.wedge(da,co))
    data_wedge32.append(w3.wedge(da,co))
    data_wedge42.append(w4.wedge(da,co))
    coef02.append(data_wedge02[-1][0]**r_pow)
    coef12.append(data_wedge12[-1][0]**r_pow)
    coef22.append(data_wedge22[-1][0]**r_pow)
    coef32.append(data_wedge32[-1][0]**r_pow)
    coef42.append(data_wedge42[-1][0]**r_pow)

f, axs = plt.subplots(5, sharex=True, figsize=(12,8))
for i in range(len(files[1:])):
    if args.sigma:
        axs[0].plot(data_wedge01[0], (data_wedge01[1] - data_wedge02[i][1]) / np.sqrt(np.diag(data_wedge01[2])), color=colors[i], linestyle=linestyles[i], label=labels[i])
        axs[1].plot(data_wedge11[0], (data_wedge11[1] - data_wedge12[i][1]) / np.sqrt(np.diag(data_wedge11[2])), color=colors[i], linestyle=linestyles[i])
        axs[2].plot(data_wedge21[0], (data_wedge21[1] - data_wedge22[i][1]) / np.sqrt(np.diag(data_wedge21[2])), color=colors[i], linestyle=linestyles[i])
        axs[3].plot(data_wedge31[0], (data_wedge31[1] - data_wedge32[i][1]) / np.sqrt(np.diag(data_wedge31[2])), color=colors[i], linestyle=linestyles[i])
        axs[4].plot(data_wedge41[0], (data_wedge41[1] - data_wedge42[i][1]) / np.sqrt(np.diag(data_wedge41[2])), color=colors[i], linestyle=linestyles[i])
    elif args.error_bar:
        axs[0].errorbar(data_wedge01[0], (data_wedge01[1] - data_wedge02[i][1])*coef01, yerr=np.sqrt(np.diag(data_wedge01[2]))*coef01, color=colors[i], linestyle=linestyles[i], label=labels[i])
        axs[1].errorbar(data_wedge11[0], (data_wedge11[1] - data_wedge12[i][1])*coef11, yerr=np.sqrt(np.diag(data_wedge11[2]))*coef11, color=colors[i], linestyle=linestyles[i])
        axs[2].errorbar(data_wedge21[0], (data_wedge21[1] - data_wedge22[i][1])*coef21, yerr=np.sqrt(np.diag(data_wedge21[2]))*coef21, color=colors[i], linestyle=linestyles[i])
        axs[3].errorbar(data_wedge31[0], (data_wedge31[1] - data_wedge32[i][1])*coef31, yerr=np.sqrt(np.diag(data_wedge31[2]))*coef31, color=colors[i], linestyle=linestyles[i])
        axs[4].errorbar(data_wedge41[0], (data_wedge41[1] - data_wedge42[i][1])*coef41, yerr=np.sqrt(np.diag(data_wedge41[2]))*coef41, color=colors[i], linestyle=linestyles[i])
    else:
        axs[0].plot(data_wedge01[0], (data_wedge01[1] - data_wedge02[i][1])*coef01, color=colors[i], linestyle=linestyles[i], label=labels[i])
        axs[1].plot(data_wedge11[0], (data_wedge11[1] - data_wedge12[i][1])*coef11, color=colors[i], linestyle=linestyles[i])
        axs[2].plot(data_wedge21[0], (data_wedge21[1] - data_wedge22[i][1])*coef21, color=colors[i], linestyle=linestyles[i])
        axs[3].plot(data_wedge31[0], (data_wedge31[1] - data_wedge32[i][1])*coef31, color=colors[i], linestyle=linestyles[i])
        axs[4].plot(data_wedge41[0], (data_wedge41[1] - data_wedge42[i][1])*coef41, color=colors[i], linestyle=linestyles[i])

for i, ax in enumerate(axs):
    ax.axhline(0, linestyle='--', color='black')
    ax.grid()
    if args.sigma:
        ax.set_ylim(-4,4)

axs[0].set_title(r'${}<\mu<{}$'.format(mu0,mu4), fontsize=TINY_SIZE, loc='right')
axs[1].set_title(r'${}<\mu<{}$'.format(mu0,mu1), fontsize=TINY_SIZE, loc='right')
axs[2].set_title(r'${}<\mu<{}$'.format(mu1,mu2), fontsize=TINY_SIZE, loc='right')
axs[3].set_title(r'${}<\mu<{}$'.format(mu2,mu3), fontsize=TINY_SIZE, loc='right')
axs[4].set_title(r'${}<\mu<{}$'.format(mu3,mu4), fontsize=TINY_SIZE, loc='right')

# axs[0].text(0., 1, r'${}<\mu<{}$'.format(mu0,mu4), fontsize=TINY_SIZE, verticalalignment='top')
# axs[1].text(0.2, 0.65, r'${}<\mu<{}$'.format(mu0,mu1), fontsize=TINY_SIZE, verticalalignment='top')
# axs[2].text(0.4, 0.35, r'${}<\mu<{}$'.format(mu1,mu2), fontsize=TINY_SIZE, verticalalignment='top')
# axs[3].text(0.6, 0.15, r'${}<\mu<{}$'.format(mu2,mu3), fontsize=TINY_SIZE, verticalalignment='top')
# axs[4].text(0.8, 0.05, r'${}<\mu<{}$'.format(mu3,mu4), fontsize=TINY_SIZE, verticalalignment='top')

axs[4].set_xlabel('r Mpc/h')
# axs[0].set_title(args.title)
axs[0].legend()
plt.show()
