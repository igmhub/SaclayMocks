import sys
import matplotlib.pyplot as plt
import fitsio
import picca.wedgize
import scipy as sp
import argparse
from LyaMocks import util

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str)
parser.add_argument("-pred", type=str, default=None)
parser.add_argument("--plot-pred", type=str, default="False")
parser.add_argument("--rt-max", type=float, default=200)
parser.add_argument("--rp-max", type=float, default=200)
parser.add_argument("--r-pow", type=int, default=2)
parser.add_argument("--nrt", type=int, default=50)
parser.add_argument("--nrp", type=int, default=50)
parser.add_argument("--title", default="")
args = parser.parse_args()

# Params
infile = args.i
rpmax = args.rp_max
rtmax = args.rt_max
nrt=args.nrt
nrp=args.nrp
r_pow = args.r_pow
pred = util.str2bool(args.plot_pred)
title = args.title

# Load CF
print("Reading {}".format(infile))
print("rt_max : {}  rp_max : {}".format(rtmax, rpmax))
data = fitsio.FITS(infile)
da = data[1]['DA'][:]
co = data[1]['CO'][:]
data.close()
data2 = fitsio.FITS("Out/debug_v4.6_22/from_transmissions/Correlations/e_cf.fits") # prov
da2 = data2[1].read()['DA']
co2 = data2[1].read()['CO']
data2.close()
da = da - da2

# Load prediction
if pred:
    infile = "../../LyaMocks/bin/data/xi_FGPA_nosc_pred.fits"
    data = fitsio.read(infile, ext=1)
    rpred = data["r"]
    xipred = 5*data["xi"]
    xiprederr = 5*data["xierr"]

# Load xi pred
if args.pred is not None:
    infile = args.pred
    data_pred = fitsio.FITS(infile)
    da_pred = data_pred[1]['DA'][:]
    # da_pred *= 19 / 13.  #prov
    co_pred = data_pred[1]['CO'][:]
    data_pred.close()

# Plot Xi0 and Xipred
w = picca.wedgize.wedge(mumin=0.,mumax=1., rtmax=rtmax, rpmax=rpmax, nrt=nrt, nrp=nrp)
data_wedge = w.wedge(da,co)
coef = data_wedge[0]**r_pow
# data_wedge2 = w.wedge(da2,co2)  # prov
if args.pred is not None:
    data_pred_wedge = w.wedge(da_pred, co_pred)
    coef_pred = data_pred_wedge[0]**r_pow

f, ax = plt.subplots(figsize=(12,8))
ax.errorbar(data_wedge[0],coef*data_wedge[1],yerr=coef*sp.sqrt(sp.diag(data_wedge[2])),marker='+', label='mock', color='b')
# ax.errorbar(data_wedge2[0],coef*data_wedge2[1],yerr=coef*sp.sqrt(sp.diag(data_wedge2[2])),marker='+', label='v4.6_22')  # prov
if pred:
    ax.errorbar(rpred, rpred**2 * xipred, yerr=rpred**2 * xiprederr, marker='.', label='prediction')
if args.pred is not None:
    ax.plot(data_pred_wedge[0],coef_pred*data_pred_wedge[1], '--', label='prediction', color='b')
ax.grid()
ax.legend()
ax.set_title("{}".format(title), fontsize=20)
ax.set_xlabel(r"$r \, [\mathrm{Mpc \, h^{-1}}]$",fontsize=20)
if r_pow == 2:
    ax.set_ylabel(r"$r^{2}\xi(r) \, [\mathrm{Mpc \, h^{-1}}]$",fontsize=20)
if r_pow == 1:
    ax.set_ylabel(r"$r\xi(r) \, [\mathrm{Mpc \, h^{-1}}]$",fontsize=20)
if r_pow == 0:
    ax.set_ylabel(r"$\xi(r) \, [\mathrm{Mpc \, h^{-1}}]$",fontsize=20)

# Plot wedges
mu0, mu1, mu2, mu3, mu4 = 0, 0.5, 0.8, 0.95, 1
w1 = picca.wedgize.wedge(mumin=mu0,mumax=mu1, rtmax=rtmax, rpmax=rpmax, nrt=nrt, nrp=nrp)
w2 = picca.wedgize.wedge(mumin=mu1,mumax=mu2, rtmax=rtmax, rpmax=rpmax, nrt=nrt, nrp=nrp)
w3 = picca.wedgize.wedge(mumin=mu2,mumax=mu3, rtmax=rtmax, rpmax=rpmax, nrt=nrt, nrp=nrp)
w4 = picca.wedgize.wedge(mumin=mu3,mumax=mu4, rtmax=rtmax, rpmax=rpmax, nrt=nrt, nrp=nrp)
data_wedge1 = w1.wedge(da,co)
coef1 = data_wedge1[0]**r_pow
data_wedge2 = w2.wedge(da,co)
coef2 = data_wedge2[0]**r_pow
data_wedge3 = w3.wedge(da,co)
coef3 = data_wedge3[0]**r_pow
data_wedge4 = w4.wedge(da,co)
coef4 = data_wedge4[0]**r_pow
# data_wedge1_2 = w1.wedge(da2,co2)  # prov
# data_wedge2_2 = w2.wedge(da2,co2)
# data_wedge3_2 = w3.wedge(da2,co2)
# data_wedge4_2 = w4.wedge(da2,co2)
if args.pred is not None:
    data_wedge1_pred = w1.wedge(da_pred,co_pred)
    coef1_pred = data_wedge1_pred[0]**r_pow
    data_wedge2_pred = w2.wedge(da_pred,co_pred)
    coef2_pred = data_wedge2_pred[0]**r_pow
    data_wedge3_pred = w3.wedge(da_pred,co_pred)
    coef3_pred = data_wedge3_pred[0]**r_pow
    data_wedge4_pred = w4.wedge(da_pred,co_pred)
    coef4_pred = data_wedge4_pred[0]**r_pow

f, ax = plt.subplots(figsize=(12,8))
ax.errorbar(data_wedge1[0],coef1*data_wedge1[1],yerr=coef1*sp.sqrt(sp.diag(data_wedge1[2])), marker='+', label=r"${}<\mu<{}$".format(mu0, mu1), color='b')
ax.errorbar(data_wedge2[0],coef2*data_wedge2[1],yerr=coef2*sp.sqrt(sp.diag(data_wedge2[2])), marker='+', label=r"${}<\mu<{}$".format(mu1, mu2), color='g')
ax.errorbar(data_wedge3[0],coef3*data_wedge3[1],yerr=coef3*sp.sqrt(sp.diag(data_wedge3[2])), marker='+', label=r"${}<\mu<{}$".format(mu2, mu3), color='orange')
ax.errorbar(data_wedge4[0],coef4*data_wedge4[1],yerr=coef4*sp.sqrt(sp.diag(data_wedge4[2])), marker='+', label=r"${}<\mu<{}$".format(mu3, mu4), color='red')
# ax.errorbar(data_wedge1_2[0],coef1*data_wedge1_2[1],yerr=coef1*sp.sqrt(sp.diag(data_wedge1_2[2])), marker='+', label=r"${}<\mu<{}$".format(mu0, mu1), color='b')  # prov
# ax.errorbar(data_wedge2_2[0],coef2*data_wedge2_2[1],yerr=coef2*sp.sqrt(sp.diag(data_wedge2_2[2])), marker='+', label=r"${}<\mu<{}$".format(mu1, mu2), color='g')
# ax.errorbar(data_wedge3_2[0],coef3*data_wedge3_2[1],yerr=coef3*sp.sqrt(sp.diag(data_wedge3_2[2])), marker='+', label=r"${}<\mu<{}$".format(mu2, mu3), color='orange')
# ax.errorbar(data_wedge4_2[0],coef4*data_wedge4_2[1],yerr=coef4*sp.sqrt(sp.diag(data_wedge4_2[2])), marker='+', label=r"${}<\mu<{}$".format(mu3, mu4), color='red')
if args.pred is not None:
    ax.plot(data_wedge1_pred[0],coef1_pred*data_wedge1_pred[1], '--', color='b')
    ax.plot(data_wedge2_pred[0],coef2_pred*data_wedge2_pred[1], '--',color='g')
    ax.plot(data_wedge3_pred[0],coef3_pred*data_wedge3_pred[1], '--',color='orange')
    ax.plot(data_wedge4_pred[0],coef4_pred*data_wedge4_pred[1], '--',color='red')
ax.grid()
ax.set_title("{}".format(title), fontsize=20)
ax.set_xlabel(r"$r \, [\mathrm{Mpc \, h^{-1}}]$",fontsize=20)
if r_pow == 2:
    ax.set_ylabel(r"$r^{2}\xi(r) \, [\mathrm{Mpc \, h^{-1}}]$",fontsize=20)
if r_pow == 1:
    ax.set_ylabel(r"$r\xi(r) \, [\mathrm{Mpc \, h^{-1}}]$",fontsize=20)
if r_pow == 0:
    ax.set_ylabel(r"$\xi(r) \, [\mathrm{Mpc \, h^{-1}}]$",fontsize=20)

ax.legend()

plt.show()
