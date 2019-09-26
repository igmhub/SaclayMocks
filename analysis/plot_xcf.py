import sys
import matplotlib.pyplot as plt
import fitsio
import picca.wedgize
import scipy as sp
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str)
parser.add_argument("--rt-min", type=float, default=0.)
parser.add_argument("--rt-max", type=float, default=200.)
parser.add_argument("--nrt", type=int, default=50)
parser.add_argument("--rp-min", type=float, default=-200.)
parser.add_argument("--rp-max", type=float, default=200.)
parser.add_argument("--nrp", type=int, default=100)
parser.add_argument("--mu-min", type=float, default=-1.)
parser.add_argument("--mu-max", type=float, default=1.)
parser.add_argument("--r-pow", type=int, default=2)
parser.add_argument("--title", default="")
args = parser.parse_args()

# Params
infile = args.i
rpmin = args.rp_min
rpmax = args.rp_max
nrp = args.nrp
rtmin = args.rt_min
rtmax = args.rt_max
nrt = args.nrt
mumin = args.mu_min
mumax = args.mu_max
r_pow = args.r_pow
title = args.title

# Load CF
print("Reading {}".format(infile))
print("rt_max : {}  rp_max : {}".format(rtmax, rpmax))
data = fitsio.FITS(infile)
da = data[1]['DA'][:]
co = data[1]['CO'][:]
data.close()

# Plot Xi0 and Xipred
w = picca.wedgize.wedge(rpmin=rpmin,rpmax=rpmax,nrp=nrp,rtmin=rtmin,rtmax=rtmax,nrt=nrt,mumin=mumin,mumax=mumax)
data_wedge = w.wedge(da,co)
coef = data_wedge[0]**r_pow

f, ax = plt.subplots(figsize=(12,8))
ax.errorbar(data_wedge[0],coef*data_wedge[1],yerr=coef*sp.sqrt(sp.diag(data_wedge[2])),marker='+', label='mock')
ax.grid()
ax.set_title("{}".format(title), fontsize=20)
ax.set_xlabel(r"$r \, [\mathrm{Mpc \, h^{-1}}]$",fontsize=20)
if r_pow == 2:
    ax.set_ylabel(r"$r^{2}\xi(r) \, [\mathrm{Mpc \, h^{-1}}]$",fontsize=20)
if r_pow == 1:
    ax.set_ylabel(r"$r\xi(r) \, [\mathrm{Mpc \, h^{-1}}]$",fontsize=20)
if r_pow == 0:
    ax.set_ylabel(r"$\xi(r) \, [\mathrm{Mpc \, h^{-1}}]$",fontsize=20)

# Plot wedges
mu0, mu1, mu2, mu3 = 0, 0.5, 0.8, 1
w1 = picca.wedgize.wedge(rpmin=rpmin,rpmax=rpmax,nrp=nrp,rtmin=rtmin,rtmax=rtmax,nrt=nrt,mumin=mu0,mumax=mu1)
w2 = picca.wedgize.wedge(rpmin=rpmin,rpmax=rpmax,nrp=nrp,rtmin=rtmin,rtmax=rtmax,nrt=nrt,mumin=mu1,mumax=mu2)
w3 = picca.wedgize.wedge(rpmin=rpmin,rpmax=rpmax,nrp=nrp,rtmin=rtmin,rtmax=rtmax,nrt=nrt,mumin=mu2,mumax=mu3)
data_wedge1 = w1.wedge(da,co)
coef1 = data_wedge1[0]**r_pow
data_wedge2 = w2.wedge(da,co)
coef2 = data_wedge2[0]**r_pow
data_wedge3 = w3.wedge(da,co)
coef3 = data_wedge3[0]**r_pow

f, ax = plt.subplots(figsize=(12,8))
ax.errorbar(data_wedge1[0],coef1*data_wedge1[1],yerr=coef1*sp.sqrt(sp.diag(data_wedge1[2])),marker='.', label=r"${}<\mu<{}$".format(mu0, mu1))
ax.errorbar(data_wedge2[0],coef2*data_wedge2[1],yerr=coef2*sp.sqrt(sp.diag(data_wedge2[2])),marker='.', label=r"${}<\mu<{}$".format(mu1, mu2))
ax.errorbar(data_wedge3[0],coef3*data_wedge3[1],yerr=coef3*sp.sqrt(sp.diag(data_wedge3[2])),marker='.', label=r"${}<\mu<{}$".format(mu2, mu3))
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
