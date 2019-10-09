#!/usr/bin/env python
from __future__ import division, print_function
from SaclayMocks import powerspectrum
import os
import fitsio
import numpy as np
import argparse
import matplotlib.pyplot as plt
from scipy import interpolate


def P_RSD(k_par,k_perp,P, beta):
    k=np.sqrt(k_par*k_par+k_perp*k_perp)
    mu = k_par/np.maximum(k,1E-15)
    #print("k_par=",k_par)
    #print ("k_perp=",k_perp)
    #print("k=",k)
    #print("mu=",mu)
    #print ((1+beta*mu**2)**2 * P(k))
    #mu = k_perp/np.maximum(k,1E-15) # prov
    return (1+beta*mu**2)**2 * P(k)


# def PW2(k, DX) :
#     W = np.exp(- DX*DX*k*k/2)
#     return W*W*P_camb.P(k)


# def Pcut(k, kny):
#     return P_camb.P(k)*(k<kny)


def PW2cut(k):
    W = np.exp(- DX*DX*k*k/2)
    return W*W*P_camb.P(k)*(k<kny)


#def PcambFct(k) :
#    return P_camb.P(k)


#def W2Fct(k) :
#    W = (DX/np.sqrt(2*PI))* np.exp(- DX*DX*k*k/2)
#    W = W/W[0]
#    return W*W


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--out-file", type=str, default=None, required=True,
    help="path to output file")

parser.add_argument("--beta", type=float, required=True,
    help="value of beta_lya")

parser.add_argument("--voxel-size", type=float, default=2.19, required=False,
    help="value of voxel-size (delta large scale)")

parser.add_argument("--k-max", type=float, default=20., required=False,
    help="kmax to compute the missing 1D power spectrum")

parser.add_argument("--dk", type=float, default=0.001, required=False,
    help="dk to compute the missing 1D power spectrum")

parser.add_argument("--plot-p1d", action='store_true', required=False,
    help="plot the missing 1D power spectrum")

parser.add_argument("--no-rsd", action='store_false', required=False,
    help="compute the missing 1D power spectrum without RSD")

args = parser.parse_args()

outfile = args.out_file
do_plots = args.plot_p1d  # False if not specified
RSD = args.no_rsd  # True if not specified
DX = args.voxel_size
beta = args.beta
kmax = args.k_max
dk = args.dk
print("The missing 1D power spectrum will be saved in {}".format(outfile))
print("It will be computed with "+
"RSD={} and beta={}; voxcel={}; kmax={}; dk={}".format(RSD, beta, DX, kmax, dk))

# dk=0.01;  # prov
PI = np.pi
kk=np.arange(kmax/dk)*dk
filename = os.path.expandvars("$SACLAYMOCKS_BASE/etc/PlanckDR12.fits")
P_camb = powerspectrum.P_0(filename)
Pcamb = P_camb.P(kk)
kny = PI/DX
# kp = np.arange(kny/dk)*dk
kp = np.arange(PI/0.2/dk)*dk  # k for plot

#.................................     compute P1D
P1Dcamb = powerspectrum.P_1D(kk,Pcamb).P1D(kp)

if RSD:
    k_par = np.arange(kmax/dk)*dk
    k_par_t = k_par.reshape(len(k_par),-1)
    k_perp = np.arange(kmax/dk)*dk
    k_perp_t = k_perp.reshape(len(k_perp),-1)
    PRSD = P_RSD(k_par_t,k_perp,P_camb.P, beta)
    P1DcambRSD = powerspectrum.P_1D_RSD(k_par,k_perp,PRSD).P1D(kp)

#print ("0")
#plt.show()
#plt.plot(k,P1DcambRSD/P1Dcamb)
#plt.show()

# Following lines are kept just to check the intermediate steps
# #.................................     compute P1D with w^2
# # once we have W^2, the effect of k_ny cut is negligible
# # so, indistinguishable from cut at k_N and W^2
# if(False) :
#     W = np.exp(- DX*DX*kk*kk/2)
#     P1DWcamb = powerspectrum.P_1D(kk,Pcamb*W*W).P1D(kp)
#     plt.plot(kp,P1DWcamb,color="blue")
    
#     if (RSD) :
#         PRSD = P_RSD(k_par_t,k_perp,PW2)
#         P1DWcambRSD = powerspectrum.P_1D_RSD(k_par,k_perp,PRSD).P1D(kp)
#         plt.plot(kp,P1DWcambRSD,color="black")
#     print ("1")

# #.................................     compute P1D with cut at k_Nyquist
# cut = (kk<=kny)
# kk_cut=kk[cut]
# Pcamb_cut=Pcamb[cut]
# if (False):
#     P1Dcutcamb = powerspectrum.P_1D(kk_cut,Pcamb_cut).P1D(kp)
#     plt.plot(kp,P1Dcutcamb,color="blue")
#     #plt.plot(kp,P1Dcutcamb,color="green")

#     if (RSD) :
#         PRSD = P_RSD(k_par_t,k_perp,Pcut)
#         P1DcutcambRSD = powerspectrum.P_1D_RSD(k_par,k_perp,PRSD).P1D(kp)
#         plt.plot(kp,P1DcutcambRSD,color="green")
#     print ("2")

#.................................     compute P1D with cut at k_N and W^2
cut = kk <= kny
kk_cut = kk[cut]
Pcamb_cut = Pcamb[cut]
W = np.exp(- DX*DX*kk_cut*kk_cut/2)
P1DWcutcamb = powerspectrum.P_1D(kk_cut,Pcamb_cut*W*W).P1D(kp)

if RSD:
    PRSD = P_RSD(k_par_t,k_perp,PW2cut, beta)
    P1DWcutcambRSD = powerspectrum.P_1D_RSD(k_par,k_perp,PRSD).P1D(kp)

#.................................      missing P^1D(k)
P1Dmissing = np.maximum(interpolate.InterpolatedUnivariateSpline(kp,P1Dcamb - P1DWcutcamb),0)

if RSD:
    P1DmissingRSD = np.maximum(interpolate.InterpolatedUnivariateSpline(kp,P1DcambRSD - P1DWcutcambRSD),0)

# Write to fits file
print("P1D computed. Writting file...")
outfits = fitsio.FITS(outfile, 'rw', clobber=True)
table = [kp, P1Dmissing(kp)]
names = ['k', 'P1Dmiss']
if RSD:
    table.append(P1DmissingRSD(kp))
    names.append('P1DmissRSD')
outfits.write(table, names=names, extname='P1D')
outfits[-1].write_key('beta', beta)
outfits[-1].write_key('voxel', DX)
outfits[-1].write_key('kmax', kmax)
outfits[-1].write_key('dk', dk)
outfits.close()
print("Wrote {}".format(outfile))

# Plots
if do_plots:
    print("plotting...")
    filename = os.path.expandvars("$SACLAYMOCKS_BASE/etc/pkmiss_interp.fits")
    f = fitsio.FITS(filename)
    data=f[1].read()
    kk=data["k"]
    z=data["z"]
    Pk=data["Pk"]
    z0 = 2.2
    i = np.argsort(np.abs(z-z0))[0]
    f, ax = plt.subplots()
    ax.plot(kp, P1Dcamb, color="blue", label='Pcamb')
    ax.plot(kp, P1DcambRSD, color="green", label='PcambRSD')
    ax.plot(kp, P1DWcutcamb, color="blue", label='Pcamb_Wcut')
    ax.plot(kp, P1DWcutcambRSD, color="green", label='PcambRSD_Wcut')
    ax.plot(kp, P1Dmissing(kp), color="red", label='Pmissing')
    ax.plot(kp, P1DmissingRSD(kp), color="red", label='PmissingRSD')
    ax.plot(kk[i], Pk[i], color="orange", label='Pmiss_old')
    ax.set_xlabel('k [h/Mpc]')
    ax.grid()
    ax.legend()
    plt.show()

# f = fitsio.FITS("data/pkmiss_interp.fits")
# data=f[1].read()
# kk=data["k"]
# z=data["z"]
# Pk=data["Pk"]
# plt.plot(kk[80],Pk[80])
# plt.plot(k,P1DmissingRSD(k),color="red")
# plt.show()


#plt.plot(k,P1DWcutcambRSD/P1DWcambRSD,color="black")
#plt.show()

# if (RSD) :
#     plt.plot(k,P1DmissingRSD(k)/P1Dmissing(k))
#     plt.show()
