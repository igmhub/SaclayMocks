#!/usr/bin/env python
from __future__ import division, print_function
from SaclayMocks import powerspectrum, util, constant
import os
import fitsio
import numpy as np
import argparse
import matplotlib.pyplot as plt
from scipy import interpolate


def P_RSD(k_par,k_perp,P, beta):
    k=np.sqrt(k_par*k_par+k_perp*k_perp)
    mu = k_par/np.maximum(k,1E-15)
    return (1+beta*mu**2)**2 * P(k)

def D_Prats(k,mu,P,zref,Growth):
    if (abs(zref-3.0)<0.0001):
        q1 = 0.792; q2=0; k_p = 17.1; k_v = 1.16; a_v = 0.578; b_v = 1.63 # Planck Tab 7 
    if (abs(zref-2.8)<0.0001):
        q1 = 0.773; q2=0; k_p = 19.2; k_v = 1.16; a_v = 0.608; b_v = 1.65 # Planck Tab 7 
    if (abs(zref-2.6)<0.0001):
        q1 = 0.781; q2=0; k_p = 21.1; k_v = 1.15; a_v = 0.611; b_v = 1.64 # Planck Tab 7 
    if (abs(zref-2.4)<0.0001):
        q1 = 0.851; q2=0; k_p = 19.5; k_v = 1.06; a_v = 0.548; b_v = 1.61 # Planck Tab 7 
        #q1 = -0.0020; q2=0.623; k_p = 10.9; k_v = 0.517; a_v = 0.152; b_v = 1.62 # Planck Tab 5 
    if (abs(zref-2.2)<0.0001):
        q1 = 0.867; q2=0; k_p = 19.4; k_v = 1.06; a_v = 0.514; b_v = 1.60 # Planck Tab 7 
    Delta2 = k**3 * P(k) * Growth * Growth /2/np.pi**2
    return np.exp ( (q1 * Delta2 + q2 * Delta2**2) * (1-(k/k_v)**a_v * mu**b_v ) - (k/k_p)**2 )

def P_RSD_Prats(k_par,k_perp,P, beta,zref,Growth):
    k=np.sqrt(k_par*k_par+k_perp*k_perp)
    mu = k_par/np.maximum(k,1E-15)
    return (1+beta*mu**2)**2 * P(k) * D_Prats(k,mu,P,zref,Growth)


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

parser.add_argument("--zref", type=float, default=2.2, required=False,
    help="value of redshift")

parser.add_argument("--k-max", type=float, default=20., required=False,
    help="kmax to compute the missing 1D power spectrum")

parser.add_argument("--dk", type=float, default=0.001, required=False,
    help="dk to compute the missing 1D power spectrum")

parser.add_argument("--plot-p1d", action='store_true', required=False,
    help="plot the missing 1D power spectrum")

parser.add_argument("--no-rsd", action='store_false', required=False,
    help="compute the missing 1D power spectrum without RSD")

parser.add_argument("--linear", action='store_true', required=False,
    help="only linear contributions")

args = parser.parse_args()

outfile = args.out_file
do_plots = args.plot_p1d  # False if not specified
RSD = args.no_rsd  # True if not specified
DX = args.voxel_size
beta = args.beta
kmax = args.k_max
dk = args.dk
linear = args.linear # False if not specified
zref = args.zref
print("The missing 1D power spectrum will be saved in {}".format(outfile))
print("It will be computed with "+
"RSD={} and beta={}; voxcel={}; kmax={}; dk={}".format(RSD, beta, DX, kmax, dk))
if not linear:
    print("z_ref = {}".format(zref))
    if zref > 3:
        zref = 3
        print("Redshift is set to {} for all z > 3".format(zref))
PI = np.pi
kk=np.arange(kmax/dk)*dk
filename = os.path.expandvars("$SACLAYMOCKS_BASE/etc/PlanckDR12.fits")
P_camb = powerspectrum.P_0(filename)
Pcamb = P_camb.P(kk)
kny = PI/DX
# kp = np.arange(PI/0.2/dk)*dk  # k for plot

#.................................     compute P1D
P1Dcamb = powerspectrum.P_1D(kk,Pcamb).P1D(kk)

if RSD:
    k_par = np.arange(kmax/dk)*dk
    k_par_t = k_par.reshape(len(k_par),-1)
    k_perp = np.arange(kmax/dk)*dk
    k_perp_t = k_perp.reshape(len(k_perp),-1)
    if (linear) :
        PRSD = P_RSD(k_par_t,k_perp,P_camb.P, beta)
    else :
        Growth = util.fgrowth(zref, constant.omega_M_0)
        PRSD = P_RSD_Prats(k_par_t,k_perp,P_camb.P, beta, zref,Growth)
    P1DcambRSD = powerspectrum.P_1D_RSD(k_par,k_perp,PRSD).P1D(kk)

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
#     P1DWcamb = powerspectrum.P_1D(kk,Pcamb*W*W).P1D(kk)
#     plt.plot(kk,P1DWcamb,color="blue")
    
#     if (RSD) :
#         PRSD = P_RSD(k_par_t,k_perp,PW2)
#         P1DWcambRSD = powerspectrum.P_1D_RSD(k_par,k_perp,PRSD).P1D(kk)
#         plt.plot(kk,P1DWcambRSD,color="black")
#     print ("1")

# #.................................     compute P1D with cut at k_Nyquist
# cut = (kk<=kny)
# kk_cut=kk[cut]
# Pcamb_cut=Pcamb[cut]
# if (False):
#     P1Dcutcamb = powerspectrum.P_1D(kk_cut,Pcamb_cut).P1D(kk)
#     plt.plot(kk,P1Dcutcamb,color="blue")
#     #plt.plot(kk,P1Dcutcamb,color="green")

#     if (RSD) :
#         PRSD = P_RSD(k_par_t,k_perp,Pcut)
#         P1DcutcambRSD = powerspectrum.P_1D_RSD(k_par,k_perp,PRSD).P1D(kk)
#         plt.plot(kk,P1DcutcambRSD,color="green")
#     print ("2")

#.................................     compute P1D with cut at k_N and W^2
cut = kk <= kny
kk_cut = kk[cut]
Pcamb_cut = Pcamb[cut]
W = np.exp(- DX*DX*kk_cut*kk_cut/2)
P1DWcutcamb = powerspectrum.P_1D(kk_cut,Pcamb_cut*W*W).P1D(kk)

if RSD:
    PRSD = P_RSD(k_par_t,k_perp,PW2cut, beta)
    P1DWcutcambRSD = powerspectrum.P_1D_RSD(k_par,k_perp,PRSD).P1D(kk)

#.................................      missing P^1D(k)
P1Dmissing = interpolate.InterpolatedUnivariateSpline(kk, np.maximum(P1Dcamb - P1DWcutcamb, 0))

if RSD:
    P1DmissingRSD = interpolate.InterpolatedUnivariateSpline(kk, np.maximum(P1DcambRSD - P1DWcutcambRSD, 0))

# Write to fits file
print("P1D computed. Writting file...")
outfits = fitsio.FITS(outfile, 'rw', clobber=True)
table = [kk, P1Dmissing(kk)]
names = ['k', 'P1Dmiss']
if RSD:
    table.append(P1DmissingRSD(kk))
    names.append('P1DmissRSD')
outfits.write(table, names=names, extname='P1D')
outfits[-1].write_key('beta', beta)
outfits[-1].write_key('voxel', DX)
outfits[-1].write_key('kmax', kmax)
outfits[-1].write_key('dk', dk)
outfits.close()
print("Wrote {}".format(outfile))

# # Plots
# if do_plots:
#     print("plotting...")
#     filename = os.path.expandvars("$SACLAYMOCKS_BASE/etc/pkmiss_interp.fits")
#     f = fitsio.FITS(filename)
#     data=f[1].read()
#     kk=data["k"]
#     z=data["z"]
#     Pk=data["Pk"]
#     z0 = 2.2
#     i = np.argsort(np.abs(z-z0))[0]
#     f, ax = plt.subplots()
#     ax.plot(kk, P1Dcamb, color="blue", label='Pcamb')
#     ax.plot(kk, P1DcambRSD, color="green", label='PcambRSD')
#     ax.plot(kk, P1DWcutcamb, color="blue", label='Pcamb_Wcut')
#     ax.plot(kk, P1DWcutcambRSD, color="green", label='PcambRSD_Wcut')
#     ax.plot(kk, P1Dmissing(kk), color="red", label='Pmissing')
#     ax.plot(kk, P1DmissingRSD(kk), color="red", label='PmissingRSD')
#     ax.plot(kk[i], Pk[i], color="orange", label='Pmiss_old')
#     ax.set_xlabel('k [h/Mpc]')
#     ax.grid()
#     ax.legend()
#     plt.show()

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
