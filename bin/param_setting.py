#!/usr/bin/env python
# compute the FGPA parameters a,b,c and P_s(k)
#       for given bias, beta and <F>
# c= beta
# bias,beta -> P_F^3D -> P_F^1D -> sig_F^2
# <F> sig_F -> a, b*sig_g
# P_s(k) = P_{mat}^{1D} -  P_{mat,pixel,smearing}^{1D} 
# P_s -> sigma_s
# sig_g^2 = sig_l^2 + sig_s^2 + c^2 sig_eta^2 + 2 c cov(delta_L,eta) -> b

# simu1D/src/fitab.cpp
# Fmean(z) + sigF(z) -> a, b
# does not compile !
# coded in FGPApdf.py and cross ckecked with fitab.cpp written results

from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
from SaclayMocks import util
from SaclayMocks import powerspectrum
from SaclayMocks import constant
from SaclayMocks import FGPApdf as pdf


kmax = 5    # 20.
dk=0.01 # 0.001
DX = 2.19 # Mpc/h voxel size
Dperp = 2.19 # effective transverse pixel size 
# need a study to get best value of Dperp and how sensitive is the final P1Dmissing <==
# or may be we need a different effective funciton of k_perp
bias = 0.13
zref = 2.3
beta = 1.6
Fmean = 0.8
#pixel = 69 # km/s   for DESI at fixed z, we have fixed km/s ?
pixel = 0.2 # Mpc/h pixels to which we apply FGPA
N=1000   # changing from 100 to 1000 changes sig_f by 2E-6 only
Mpc2kms = 100
# taken directly forn the produced boxes
sig_l = constant.sigma_l
sig_eta = constant.sigma_eta
cov = constant.mean_delta_l_eta - constant.mean_delta_l * constant.mean_eta #cov(delta_L,eta)
print("cov(sig_l,eta) and rho",cov,cov/sig_l/sig_eta)

def mysinc(x):
	return np.sinc(x/np.pi)		#np.sinc(x)= sin(pi*x)/(pi*x)

def P_RSD(k_par,k_perp,P, beta):
    k=np.sqrt(k_par*k_par+k_perp*k_perp)
    mu = k_par/np.maximum(k,1E-15)
    return (1+beta*mu**2)**2 * P(k)

def P_RSD_W(k_par,k_perp,P, beta):
	k=np.sqrt(k_par*k_par+k_perp*k_perp)
	mu = k_par/np.maximum(k,1E-15)
	Wg = np.exp(- DX*DX*k_par) 
	Wperp = mysinc(Dperp*k_perp/2)**2   # how much for Dperp or another function ? <==
	Wpar = mysinc(DX*k_par/2)**2
	return Wg * Wpar * Wperp * (1+beta*mu**2)**2 * P(k)

def PF1D(k):    # P_camb^1D -> P_F 
    return bias*bias*growth*growth*Pcamb1D(k)

def PF1Dkms(k):    # Mpc/h -> km/s  useful ?
    return Mpc2kms*PF1D(k*Mpc2kms)
    
def Ps(k):
	return Pcamb1D(k) - Pcamb1DW(k)


kk=np.arange(kmax/dk)*dk
growth = util.fgrowth(zref, constant.omega_M_0)
print("growth=",growth)

#................. compute P_F^1D including RSD contribution

P_camb = powerspectrum.P_0("$SACLAYMOCKS_BASE/etc/PlanckDR12.fits")  # class P_0
#p1d_camb = powerspectrum.P_1D(k, P_camb.P(kk)).P1D(kk)

k_par = np.arange(kmax/dk)*dk
k_par_t = k_par.reshape(len(k_par),-1)
k_perp = np.arange(kmax/dk)*dk
k_perp_t = k_perp.reshape(len(k_perp),-1)
PRSD = P_RSD(k_par_t,k_perp,P_camb.P, beta)
Pcamb1D = powerspectrum.P_1D_RSD(k_par,k_perp,PRSD).P1D   # in MPc/h

plt.plot(kk,PF1D(kk))
plt.grid()
plt.axis(xmin=0,ymin=0)
#plt.show()

# ~ #.................. compute sigma_F in 69 km/s pixels from PF1D in km/s
# ~ L = N*69
# ~ kj = 2*np.pi / L * np.arange(1, N/2)
# ~ var_F = 2*PF1Dkms(kj).sum() / L
# ~ var_F += PF1Dkms(0) / L  # kj=0 term
# ~ var_F += PF1Dkms(np.pi/pixel) / L  # kj=k_nyquist term
# ~ sigma_F = np.sqrt(var_F)
# ~ print("sigma_F=",sigma_F)

#
#.................. compute sigma_F in 0.2 Mpc/h pixels from PF1D in Mpc/h
L = N*pixel
kj = 2*np.pi / L * np.arange(1, N/2)
var_F = 2*PF1D(kj).sum() / L
var_F += PF1D(0) / L  # kj=0 term
var_F += PF1D(np.pi/pixel) / L  # kj=k_nyquist term
sigma_F = np.sqrt(var_F)
print("sigma_F=",sigma_F)

#................... <F>, sig_F  =>  a, b sig_g
Fmean = 0.791119; sigma_F = 0.360234 # for test ln(a), b*sig_g: -5.3199529 5.71883103
ln_a,bsig_g = pdf.fitab(Fmean,sigma_F)
print("ln(a), b*sig_g:",ln_a,bsig_g)

#................... P_s(k) = P_{m}^{1D} -  P_{m,pixel,smearing}^{1D}
PRSDW = P_RSD_W(k_par_t,k_perp,P_camb.P, beta)
Pcamb1DW = powerspectrum.P_1D_RSD(k_par,k_perp,PRSDW).P1D   # in MPc/h

# sig_s
L = N*pixel
kj = 2*np.pi / L * np.arange(1, N/2)
var_s = 2*Ps(kj).sum() / L
var_s += Ps(0) / L  # kj=0 term
var_s += Ps(np.pi/pixel) / L  # kj=k_nyquist term
sigma_s = np.sqrt(var_s)
print("sigma_F=",sigma_s)


# sig_g and b
# sig_g^2 = sig_l^2 + sig_s^2 + c^2 sig_eta^2 + 2 c cov(delta_L,eta) -> b
sig_g2 = sig_l**2+ var_s + beta*beta*sig_eta + 2 * beta * cov
b = bsig_g / np.sqrt(sig_g2)
print(ln_a,b)
