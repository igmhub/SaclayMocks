#!/usr/bin/env python

# c= beta
# b,beta -> P_F^3D -> P_F^1D -> sig_F^2
# <F> sig_F -> a, b*sig_g
# P_s(k) = P_{mat}^{1D} -  P_{mat,pixel,smearing}^{1D} -> sigma_s
# sig_g^2 = sig_l^2 + sig_s^2 + c^2 sig_eta^2 + 2 c cov(delta_L,eta) -> b

# simu1D/src/fitab.cpp
# Fmean(z) + sigF(z) -> a, b
# does not compile !

from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
from SaclayMocks import powerspectrum
from SaclayMocks import FGPApdf as pdf



b = 0.1
beta = 1.6
Fmean = 0.8
pixel = 69 # km/s   for DESI at fixed z, we have fixed km/s ?
N=100   # changing from 100 to 1000 changes sig_f by 2E-6 only
Mpc2kms = 100

def P_RSD(k_par,k_perp,P, beta):
    k=np.sqrt(k_par*k_par+k_perp*k_perp)
    mu = k_par/np.maximum(k,1E-15)
    return (1+beta*mu**2)**2 * P(k)

def PF1D(k):
    return Mpc2kms*b*b*Pcamb1D(k)

kmax = 5    # 20.
dk=0.1 # 0.001
kk=np.arange(kmax/dk)*dk

#................. compute P_F^1D

P_camb = powerspectrum.P_0("/Users/legoff/ana/SaclayMocks/etc/PlanckDR12.fits")
#p1d_camb = powerspectrum.P_1D(k, P_camb.P(kk)).P1D(kk)

k_par = np.arange(kmax/dk)*dk
k_par_t = k_par.reshape(len(k_par),-1)
k_perp = np.arange(kmax/dk)*dk
k_perp_t = k_perp.reshape(len(k_perp),-1)
PRSD = P_RSD(k_par_t,k_perp,P_camb.P, beta)
Pcamb1D = powerspectrum.P_1D_RSD(k_par,k_perp,PRSD).P1D   # in MPc/h

#plt.plot(kk,b*b*Pcamb1D(kk)/PF1D(kk))
plt.plot(kk,PF1D(kk))
plt.show()

#.................. compute sigma_F

L = N*pixel
kj = 2*np.pi / L * np.arange(1, N/2)
# xx = PF1D(kj)
var_F = 2*PF1D(kj).sum() / L
var_F += PF1D(0) / L  # kj=0 term
var_F += PF1D(np.pi/pixel) / L  # kj=k_nyquist term
sigma_F = np.sqrt(var_F)

print("sigam_F=",sigma_F)

# call fitab.cpp and enter "a" and b*sig_g
ln_a,bsig_g = pdf.fitab(Fmean,sigma_F)
print("la(a), b*sig_g:",ln_a,bsig_g)

# P_s(k) 

# sig_s

# sig_g and b
