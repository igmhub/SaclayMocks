# get P(k) by reading a file and interpolating
# P(k) => xi(r) and xi(r) => P(k)
import os
import numpy as np
import scipy as sp
from scipy import interpolate, integrate
import matplotlib.pyplot as plt # prov
import pyfftw.interfaces.numpy_fft as fft
from fitsio import FITS
from pkg_resources import resource_filename

from SaclayMocks import util
from SaclayMocks import constant


#********************************************************************
class P_1D() :
    '''Computes P^1D = (1/2PI) int_{k_//}^\infty P^3D(k)kdk  '''
    def __init__(self,k,P,kmax=-1):
        if (kmax>0) :
            k=k[np.where(k<=kmax)]
            P=P[np.where(k<=kmax)]
#        print(k.shape,P.shape)
        P1D=np.zeros(len(k))
        for i in np.arange(len(k)) :      #  (1/2PI) int_{k_//}^\infty P(k)kdk
            P1D[i] = np.trapz((P*k)[i:],k[i:]) /2/np.pi
        #dk = k[1:]-k[0:-1]
        #dk = np.append(dk,dk[-1])
        #P1D = np.cumsum(P*k*dk[-1::-1]) # from last to first, in inverse order
        #P1D = P1D[-1::-1]
        #       completely fails ???
        self.pk1DInter=interpolate.InterpolatedUnivariateSpline(k,P1D)

    def P1D(self,k):
        return np.maximum(self.pk1DInter(k),0)
            # to avoid P(0) = -1.38775417635e-17

#********************************************************************
class P_1D_RSD() :
    '''Computes P^1D = (1/2PI) int_0^\infty P^3D(k_par,k_perp)k_perp dk_perp  '''
    def __init__(self,k_par,k_perp,P):  # P(k_par,k_perp) 2D array
        P1D=np.zeros(len(k_par))
                #  (1/2PI) int_0^\infty P(k_//,k_perp)k_perp dk_perp
        for i in np.arange(len(k_par)) :
            P1D[i] = np.trapz(k_perp*P[i,:],k_perp) /2/np.pi
        self.pk1DInter=interpolate.InterpolatedUnivariateSpline(k_par,P1D)

    def P1D(self,k):
        return np.maximum(self.pk1DInter(k),0)
            # to avoid P(0) = -1.38775417635e-17



#********************************************************************
class P_0() :
    '''Read P(k) from a file and interpolate it'''
    def __init__(self,filename=None,G_times_bias=1):
          #     read filename skiping skiprows lines
          #     k in column colk and p in colp
        if filename is None:
            filename = resource_filename('SaclayMocks', '/etc/PlanckDR12.fits')
        fits = FITS(filename)
        data = fits[1].read()
        k = data['K']
        P = data['PK']
        zref = fits[1].read_header()['ZREF']
        fgrowth = util.fgrowth(zref, constant.omega_M_0)
        P /= fgrowth**2  # go to z=0
        fits.close()
        zero = np.arange(1)
        k = np.append(zero,k)
        P = np.append(zero,P)
        P *= G_times_bias**2

        # # Add Gaussian smearing
        # DX = 2.19
        # P *= np.exp(-k**2 * DX**2)

        self.pkInter=interpolate.InterpolatedUnivariateSpline(k,P)
    def P(self,k):
        return np.maximum(self.pkInter(k),0)


#********************************************************************
# could be same class as P_0 with
#       if (logNormal) : k_input,P_input = LogNormalP(k_input,P_input)
class P_ln() :
    '''Read P(k) from a file, compute the lognormal P(k), interpolate it'''
    def __init__(self,filename=None,G_times_bias=1):
          #     read filename skiping skiprows lines
          #     k in column colk and p in colp
        if filename is None:
            filename = resource_filename('SaclayMocks', '/etc/PlanckDR12.fits')
        fits = FITS(filename)
        data = fits[1].read()
        k_input = data['K']
        P_input = data['PK']
        zref = fits[1].read_header()['ZREF']
        fgrowth = util.fgrowth(zref, constant.omega_M_0)
        P_input /= fgrowth**2  # go to z=0
        fits.close()
        zero = np.arange(1)
        k_input = np.append(zero,k_input)
        P_input = np.append(zero,P_input)
        P_input *= G_times_bias**2

        # # Add Gaussian smearingxb
        # DX = 2.19
        # P_input *= np.exp(-k_input**2 * DX**2)

        kln, Pln = LogNormalP(k_input,P_input)
        self.pkInter=interpolate.InterpolatedUnivariateSpline(kln,Pln)
    def P(self,k):
        return np.maximum(self.pkInter(k),0)


#********************************************************************
#def P_1D(k) :
# should be part of class P_0 with an initialization and then interpolation

#********************************************************************
#  from P(k) to xi(r) for uneven spaced k points
#********************************************************************
#  r = 2 PI m / kmax
#  nr = nk / 2
#  r_max = PI nk / kmax
# xi(r) = - 1/(2r PI^2) \im FFT kP(k)
# dr = pi / kmax    r_max = N dr / 2 = N pi/ 2 k_max
#  kmax = 10 -> dr = 0.314 and N >= 1000 advisable
def xi_from_pk(k,pk,nk=32*1024,direct=True):
    if (k[0] != 0) :
        k = np.append([0],k)
        pk = np.append([0],pk)
    pkInter=interpolate.InterpolatedUnivariateSpline(k,pk) #,kind='cubic')
    kmax=np.max(k)
    if direct:
        if (kmax < 10) :
            print("**** WARNING kmax=",kmax," while kmax >= 10 advised ****")
        rmax = np.pi * nk / kmax
        if (rmax < 1000) :
            print("**** WARNING rmax=",rmax," while rmax >= 1000 advised ****")
    kIn=np.linspace(0,kmax,nk)
    pkIn=pkInter(kIn)
    r=2.*np.pi*np.arange(nk)/kmax
    pkk=kIn*pkIn
    r[0]=1E-10
    cric=-np.imag(np.fft.fft(pkk)/nk)/r/2./np.pi**2*kmax    # r[0]=0 => RuntimeWarning
    r[0]=0
    pkInter=interpolate.InterpolatedUnivariateSpline(k,pk*k*k) #,kind='cubic')
    cric[0]=pkInter.integral(0,kmax) /2/np.pi**2    #   (1/2 PI^2) int P(k) k^2 dk
    r=r[0:nk/2]
    cric=cric[0:nk/2]
    return r,cric

#********************************************************************
#  from xi(r) to P(k) for uneven spaced k points
#********************************************************************
#  xi(r) = (2 pi)^{-3} \int exp(ikr) Pk) dk
#  P(k) = \int exp(ikr) Pk) dk
#  so for xi -> P same function as for P -> xi but multiply by (2 pi)^3
def pk_from_xi(r,xi,nr=32*1024,direct=False):
    k, Pk = xi_from_pk(r,xi,nr)
    Pk *= (2 * np.pi)**3
    return k, Pk


#********************************************************************
# P(k) => xi(r) => cln(r) = log [1 + \xi(r) ] => Pln(k)
def LogNormalP(k,P,nk=1024*1024):
    r , xi = xi_from_pk(k,P,nk=nk)
    cln = np.log(1+xi)
    nr = nk/2
    kln , Pln = pk_from_xi(r,cln,nr=nr)
    Pln = np.maximum(Pln,0)
    return kln, Pln


#*************************************************************
def P1D_1spectrum(delta,DX) :
    ''' compute 1D power spectrum of one spectrum, to be averaged over many spectra '''
    z = fft.rfft(delta)
    x = z.real**2 +z.imag**2
    return x * DX / len(delta)   # (LZ/NZ/NZ) |FFT|^2


#*************************************************************
class ComputeP1D() :

    def __init__(self,DeltaR):
        # self.Pkvec = np.array([], dtype=np.float32)
        # self.kvec = np.array([], dtype=np.float32)
        self.kvec = []
        self.Pkvec = []
        self.DeltaR = DeltaR
        self.k_ny = np.pi / DeltaR
    def add_spectrum(self,delta):
        # nz = 16
        # while (nz < len(delta)) : nz *= 2
        # mydelta = delta[0:nz/2]     # to avoid losing time in FFT
        # freq = np.fft.rfftfreq(len(mydelta)) * 2 * self.k_ny
        # cut = np.where(freq<2)
        # self.kvec = np.append(self.kvec , freq[cut] )
        # if len(np.where(np.isnan(self.kvec))[0]) > 0:
        #     print(np.where(np.isnan(self.kvec))[0])
        # self.Pkvec = np.append(self.Pkvec , P1D_1spectrum(mydelta,self.DeltaR)[cut])
        k = np.fft.rfftfreq(len(delta)) * 2 * self.k_ny
        self.kvec.append(k)
        self.Pkvec.append(P1D_1spectrum(delta, self.DeltaR))
    def P1D(self,bins=50):
        # return util.MakeProfileHisto(self.kvec,self.Pkvec,bins=bins) # k,P1D,errP1D
        kvec = np.concatenate(self.kvec)
        Pkvec = np.concatenate(self.Pkvec)
        return util.MakeProfileHisto(kvec, Pkvec, bins=bins)  #k, P1D, P1Derr


#********************************************************************
def xi_predicted(xi_g, f, dfbar, delta_s=False, sigma=np.inf):
    # This function returns the predicted xi after applying f to a GRF g
    # g should be normalized such that sigma_g=1, which implies |xi_g| <1
    # see eq 2.6 in Font et al. 2012
    # the numerical integral is not reliable for 0.005 < x_g < 1
    # delta_s = True => add small scales in prediction
    if delta_s == True:
        def integrand(dl1, dl2, ds1, ds2):
            num = np.exp(-(dl1**2 + dl2**2 - 2*dl1*dl2*xi_g)/(2*(1-xi_g**2)))
            denum = 2*np.pi*np.sqrt(1 - xi_g**2)
            return ((f(dl1, ds1)/dfbar-1)*(f(dl2, ds2)/dfbar-1)
                    * num/denum * np.exp(-0.5*(ds1**2+ds2**2))/(2*np.pi))
        return integrate.nquad(integrand, [[-sigma, sigma], [-sigma, sigma], [-sigma, sigma], [-sigma, sigma]])
    else:
        if ( (abs(xi_g) > 1) | (xi_g == -1) ) :
            print("******** error in xi_predicted, xi_g=", xi_g)
            exit(0)
        if ( (xi_g > 0.995) & (xi_g < 1) ) :
            print("**** Warning, computation of xi_predicted not reliable for 0.995 < xi_g < 1 *****")
            print("should rather interpolate linearly bewteen 0.995 and 1.")
        if (xi_g==1) :
            def integrand(dl):
                num = np.exp(-dl**2/2)
                denum = np.sqrt(2*np.pi)
                return (f(dl)/dfbar-1)**2 * num/denum
            return integrate.quad(integrand, -np.inf, np.inf)
        def integrand(dl1, dl2):
            num = np.exp(-(dl1**2 + dl2**2 - 2*dl1*dl2*xi_g)/(2*(1-xi_g**2)))
            denum = 2*np.pi*np.sqrt(1 - xi_g**2)
            return (f(dl1)/dfbar-1)*(f(dl2)/dfbar-1) * num/denum
        return integrate.dblquad(integrand, -np.inf, np.inf, lambda x: -np.inf, lambda x: np.inf)

#*************************************************************
class xi_prediction() :
# given a FGPA transformation F = exp[ -a exp (b G g) ]
# and a GRF g = delta + c eta, characterized by sigma_g and c
# computes xi_g and the resulting xi_F
# see Font Ribera at al. (2012) eq 2.7
# DX is the width of the Gaussian smearing applied to the GRF field

    def F(self,delta):  # FGPA function
        return np.exp(-self.a * np.exp (self.b * delta) )

    def F_pdf(self,F):   # returns F x pdf(F)
#	F=exp[ -a exp(b g)] = exp[-exp(bg+ln(a))]
#	where g = gaus(0,1) and bg+ln(a) = Gaus(ln(a),b)
#	pdf(F) = 1/(F tau \sqrt{2PI} b) exp [-(ln(tau)-ln(a))^2/(2b^2)]
#       with tau = -ln(F)
        b=self.b
        a=self.a
        if (b==0):
            print("b=0 in F_pdf\n")
            exit(0)
        if ((F<=0) | (F>=1)):
            print("F=",F, "in F_pdf\n")
            exit(0)
        tau = - np.log(F)
        xx = np.log(tau/a)
        xx = np.exp(- xx*xx/2./b/b)
        # xx /= F * tau * np.sqrt(2 * np.pi) * b  # pdf of F
        xx /= tau * np.sqrt(2 * np.pi) * b      # F x pdf of F
        return xx

    def pdf_integrale(self,Fmax) :
        # compute int_0^Fmax pdf(F)dF using analytical primitive of pdf
        b=self.b
        a=self.a
        u = np.log(a) -np.log(-np.log(Fmax))
        xx = 0.5 * ( 1 + sp.special.erf( u/b/np.sqrt(2) ) )
        print(1-Fmax,a,b,1-xx)
        return xx

    def ComputeFmean(self):     # checked wrt simu1D.cpp
        # F pdf diverges in 1, making numerical integration unstable
        # compute numerical integral [0,1-epsilon]
        # approximate F ~ 1 over [1-epsilon,1]
        # analytically compute integral pdf over [1-epsilon,1]
        epsilon = 1E-8
        Fmean , error = sp.integrate.quad(self.F_pdf,0,1-epsilon,limit=50)
        #print(epsilon,a,b,Fmean)
        Fmean += 1- self.pdf_integrale(1-epsilon) # \int_{1-epsilon}^1
        return Fmean

    def __init__(self,a,b,G,sigma_g,c,DX=2.19):
        self.a = a
        self.b = b * G * sigma_g
        self.sigma_g = sigma_g
        self.c = c
        self.DX = DX
            #...........................       P(k)
        k = np.linspace(0,10,10000)
        #Pcamb = powerspectrum.P_0()
        Pcamb = P_0()
        P = Pcamb.P(k)
        W = np.exp(-k*k*DX*DX/2)
        P *= W*W
            #...........................       xi_g(r)
        r,xi = xi_from_pk(k,P)
        rmax=300
        cut=np.where(r<rmax)
        self.r=r[cut]
        self.xi=xi[cut]
        self.xi_Ham = xi_Hamilton(r,xi,rmax)
            #...........................        xi_g -> xi_F
        Fmean = self.ComputeFmean()
        xx=np.linspace(-0.99,0,41)  # make sure to have zero
        yy=np.linspace(0,1,41)
        xi_g=np.hstack((xx,yy[1:])) # do not have zero twice
        xi_F = np.zeros(len(xi_g))
        for i in range(len(xi_g)):
            xi_F[i], err = xi_predicted(xi_g[i], self.F, Fmean)
            #print(i, xi_g[i], xi_F[i])
        self.xig2xiF=interpolate.InterpolatedUnivariateSpline(xi_g,xi_F)
        self.xi_g_array=xi_g
        self.xi_F_array=xi_F

    def xig_xiF(self):      # returns xi_g and xi_F array
        return self.xi_g_array,self.xi_F_array

    def xi_g(self,r,mu):
        return self.xi_Ham.xi(self.c,r,mu)


#................. main function, returns the predicted xi_F
    def xi_F(self,r,mu):
        xi_g = self.xi_Ham.xi(self.c,r,mu) / self.sigma_g**2
        return self.xig2xiF( xi_g )

#.... BEWARE:  the next three functions return xig2xiF applied
# to the multipoles of xi_g. But since F(g) is non linear, these
# are not exactly the multipoles of xi_F.
    def xi0_F(self,r):
        xi_g = self.xi_Ham.xi0(self.c,r) / self.sigma_g**2
        return self.xig2xiF( xi_g )

    def xi2_F(self,r):
        xi_g = self.xi_Ham.xi2(self.c,r) / self.sigma_g**2
        return self.xig2xiF( xi_g )

    def xi4_F(self,r):
        xi_g = self.xi_Ham.xi4(self.c,r) / self.sigma_g**2
        return self.xig2xiF( xi_g )


#********************************************************************
class xi_Hamilton() :
    '''Computes xi^bar, xi^bar-bar, xi_0, xi_2 xi_4 and xi(r,mu) Hamilton 1992 '''
    def __init__(self,r,xi,rmax=-1):
        if (rmax>0) :
            r=r[np.where(r<=rmax)]
            xi=xi[np.where(r<=rmax)]
        xibar = np.zeros(len(r))
        xibarbar = np.zeros(len(r))
        '''Computes xi^bar (r) = (1/3r^3) =\int_0^r s^2 xi(s) ds
        and  xi^bar-bar (r) = (1/5r^5) =\int_0^r s^4 xi(s) ds   '''
        for i in np.arange(len(r)) :
            if (r[i] == 0 ) :
                xibar[i]=xi[i]
                xibarbar[i]=xi[i]
                continue
            xibar[i] = np.trapz((r*r*xi)[:i],r[:i]) * 3 / r[i]**3
            xibarbar[i] = np.trapz((r**4 * xi)[:i],r[:i]) * 5 / r[i]**5
        self.xibarInter=interpolate.InterpolatedUnivariateSpline(r,xibar)
        self.xibarbarInter=interpolate.InterpolatedUnivariateSpline(r,xibarbar)
        self.xiInter=interpolate.InterpolatedUnivariateSpline(r,xi)

    def xicamb(self,r):
        return self.xiInter(r)
    def xibar(self,r):
        return self.xibarInter(r)
    def xibarbar(self,r):
        return self.xibarbarInter(r)
    def xi0(self,f,r):
        return (1+2*f/3.+f*f/5.)*self.xicamb(r)
    def xi2(self,f,r):
        return (4*f/3.+4*f*f/7.)*(self.xicamb(r)-self.xibar(r))
    def xi4(self,f,r):
        return (8*f*f/35.)*(self.xicamb(r)+2.5*self.xibar(r)-3.5*self.xibarbar(r))
    def xi(self,f,r,mu):
        return self.xi0(f,r) + (3*mu*mu-1)*self.xi2(f,r)/2 + (35*mu**4-30*mu*mu+3)*self.xi4(f,r)/8
#   for XCF see Mountrichas et al. 2009 Y.P. Jing priv comm
    def xcf0(self,f1,f2,r):
        return (1+(f1+f2)/3.+f1*f2/5.)*self.xicamb(r)
    def xcf2(self,f1,f2,r):
        return (2*(f1+f2)/3.+4*f1*f2/7.)*(self.xicamb(r)-self.xibar(r))
    def xcf4(self,f1,f2,r):
        return (8*f1*f2/35.)*(self.xicamb(r)+2.5*self.xibar(r)-3.5*self.xibarbar(r))
    def xcf(self,f1,f2,r,mu):
        return self.xcf0(f1,f2,r) + (3*mu*mu-1)*self.xcf2(f1,f2,r)/2 + (35*mu**4-30*mu*mu+3)*self.xcf4(f1,f2,r)/8



#********************************************************************
def IntegratePKaiser(P,beta,k_max):
    ''' compute 3D integral of P_Kaiser up to k_max
    i.e. 4*pi (1+2*beta/3+beta^2/5) int_0^kmax k^2P(k)
    was aimed at computing sigma^2 but does not work
    '''
    def k2P(k):
        return k*k*P(k)
    intk2P , error = integrate.quad(k2P,0,k_max,limit=100)
    return 4*np.pi*(1+2*beta/3.+beta*beta/5.) * intk2P



#********************************************************************
#********************************************************************
#********************************************************************
#********************************************************************
#********************************************************************

'''
#------------------------------------------------------------
#  from P(k) to xi(r) for uneven spaced k points (version JCH)
#------------------------------------------------------------
def xi_from_pk_JCH(k,pk):
   pkInter=interpolate.InterpolatedUnivariateSpline(k,pk) #,kind='cubic')
   nk=32768
   kmax=np.max(k)
   kmin=np.min(k)
   kIn=np.linspace(kmin,kmax,nk)
   pkIn=pkInter(kIn)
   kIn[0]=0.
   pkIn[0]=0.
   r=2.*np.pi/kmax*np.arange(nk)
   pkk=kIn*pkIn
   cric=-np.imag(np.fft.fft(pkk)/nk)/r/2./np.pi**2*kmax
   cric[0]=0
   return r,cric


def pk_from_xi_old(k,pk):
   pkInter=interpolate.InterpolatedUnivariateSpline(k,pk) #,kind='cubic')
   nk=32768
   kmax=np.max(k)
   kmin=np.min(k)
   kIn=np.linspace(kmin,kmax,nk)
   pkIn=pkInter(kIn)
   kIn[0]=0.
   pkIn[0]=0.
   r=2.*np.pi/kmax*np.arange(nk)
   pkk=kIn*pkIn
   cric=-np.imag(np.fft.fft(pkk)/nk)/r/2./np.pi**2*kmax
   cric[0]=0
   return r,cric


def xi_from_pk_old(k,pk):
    pkInter=interpolate.InterpolatedUnivariateSpline(k,pk) #,kind='cubic')

    nk=32768
    #   nk=1024
    kmax=np.max(k)
    kmin=np.min(k)
    print(kmin,kmax)
    kIn=np.linspace(kmin,kmax,nk)
    pkIn=pkInter(kIn)
#    plt.plot(kIn,pkIn) # prov
#    plt.show() # prov
    kIn[0]=0.
    pkIn[0]=0.
    r=2.*np.pi/kmax*np.arange(nk)
    pkk=kIn*pkIn
    cric=-np.imag(np.fft.fft(pkk)/nk) * kmax /r/2./np.pi**2
    print(cric)
    cric[0]=0  #   should be (1/2 PI^2) int P(k) k^2 dk
    return r,cric

'''
