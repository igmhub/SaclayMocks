import scipy as sp
import scipy.stats as stats
from scipy import interpolate
from scipy import integrate
import numpy as np
import numpy.ma as ma
import healpy as hp
import subprocess
import os
import sys
from matplotlib import pyplot as plt
import fitsio
from fitsio import FITS
from SaclayMocks import constant
import h5py


PI = np.pi

# *************************************************************
#   average m pixels together, leaving away the N%m remaining pixels


def regroup(spectrum, m):
    m = int(m)
    N=len(spectrum)
    p = N / m
#    print N,p*m # prov
    spectrum = spectrum[0:p*m]
    xx = spectrum.reshape(p,m)
    xx = xx.mean(1)
    return xx

# different way of rebinning pixel, for 2D arrays
# shape need to be a divisor or a.shape
def rebin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)

#*************************************************************
def MakeProfileHisto(x,y,bins=50,err=True) :
    ''' returns the mean value of y and its error in several bins in x
    '''
    if len(np.where(np.isnan(x))[0]) > 0:
        print("x is nan here: {}".format(np.where(np.isnan(x))[0]))
        msk = (np.isnan(x) == False) & (np.isnan(y) == False)
        x = x[msk]
        y = y[msk]

    means_result = stats.binned_statistic(x, [y, y**2, x], bins=bins, statistic='mean')
    meany, meany2, meanx = means_result.statistic
    if len(np.where(np.isnan(meanx))[0]) > 0:
        msk = (np.isnan(meanx) == False) & (np.isnan(meany) == False) & (np.isnan(meany2) == False)
        meanx = meanx[msk]
        meany = meany[msk]
        meany2 = meany2[msk]
    erry = np.sqrt(np.maximum(meany2 - meany**2 , 0))
    if (err):   # if err = False, just return rms
        means_result = stats.binned_statistic(x, [y], bins=bins, statistic='count')
        N_in_bin = means_result.statistic  # number of x value in each bin
        N_in_bin = N_in_bin[0]  #  N_in_bin.shape = (1, bins)
        try:
            erry /= np.sqrt(N_in_bin)    # when N_in_bin = 0, anyway erry = nan already
        except ValueError:
            erry /= np.sqrt(N_in_bin[msk])
    if len(np.where(np.isnan(meanx))[0]) > 0:
        print("mean x is nan here: {}".format(np.where(np.isnan(meanx))[0]))
        msk = (np.isnan(meanx) == False) & (np.isnan(meany) == False) & (np.isnan(erry) == False)
        meanx = meanx[msk]
        meany = meany[msk]
        erry = erry[msk]

    return meanx,meany,erry


#*************************************************************
def ChunkHalfWidth(dec,angle_max):
    ''' minimum half width in RA of a square chunk (+/- angle_max) located at a given dec '''
    a1 = (dec-angle_max)*PI/180
    a2 = (dec+angle_max)*PI/180
    if (a1*a2 < 0) :
        return angle_max # the chunk covers dec=0
    else :
        cos1 = np.cos(a1)
        cos2 = np.cos(a2)
        return angle_max / np.maximum(cos1,cos2)


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        print("Error: boolean value expected.")
        exit()


def radec2pix(nside, ra, dec, nest=True):
    '''ra dec in degrees'''
    phi = ra*np.pi/180
    theta = np.pi/2 - dec*np.pi/180
    return hp.ang2pix(nside, theta, phi, nest=nest)


def pix2radec(nside, pix, nest=True):
    '''ra dec in degrees'''
    theta, phi = hp.pix2ang(nside, pix, nest=nest)
    ra = phi*180/np.pi
    dec = (np.pi/2 - theta) * 180/np.pi
    return ra, dec


def diffmod(a,b,c) :
    ''' absolute value of difference of a and b mod c, e.g. c = 360 deg
    '''
    d = (a-b) %c
    return np.minimum(d,c-d)


def fgrowth(z, Om0, unnormed=False):
    # Assume flat lambda CDM cosmo, with only Om and Ol
    # Comes from cosmolopy.perturbation
    Om = 1 / (1 + (1 - Om0)/(Om0*(1+z)**3))
    Ol = 1 - Om
    a = 1 / (1+z)
    if unnormed:
        norm = 1.0
    else:
        norm = 1.0 / fgrowth(0.0, Om0, unnormed=True)
    return (norm * (5./2.) * a * Om /
            (Om**(4./7.) - Ol + (1. + Om/2.) * (1. + Ol/70.)))


def primes(n):
    ''' decomposition in prime factors
    '''
    # more efficient: https://pypi.python.org/pypi/primefac
    primfac = []
    d = 2
    while d*d <= n:
        while (n % d) == 0:
            primfac.append(d)
            n //= d
        d += 1
    if n > 1:
       primfac.append(n)
    return primfac


def P1D_datafit(k, z):
    # Return the P1Dk, fitted on data (Given by Christophe Yeche)
    # fit is ok for 2.1 < z < 4.5, below is just interpolation
    # k is in (km/s)^-1
    # k can be an array

    z0 = 2.1
    dz = 0.2
    i = (z - z0)/dz
    if ((i < 0) | (i >= 12)):
        print("P1Dk: z = {} out of range\n".format(z))
        sys.exit()
    P = np.zeros_like(k)
    msk = k > 0
    P[msk] = -3.4 + 0.175*i - 0.35*np.log(k[msk]) - 0.075*np.log(k[msk])**2
    P[msk] = np.pi * np.exp(P[msk]) / k[msk]
    return P


def read_P1D(redshift):
    """
    Read Pk file from Nathalie. Si III oscillations have been removed.
    Format is z, k, pk, pkerr, 0, 0, 0
    k in km/s and pk and pkerr in (km/s)**-1
    """
    filename = os.path.expandvars("$SACLAYMOCKS_BASE/etc/pk_fft35bins_noSi.out")
    data = np.loadtxt(filename)
    z = np.round(data[:, 0], 3)
    msk = np.where(z == np.round(redshift, 3))[0]
    if len(msk) == 0:
        # z from 2.2 to 4.4
        print("ERROR -- You entered a wrong redshift: {}. Here is the list of redshift : {}".format(redshift, np.unique(z)))
        sys.exit()
    k = data[:, 1][msk]
    Pk = data[:, 2][msk]
    Pkerr = data[:, 3][msk]
    return k, Pk, Pkerr


def read_P1D_fit(redshift):
    """
    Read Pk fits file, fitted on eBOSS data
    kPk column is k*Pk/pi
    output is k in km/s and pk in s/km
    """
    if redshift < 3.9:
        filename = os.path.expandvars("$SACLAYMOCKS_BASE/etc/models_eBOSS_lowz.fits")
    else:
        filename = os.path.expandvars("$SACLAYMOCKS_BASE/etc/models_eBOSS_highz.fits")
    data = fitsio.read(filename, ext=1)
    z = np.round(data['z'], 3)
    msk = np.where(z == np.round(redshift, 3))[0]
    if len(msk) == 0:
        # z from 2.2 to 4.4
        print("ERROR -- You entered a wrong redshift: {}. Here is the list of redshift : {}".format(redshift, np.unique(z)))
        sys.exit()
    k = data['k'][msk]
    kPk = data['kPk'][msk]
    Pk = kPk / k * np.pi
    return k, Pk


def computechi2(mod, data, dataerr):
    return (((mod-data) / dataerr)**2).sum()


def hist2D(x, y, z=None, bins=100, xlabel=None, ylabel=None, zlabel=None, title=None, vmin=None, vmax=None):
    f, ax = plt.subplots()
    if z is not None:
        res = stats.binned_statistic_2d(x, y, z, bins=bins)
        X, Y = np.meshgrid(res.x_edge, res.y_edge)
        cm = ax.pcolormesh(X, Y, np.nan_to_num(res.statistic.T), vmin=vmin, vmax=vmax)
        cbar = f.colorbar(cm)
    else:
        h = ax.hist2d(x, y, bins=bins)
        cbar = f.colorbar(h[3], ax=ax)
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=20)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=20)
    if title is not None:
        ax.set_title(title, fontsize=20)
    if zlabel is not None:
        cbar.set_label(zlabel, fontsize=20, rotation=90)
    f.tight_layout()
    plt.show()


def read_transmission(inDir, ext, field=None, debug=False):
    x = []
    pix = []
    for d in os.listdir(inDir):
        if os.path.isdir(inDir+d):
            for e in os.listdir(inDir+d):
                if os.path.isdir(inDir+d+'/'+e):
                    for f in os.listdir(inDir+d+'/'+e):
                        if f[:12] == 'transmission':
                            print("Reading : {}".format(inDir+d+'/'+e+'/'+f))
                            fits = FITS(inDir+d+'/'+e+'/'+f)
                            if field is None:
                                try:
                                    y = fits[ext].read()
                                    x.append(y)
                                except:
                                    print("Cannot read {}".format(f))
                                    pix.append(e)
                            else:
                                x = np.concatenate((x, fits[ext][field].read()))
                            fits.close()
    if debug:
        return x,pix
    else:
        return x


class InterpFitsTable():
    '''Read table from a fits file, and compute the interpolate function'''
    def __init__(self, inDir, field1, field2):
        fits = FITS(inDir)
        data = fits[1].read()
        z = data[field1]
        y = data[field2]
        fits.close()
        self.zmin = z.min()
        self.zmax = z.max()
        self.f_of_z = interpolate.interp1d(z, y)
    def interp(self, redshift):
        if np.array(redshift).min() < self.zmin or np.array(redshift).max() > self.zmax:
            raise ValueError("ERROR: Redshift {} is out of range !\nPlease enter a redshift between {} and {}".format(redshift, self.zmin, self.zmax))
        return self.f_of_z(redshift)


def array_interp(xa, ya, xb, yb):
    '''Interpolate ya and yb and return x_interp, ya_interp, yb_interp'''
    xmin = np.max([xa.min(), xb.min()])
    xmax = np.min([xa.max(), xb.max()])
    dx = np.max([xa[1]-xa[0], xb[1]-xb[0]])
    x_interp = np.arange(xmin, xmax, dx)
    f = interpolate.interp1d(xa, ya)
    ya_interp = f(x_interp)
    f = interpolate.interp1d(xb, yb)
    yb_interp = f(x_interp)
    return x_interp, ya_interp, yb_interp


def iterfiles(root, prefix):
    '''
    Returns iterator over files starting with 'prefix' found under 'root' dir
    from desispec.io.util
    '''
    for dirpath, dirnames, filenames in os.walk(root, followlinks=True):
        for filename in filenames:
            if filename.startswith(prefix):
                yield os.path.join(dirpath, filename)


def h_of_z(redshift, Om, H0):
    return H0 * np.sqrt(Om*(1+redshift)**3 + (1-Om))


def zeff(filename, rmin=80., rmax=120.):
    data = fitsio.read(filename, ext=1)
    z = data['Z']
    we = np.diagonal(data['CO'])
    rr = np.sqrt(data['RT']**2+data['RP']**2)
    msk = (rr > rmin) & (rr < rmax)
    return np.sum(z[msk]*we[msk]) / np.sum(we[msk])


class desi_footprint():
    def __init__(self, filename="$SACLAYMOCKS_BASE/etc/desi-healpix-weights.fits"):
        '''
        return a mask to select only desi footprint
        from quickquasars script (desisim)
        '''
        desi_nside =256  # same resolution as original map (cf quickquasars)
        self.desi_nside = desi_nside
        if filename is None:
            filename = os.path.expandvars("$SACLAYMOCKS_BASE/etc/desi-healpix-weights.fits")
        pixmap = fitsio.read(filename, ext=0)
        npix = len(pixmap)
        truenside = hp.npix2nside(npix)
        if truenside < desi_nside:
            print("Warning downsampling is fuzzy...Passed nside={}, but file {} is stored at nside={}".format(truenside, filename, desi_nside))
        healpix_weight = hp.pixelfunc.ud_grade(pixmap, desi_nside, order_in='NESTED', order_out='NESTED')
        self.healpix_weight = healpix_weight
    def selection(self, ra, dec):
        healpix = radec2pix(self.desi_nside, ra, dec)
        return np.where(self.healpix_weight[healpix] > 0.99)[0]


def sigma_p1d(p1d, pixel=0.2, N=10000):
    '''
    p1d is an interpolated function of P1Dmissing (spline basically)
    '''
    L = N*pixel
    kj = 2*np.pi / L * np.arange(1, N/2)
    sigma_2 = 2*p1d(kj).sum() / L
    sigma_2 += p1d(0) / L  # kj=0 term
    sigma_2 += p1d(np.pi/pixel) / L  # kj=k_nyquist term
    return np.sqrt(sigma_2)


def taubar_over_a(sigma, growth, bb=1.58):
    tau = np.exp(0.5*bb**2*growth**2*sigma**2)
    return tau


def fgpa(delta, eta_par, growthf, a, b, c):  #, redshift, taubar_over_a=None):
    # # method 1:
    # tau_over_a = np.exp(b*growthf*delta) - c*(1+redshift)*taubar_over_a*eta_par
    # method 2:
    # tau_over_a = np.exp(b*(growthf*delta - c*(1+redshift)*eta_par))
    # FGPA:
    tau_over_a = np.exp(b*growthf*(delta + c*eta_par))
    flux = np.exp(-a*tau_over_a)
    # flux = np.exp(-1*tau_over_a)  # prov

    # flux = 1 + a*b*growthf*(delta + c*eta_par)

    return flux


def dof(rmin, rmax, binsize=4):
    return 1/(4.*binsize**2) * np.pi * (rmax**2 - rmin**2)


def pol(x, p):
    deg = len(p)-1
    res = 0
    for i, pp in enumerate(p):
        res += pp*x**(deg-i)
    return res


def read_p1dmiss(z, k, filename=None):
    '''
    This function reads the coefficients of the polynomial fit on P1Dmissing
    and uses them to interpolate it to the redshift z
    The fit was done on P1D/(1+z)**3.8
    '''
    if filename is None:
        filename = os.path.expandvars("$SACLAYMOCKS_BASE/etc/pkmiss_fit.fits")
    data = fitsio.read(filename, ext=1)
    f = interpolate.interp1d(data['k'], data['p'], axis=0)
    pkinterp = np.array([pol(z, f(kk)) for kk in k])
    pkinterp *= (1+z)**3.8
    return pkinterp


class InterpP1Dmissing():
    '''Read P1D from fits file, and compute the interpolate function'''
    def __init__(self, infile):
        data = fitsio.read(infile, ext=1)
        z = data['z']
        k = data['k']
        pk = data['Pk']
        self.z = z
        self.k = k
        self.pk = pk
        self.zmin = z.min()
        self.zmax = z.max()
        self.pk_interp = {}

    def interp(self, redshift):
        if redshift < self.zmin or redshift > self.zmax:
            raise ValueError("ERROR: Redshift {} is out of range !\nPlease enter a redshift between {} and {}".format(redshift, self.zmin, self.zmax))
        iz = np.argsort(np.abs(self.z - redshift))[0]
        z = self.z[iz]
        self.pk_interp[str(z)] = interpolate.interp1d(self.k[iz], self.pk[iz])

    def __call__(self, redshift, k):
        iz = np.argsort(np.abs(self.z - redshift))[0]
        z = self.z[iz]
        if str(z) not in self.pk_interp.keys():
            self.interp(redshift)
        return self.pk_interp[str(z)](k)


def find_A_in_B(A,B):
    '''This function looks for the A elements the B array.
    it assumes that each element of A has one matching element in B
    (several A elements can have the same B corresponding element)
    Return the array containing the indices such that
    B[idx] = A
    Function from Charles Antoine Claveau - CEA Saclay
    '''
    B_sort_ind = np.argsort(B)
    B_max = B.max()
    B_extended = np.ones(B_max+1, dtype = np.uint64) * (B_max+1)
    B_extended[B[B_sort_ind]] = B_sort_ind
    A_B_matched_ind = B_extended[A]
    return A_B_matched_ind


def convert1DTo2D(array1D,nbX=None,nbY=None):
    '''
        convert a 1D array to a 2D array
    '''
    if nbX is None and nbY is None:
        nbX = int(np.sqrt(len(array1D)))
        nbY = nbX
    array2D = np.zeros(shape=(nbX,nbY))
    for k in range(array1D.size):
        i = k//nbY
        j = k%nbY
        array2D[i][j] = array1D[k]
    return array2D


def bias_qso(redshift):
    '''
    This function return the bias of QSO for a given redshift
    The parametrisation comes from P. Laurent et al (2017)
    '''
    return 3.7 * ((1+redshift)/(1+2.33))**1.7


def qso_a_of_z(redshift, z_qso_bias):
    return bias_qso(redshift)*(1+z_qso_bias)/(bias_qso(z_qso_bias)*(1+redshift))


def qso_lognormal_coef(filename='$SACLAYMOCKS_BASE/etc/qso_lognormal_coef.txt'):
    '''
    This function returns the interpolated coefficient for computing the
    interpolation between the 3 different lognormal fields
    '''
    filename = os.path.expandvars(filename)
    data = np.loadtxt(filename)
    z = data[:,0]
    coef = data[:,1]
    # Set coef=0 for redshift values above z.max() in txt file
    z = np.append(z, 10)
    coef = np.append(coef, 0)
    # Set coef=1 for redshift values below z.min() in txt file
    z = np.append(0, z)
    coef = np.append(1, coef)
    f = interpolate.interp1d(z, coef)
    return f


def extract_h5file(fname):
    '''
    This function is taken from picca
    https://github.com/igmhub/picca/blob/master/py/picca/fitter2/effective-bins.py
    '''
    f = h5py.File(os.path.expandvars(fname),'r')

    free_p = [ el.decode('UTF-8') for el in f['best fit'].attrs['list of free pars'] ]
    fixed_p = [ el.decode('UTF-8') for el in f['best fit'].attrs['list of fixed pars'] ]
    pars = { el:f['best fit'].attrs[el][0] for el in free_p }
    err_pars = { el:f['best fit'].attrs[el][1] for el in free_p }
    pars.update({ el:f['best fit'].attrs[el][0] for el in fixed_p })
    err_pars.update({ el:0. for el in fixed_p })

    f.close()

    return free_p, fixed_p, pars, err_pars

class cosmo:
    '''
    From picca.constant.py
    https://github.com/igmhub/picca/blob/master/py/picca/constants.py
    '''
    def __init__(self,Om,Ok=0.,Or=0.,wl=-1.,H0=100.):

        ### Ignore evolution of neutrinos from matter to radiation
        ### H0 in km/s/Mpc
        c = constant.c
        Ol = 1.-Ok-Om-Or

        nbins = 10000
        zmax = 10.
        dz = zmax/nbins
        z = sp.arange(nbins)*dz
        hubble = H0*sp.sqrt( Ol*(1.+z)**(3.*(1.+wl)) + Ok*(1.+z)**2 + Om*(1.+z)**3 + Or*(1.+z)**4 )

        chi = sp.zeros(nbins)
        for i in range(1,nbins):
            chi[i] = chi[i-1]+c*(1./hubble[i-1]+1./hubble[i])/2.*dz

        self.r_comoving = interpolate.interp1d(z,chi)

        ### dm here is the comoving angular diameter distance
        if Ok==0.:
            dm = chi
        elif Ok<0.:
            dm = sp.sin(H0*sp.sqrt(-Ok)/c*chi)/(H0*sp.sqrt(-Ok)/c)
        elif Ok>0.:
            dm = sp.sinh(H0*sp.sqrt(Ok)/c*chi)/(H0*sp.sqrt(Ok)/c)

        self.hubble = interpolate.interp1d(z,hubble)
        self.r_2_z = interpolate.interp1d(chi,z)

        ### D_H
        self.dist_hubble = interpolate.interp1d(z,c/hubble)
        ### D_M
        self.dm = interpolate.interp1d(z,dm)
        ### D_V
        y = sp.power(z*self.dm(z)**2*self.dist_hubble(z),1./3.)
        self.dist_v = interpolate.interp1d(z,y)
