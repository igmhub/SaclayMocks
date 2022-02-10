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
from SaclayMocks import constant
import h5py
try:
    import picca.wedgize
    use_picca = True
except:
    print("/!\ Unable to import picca.wedgize !")
    use_picca = False


PI = np.pi

# *************************************************************
#   average m pixels together, leaving away the N%m remaining pixels


def regroup(spectrum, m):
    m = int(m)
    N=len(spectrum)
    p = N / m
#    print(N,p*m) # prov
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
    ''' input: array of data x and y, bins is the number of bins of the profile histogram,  
       bins should not be too large so that all (most) bins contain at least a data point
       returns the mean value of y and its error in bins in x
       if err == False returns rms, else returns rms/sqrt(N)
    '''
    if (len(x)!=len(y)) :
        print("len(x)!=len(y)")
        return
    if len(np.where(np.isnan(x))[0]) > 0:  # this should not happen
       print("x is nan here: {}".format(np.where(np.isnan(x))[0]))
       msk = (np.isnan(x) == False) & (np.isnan(y) == False)
       x = x[msk]
       y = y[msk]

    means_result = stats.binned_statistic(x, [y, y**2, x], bins=bins, statistic='mean')
    meany, meany2, meanx = means_result.statistic
    if len(np.where(np.isnan(meanx))[0]) > 0:   # if there are empty bins, remove them
        msk = (np.isnan(meanx) == False) & (np.isnan(meany) == False) & (np.isnan(meany2) == False)
        meanx = meanx[msk]
        meany = meany[msk]
        meany2 = meany2[msk]
    erry = np.sqrt(np.maximum(meany2 - meany**2 , 0))  #  make sure not to have sqrt(<0) due to rounding error
    if (err):   # if err = False, just return rms
        means_result = stats.binned_statistic(x, [y], bins=bins, statistic='count')
        N_in_bin = means_result.statistic  # number of x value in each bin
        N_in_bin = N_in_bin[0]  #  N_in_bin.shape = (1, bins)
        erry /= np.sqrt(N_in_bin[msk])  # if N_in_bin[i] == 0 then meanx[i] == nan and msk does not include i

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


def fgrowth(z, Om0=constant.omega_M_0, unnormed=False):
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


def read_P1D(redshift, filename="$SACLAYMOCKS_BASE/etc/pk_1d_DR12_13bins_noSi.out", mpc=False):
    """
    Read Pk file from Nathalie. Si III oscillations have been removed.
    Format is z, k, pk, pkerr, 0, 0, 0
    k in km/s and pk and pkerr in (km/s)**-1
    """
    data = np.loadtxt(os.path.expandvars(filename))
    z = np.round(data[:, 0], 3)
    msk = np.where(z == np.round(redshift, 3))[0]
    if len(msk) == 0:
        # z from 2.2 to 4.4
        print("ERROR -- You entered a wrong redshift: {}. Here is the list of redshift : {}".format(redshift, np.unique(z)))
        sys.exit()
    k = data[:, 1][msk]
    Pk = data[:, 2][msk]
    Pkerr = data[:, 3][msk]
    if mpc:
        print("Output is in Mpc/h")
        k *= kms2mpc(redshift)
        Pk /= kms2mpc(redshift)
        Pkerr /= kms2mpc(redshift)
    else:
        print("Output is in km/s")
    return k, Pk, Pkerr


def read_P1D_fit(redshift, mpc=False):
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
    if mpc:
        print("Output is in Mpc/h")
        k *= kms2mpc(redshift)
        Pk /= kms2mpc(redshift)
    else:
        print("Output is in km/s")
    return k, Pk


def read_P1D_model(redshift, filename="$SACLAYMOCKS_BASE/etc/P1DmodelPrats.fits", mpc=False, z_corr=True):
    '''
    This function reads the P1D used as model to tune the P1D shape in mocks
    It returns k, pk for a given redshift.
    Units are in km/s by default (mpc=False)
    The z dependency is smoothed by default (z_corr=True)
    '''
    corrV = np.array([1.0037348 , 0.99913817, 0.99267796, 0.99769447, 1.00333914, 1.00691257, 1.00574199, 0.99331247, 0.99156929, 1.00622274])   #  harcoded correction to smooth sig_F(z)
    i = int(round((redshift-1.8)/0.2))
    if (z_corr): 
        cor = corrV[i]**2  # P ~ sig^2
    else : 
        cor =1
    fits = fitsio.FITS(os.path.expandvars(filename))
    z = fits[0].read()
    k = fits[1].read()
    pk = fits[2].read()
    ### select the redshift bin
    # extrapole z=2.2 to z=2.0
    if np.abs(redshift - 2) < 1e-2:
        msk1 = np.abs(z - 2.2) < 1e-2
        msk2 = np.abs(z - 2.4) < 1e-2
        k = k[msk1]
        pk = pk[msk1]**2 / pk[msk2]
    # extrapole z=2.2 to z=1.8
    elif np.abs(redshift - 1.8) < 1e-2:
        msk1 = np.abs(z - 2.2) < 1e-2
        msk2 = np.abs(z - 2.4) < 1e-2
        k = k[msk1]
        pk = pk[msk1]**3 / pk[msk2]**2
    # read P1D(z)
    else:
        msk = np.abs(z - redshift) < 1e-2
        if msk.sum() == 0:
            print("ERROR -- You entered a wrong redshift: {}. Here is the list of redshifts : {}".format(redshift, np.unique(z)))
            sys.exit(1)
        k = k[msk]
        pk = pk[msk]
    if mpc:
        print("Output is in Mpc/h")
        k *= kms2mpc(redshift)
        pk /= kms2mpc(redshift)
    else:
        print("Output is in km/s")
    return k, pk*cor


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


class InterpFitsTable():
    '''Read table from a fits file, and compute the interpolate function'''
    def __init__(self, inDir, field1, field2):
        fits = fitsio.FITS(inDir)
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


def sigma_p1d(redshift=None, filename="$SACLAYMOCKS_BASE/etc/pkmiss_interp.fits.gz", p1dmiss=None, pixel=0.2, N=10000):
    '''
    Return the sigma of delta_s for a given redshift and a given P1Dmissing
    the p1d can be directly given via p1dmiss argument (it must be a function)
    '''
    L = N*pixel
    kj = 2*np.pi / L * np.arange(1, N/2)
    if p1dmiss is None:
        if redshift is None:
            print("Please enter a valid redshift")
            sys.exit(1)
        redshift = np.array(redshift).reshape(-1)
        sigma_s = np.zeros_like(redshift)
        filename = os.path.expandvars(filename)
        p1dmiss = InterpP1Dmissing(filename)
        for i, z in enumerate(redshift):
            var_s = 2*p1dmiss(z, kj).sum() / L
            var_s += p1dmiss(z, 0) / L  # kj=0 term
            var_s += p1dmiss(z, np.pi/pixel) / L  # kj=k_nyquist term
            sigma_s[i] = np.sqrt(var_s)
    else:
        var_s = 2*p1dmiss(kj).sum() / L
        var_s += p1dmiss(0) / L  # kj=0 term
        var_s += p1dmiss(np.pi/pixel) / L  # kj=k_nyquist term
        sigma_s = np.sqrt(var_s)

    return sigma_s


def sigma_g(redshift, pkfile="$SACLAYMOCKS_BASE/etc/pkmiss_interp.fits.gz", paramfile="$SACLAYMOCKS_BASE/etc/params.fits", p1dmiss=None, c=None, pixel=0.2, N=10000):
    '''
    Returns the sigma of g = delta_l + delta_s + c*eta_par field
    '''
    if c is None:
        c_of_z = InterpFitsTable(paramfile, 'z', 'c')
        c = c_of_z.interp(redshift)
    var_g = constant.sigma_l**2 + c*constant.sigma_eta**2
    var_g += sigma_p1d(redshift, pkfile, p1dmiss, pixel, N)**2
    var_g +=c*(constant.mean_delta_l_eta - constant.mean_delta_l*constant.mean_eta)  # covariance between delta_l and eta_par
    sigma_g = np.sqrt(var_g)
    sigma_g *= constant.sigma_g_tuning  # prov
    return sigma_g


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


class InterpP1Dmissing():
    '''Read P1D from fits file, and compute the interpolate function'''
    def __init__(self, infile):
        fits = fitsio.FITS(infile)
        z = fits['z'].read()
        k = fits['k'].read()
        pk = fits['pk'].read()
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
        self.pk_interp[str(z)] = interpolate.interp1d(self.k, self.pk[iz])

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
    This function creates a dictionnary that contains all the parameters
    stored in the h5file produced by picca
    '''
    f = h5py.File(os.path.expandvars(fname),'r')

    free_p = [ el.decode('UTF-8') for el in f['best fit'].attrs['list of free pars'] ]
    fixed_p = [ el.decode('UTF-8') for el in f['best fit'].attrs['list of fixed pars'] ]
    pars = { el:f['best fit'].attrs[el][0] for el in free_p }
    err_pars = { el:f['best fit'].attrs[el][1] for el in free_p }
    pars.update({ el:f['best fit'].attrs[el][0] for el in fixed_p })
    err_pars.update({ el:0. for el in fixed_p })
    pars['zeff'] = f['best fit'].attrs['zeff']
    pars['chi2'] = f['best fit'].attrs['fval']
    cov_pars = { 'cov[{}, {}]'.format(el1, el2):f['best fit'].attrs['cov[{}, {}]'.format(el1, el2)] for el1 in free_p for el2 in free_p }
    if 'bias_eta_LYA' in f['best fit'].attrs.keys():
        if 'cov[bias_eta_LYA, beta_LYA]' not in f['best fit'].attrs.keys():
            cov_pars['cov[beta_LYA, bias_eta_LYA]'] = 0
        pars['bias_LYA'] = f['best fit'].attrs['bias_eta_LYA'][0] * f['best fit'].attrs['growth_rate'][0] / f['best fit'].attrs['beta_LYA'][0]
        pars['beff_LYA'] = pars['bias_LYA'] * np.sqrt(1+2/3*f['best fit'].attrs['beta_LYA'][0]+1/5*f['best fit'].attrs['beta_LYA'][0]**2)
        # if 'cov[bias_eta_LYA, beta_LYA]' in f['best fit'].attrs.keys():
        err_pars['bias_LYA'] = bias_err(f['best fit'].attrs['bias_eta_LYA'][0], f['best fit'].attrs['bias_eta_LYA'][1], f['best fit'].attrs['beta_LYA'][0], f['best fit'].attrs['beta_LYA'][1], cov_pars['cov[beta_LYA, bias_eta_LYA]'])
        err_pars['beff_LYA'] = beff_err(f['best fit'].attrs['bias_eta_LYA'][0], f['best fit'].attrs['bias_eta_LYA'][1], f['best fit'].attrs['beta_LYA'][0], f['best fit'].attrs['beta_LYA'][1], cov_pars['cov[beta_LYA, bias_eta_LYA]'], f=f['best fit'].attrs['growth_rate'][0])
        # else:
        #     err_pars['bias_LYA'] = 0
        #     err_pars['beff_LYA'] = 0
        free_p += ['bias_LYA', 'beff_LYA']
    if 'bias_LYA' in f['best fit'].attrs.keys():
        if 'cov[bias_LYA, beta_LYA]' not in f['best fit'].attrs.keys():
            cov_pars['cov[beta_LYA, bias_LYA]'] = 0
        pars['b_LYA'] = f['best fit'].attrs['bias_LYA'][0] * f['best fit'].attrs['growth_rate'][0] / f['best fit'].attrs['beta_LYA'][0]
        pars['b_eff_LYA'] = pars['b_LYA'] * np.sqrt(1+2/3*f['best fit'].attrs['beta_LYA'][0]+1/5*f['best fit'].attrs['beta_LYA'][0]**2)

        # if 'cov[bias_LYA, beta_LYA]' in f['best fit'].attrs.keys():
        err_pars['b_LYA'] = bias_err(f['best fit'].attrs['bias_LYA'][0], f['best fit'].attrs['bias_LYA'][1], f['best fit'].attrs['beta_LYA'][0], f['best fit'].attrs['beta_LYA'][1], cov_pars['cov[beta_LYA, bias_LYA]'])
        err_pars['b_eff_LYA'] = beff_err(f['best fit'].attrs['bias_LYA'][0], f['best fit'].attrs['bias_LYA'][1], f['best fit'].attrs['beta_LYA'][0], f['best fit'].attrs['beta_LYA'][1], cov_pars['cov[beta_LYA, bias_LYA]'], f=f['best fit'].attrs['growth_rate'][0])
        # else:
        #     err_pars['b_LYA'] = 0
        #     err_pars['b_LYA'] = 0        
        free_p += ['b_LYA', 'b_eff_LYA']
    if 'bias_eta_QSO' in f['best fit'].attrs.keys():
        if 'cov[growth_rate, beta_QSO]' not in f['best fit'].attrs.keys():
            cov_pars['cov[beta_QSO, growth_rate]'] = 0
        pars['bias_QSO'] = f['best fit'].attrs['bias_eta_QSO'][0] * f['best fit'].attrs['growth_rate'][0] / f['best fit'].attrs['beta_QSO'][0]
        err_pars['bias_QSO'] = bias_err(f['best fit'].attrs['growth_rate'][0], f['best fit'].attrs['growth_rate'][1], f['best fit'].attrs['beta_QSO'][0], f['best fit'].attrs['beta_QSO'][1], cov_pars['cov[beta_QSO, growth_rate]'])
        free_p += ['bias_QSO']
    if 'bias_eta_HCD' in f['best fit'].attrs.keys():
        if 'cov[growth_rate, beta_HCD]' not in f['best fit'].attrs.keys():
            cov_pars['cov[beta_HCD, growth_rate]'] = 0
        pars['bias_HCD'] = f['best fit'].attrs['bias_eta_HCD'][0] * f['best fit'].attrs['growth_rate'][0] / f['best fit'].attrs['beta_HCD'][0]
        err_pars['bias_HCD'] = bias_err(f['best fit'].attrs['growth_rate'][0], f['best fit'].attrs['growth_rate'][1], f['best fit'].attrs['beta_HCD'][0], f['best fit'].attrs['beta_HCD'][1], cov_pars['cov[beta_HCD, growth_rate]'])
        free_p += ['bias_HCD']
    f.close()
    return free_p, fixed_p, pars, err_pars, cov_pars


def print_h5file(fname, cor=False):
    '''
    This function print the h5 output file from picca fitter2
    '''
    pars = extract_h5file(fname)
    print("- Free params:")
    print("zeff = {}".format(pars[2]['zeff']))
    print("chi2 = {}".format(pars[2]['chi2']))
    for h in pars[0]:
        print("{} = {} +/- {}".format(h, pars[2][h], pars[3][h]))

    print("\n- Fixed params:")    
    for h in pars[1]:
        print("{} = {}".format(h, pars[2][h]))

    print("\nCov:")
    for h in pars[4].keys():
        print("{} = {}".format(h, pars[4][h]))

    if cor:
        print("\nCor:")
        for h in pars[4].keys():
            idx1 = h.find('[')
            idx2 = h.find(',')
            idx3 = h.find(']')
            cov = pars[4][h]
            h1 = h[idx1+1:idx2]
            h2 = h[idx2+2:idx3]
            err1 = pars[3][h1]
            err2 = pars[3][h2]
            print("cor[{}, {}] = {}".format(h1, h2, cov/err1/err2))

    return


def h5file_to_latex(file_list, ap_digits=3, at_digits=3, b_eta_lya_digits=4, beta_lya_digits=3, b_lya_digits=4, beff_lya_digits=4, beta_qso_digits=3, b_qso_digits=3, beta_HCD_digits=3, b_HCD_digits=3, f_digits=3, b_si1190_digits=2, b_si1193_digits=2, b_si1207_digits=2, b_si1260_digits=2, b_cv_digits=2, b_hcd_digits=4, beta_hcd_digits=3, a_sky_digits=3, sigma_sky_digits=1, chi2_digits=0, zeff_digits=3, drp_digits=3, header=True):
    '''
    This function print results of picca fitter2 stored in h5 files
    in the latex table format
    '''
    pars = extract_h5file(file_list[0])
    if header:
        print("\\toprule")
        if len(file_list) == 4:
            print("Param\\`etre  & $\\num{0} < z < \\num{2.35}$ & $\\num{2.35} < z < \\num{2.65}$ & $\\num{2.65} < z < \\num{3.05}$ & $\\num{3.05} < z < \\num{10}$ \\\\")
        if len(file_list) == 5:
            print("Param\\`etre  & $\\num{0} < z < \\num{2.35}$ & $\\num{2.35} < z < \\num{2.65}$ & $\\num{2.65} < z < \\num{3.05}$ & $\\num{3.05} < z < \\num{10}$  & $\\num{0} < z < \\num{10}$ \\\\")
        if len(file_list) == 1:
            print("Param\\`etre  & $\\num{0} < z < \\num{10}$ \\\\")
        print("\\midrule")
    if 'ap' in pars[0]:
        row = "$\\apar{} $"
        row += loop_on_h5file(file_list, 'ap', ap_digits)
        row += " \\\\"
        print(row)
    if 'at' in pars[0]:
        row = "$\\aperp{} $"
        row += loop_on_h5file(file_list, 'at', at_digits)
        row += " \\\\"
        print(row)
    if 'bias_eta_LYA' in pars[0]:
        row = "$b_{\\eta, \\mathrm{Ly}\\alpha} $"
        row += loop_on_h5file(file_list, 'bias_eta_LYA', b_eta_lya_digits)
        row += " \\\\"
        print(row)
    if 'beta_LYA' in pars[0]:
        row = "$\\beta_{\\mathrm{Ly}\\alpha} $"
        row += loop_on_h5file(file_list, 'beta_LYA', beta_lya_digits)
        row += " \\\\"
        print(row)
    if 'beta_QSO' in pars[0]:
        row = "$\\beta_{\\mathrm{QSO}} $"
        row += loop_on_h5file(file_list, 'beta_QSO', beta_qso_digits)
        row += " \\\\"
        print(row)
    if 'beta_HCD' in pars[0]:
        row = "$\\beta_{\\mathrm{HCD}} $"
        row += loop_on_h5file(file_list, 'beta_HCD', beta_HCD_digits)
        row += " \\\\"
        print(row)
    if 'growth_rate' in pars[0]:
        row = "$f$"
        row += loop_on_h5file(file_list, 'growth_rate', f_digits)
        row += " \\\\"
        print(row)
    if 'bias_eta_SiII(1190)' in pars[0]:
        row = "$10^3 b_{\\eta, SiII(1190)} $"
        row += loop_on_h5file(file_list, 'bias_eta_SiII(1190)', b_si1190_digits)
        row += " \\\\"
        print(row)
    if 'bias_eta_SiII(1193)' in pars[0]:
        row = "$10^3 b_{\\eta, SiII(1193)} $"
        row += loop_on_h5file(file_list, 'bias_eta_SiII(1193)', b_si1193_digits)
        row += " \\\\"
        print(row)
    if 'bias_eta_SiII(1260)' in pars[0]:
        row = "$10^3 b_{\\eta, SiII(1260)} $"
        row += loop_on_h5file(file_list, 'bias_eta_SiII(1260)', b_si1260_digits)
        row += " \\\\"
        print(row)
    if 'bias_eta_SiIII(1207)' in pars[0]:
        row = "$10^3 b_{\\eta, SiIII(1207)} $"
        row += loop_on_h5file(file_list, 'bias_eta_SiIII(1207)', b_si1207_digits)
        row += " \\\\"
        print(row)
    if 'bias_eta_CIV(eff)' in pars[0]:
        row = "$10^3 b_{\\eta, CIV(\\mathrm{eff})} $"
        row += loop_on_h5file(file_list, 'bias_eta_CIV(eff)', b_cv_digits)
        row += " \\\\"
        print(row)
    if 'bias_hcd' in pars[0]:
        row = "$b_{\\textsc{HCD}} $"
        row += loop_on_h5file(file_list, 'bias_hcd', b_hcd_digits)
        row += " \\\\"
        print(row)
    if 'beta_hcd' in pars[0]:
        row = "$\\beta_{\\textsc{HCD}} $"
        row += loop_on_h5file(file_list, 'beta_hcd', beta_hcd_digits)
        row += " \\\\"
        print(row)
    if 'BB-cf_z_0_10-0-broadband_sky-scale-sky' in pars[0]:
        row = "$10^2 A_{sky} $"
        row += loop_on_h5file(file_list, 'BB-cf_z_0_10-0-broadband_sky-scale-sky', a_sky_digits)
        row += " \\\\"
        print(row)
    if 'BB-cf_z_0_10-0-broadband_sky-sigma-sky' in pars[0]:
        row = "$\\sigma_{sky} $"
        row += loop_on_h5file(file_list, 'BB-cf_z_0_10-0-broadband_sky-sigma-sky', sigma_sky_digits)
        row += " \\\\"
        print(row)
    if 'drp_QSO' in pars[0]:
        row = "$\\Delta_{\\rpar{}, \\textsc{QSO}}$"
        row += loop_on_h5file(file_list, 'drp_QSO', drp_digits)
        row += " \\\\"
        print(row)
    print("\\midrule")
    row = "$\\chi^2$"
    row += loop_on_h5file(file_list, 'chi2', chi2_digits)
    row += " \\\\"
    print(row)
    row = "$z_{\\mathrm{eff}}$"
    row += loop_on_h5file(file_list, 'zeff', zeff_digits)
    row += " \\\\"
    print(row)
    print("\\midrule")
    if 'bias_LYA' in pars[0]:
        row = "$b_{\\mathrm{Ly}\\alpha} $"
        row += loop_on_h5file(file_list, 'bias_LYA', b_lya_digits)
        row += " \\\\"
        print(row)
    if 'beff_LYA' in pars[0]:
        row = "$b_{\\mathrm{eff}, \\mathrm{Ly}\\alpha} $"
        row += loop_on_h5file(file_list, 'beff_LYA', beff_lya_digits)
        row += " \\\\"
        print(row)
    if 'bias_QSO' in pars[0]:
        row = "$b_{\\mathrm{QSO}} $"
        row += loop_on_h5file(file_list, 'bias_QSO', b_qso_digits)
        row += " \\\\"
        print(row)
    if 'bias_HCD' in pars[0]:
        row = "$b_{\\mathrm{HCD}} $"
        row += loop_on_h5file(file_list, 'bias_HCD', b_HCD_digits)
        row += " \\\\"
        print(row)
    print("\\bottomrule")
    # if '' in pars[0]:
    #     row = ""
    #     row += loop_on_h5file(file_list, '', _digits)
    #     row += " \\\\"
    #     print(row)
    return

def loop_on_h5file(file_list, param, digits=5):
    res = ""
    for f in file_list:
        pars = extract_h5file(f)
        val = pars[2][param]
        if param == 'zeff':
            res += " & $ {} $".format(np.round(val,digits))
            continue
        if param == 'chi2':
            res += " & $ {} $".format(np.int32(np.round(val,digits)))
            continue
        err = np.abs(pars[3][param])
        if 'Si' in param or 'CIV' in param:
            val *= 1e3
            err *= 1e3
        if param == 'BB-cf_z_0_10-0-broadband_sky-scale-sky':
            val *= 1e2
            err *= 1e2
        val = np.round(val, digits)
        err = np.round(err,digits)
        res += " & $ {} \\pm {}$".format(val, err)

    return res

class cosmo:
    '''
    From picca.constant.py
    https://github.com/igmhub/picca/blob/master/py/picca/constants.py
    '''
    def __init__(self, Om, Ok=0., Or=0., wl=-1., H0=100.):
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


def kms2mpc(redshift, omega_m=None, omega_k=None, h=None):
    '''
    This function returns the factor to convert km/s to Mpc/h:
    [Mpc/h] = H(z) / [(1+z)*h] [km/s]
    '''
    redshift = np.array(redshift)
    if omega_m is None:
        omega_m = constant.omega_M_0
    if omega_k is None:
        omega_k = constant.omega_k_0
    if h is None:
        h = constant.h
    cosmo_fid = cosmo(Om=omega_m, Ok=omega_k, H0=h*100)
    hubble_z = cosmo_fid.hubble(redshift)
    factor = hubble_z / (1+redshift) / h
    return factor


def bias_err(bias_eta, bias_eta_err, beta, beta_err, cov):
    return bias_eta / beta * np.sqrt((bias_eta_err/bias_eta)**2 + (beta_err/beta)**2 - 2*cov/bias_eta/beta)


def beff_err(bias_eta, bias_eta_err, beta, beta_err, cov, f=0.97):
    db_dbiaseta = f*np.sqrt(1+2/3*beta+1/5*beta**2)/beta
    db_dbeta = (bias_eta*f/beta)*np.sqrt(1+2/3*beta+1/5*beta**2)*((1/3+1/5*beta)/(1+2/3*beta+1/5*beta**2) - 1/beta)
    beff_err = np.sqrt((db_dbiaseta*bias_eta_err)**2 + (db_dbeta*beta_err)**2 + 2*db_dbiaseta*db_dbeta*cov)
    return beff_err

def plot_wedge(da_list, co_list, label_list=None, color_list=constant.colors, title='', marker_list=['.','.','.','.'],mumin=0, mumax=1, rtmin=0, rtmax=200, rpmin=0, rpmax=200, nrt=50, nrp=50, absoluteMu=True, rpow=2):
    """
    da_list is a list of 1D array that contains correlation funtions
    co_list is a list of 2D array that contains covariance matrices
    """
    if not use_picca:
        print("Unable to use this function. Install picca first.")
        sys.exit(1)
    if label_list == None:
        label_list = np.arange(len(da_list)) + 1

    fig, ax = plt.subplots(figsize=(12,8))
    w = picca.wedgize.wedge(mumin=mumin,mumax=mumax, rtmax=rtmax, rpmax=rpmax, rtmin=rtmin, rpmin=rpmin, nrt=nrt, nrp=nrp,absoluteMu=absoluteMu)
    for da, co, lab, col, mrk in zip(da_list, co_list, label_list, color_list, marker_list):
        data_wedge = w.wedge(da,co)
        coef = data_wedge[0]**rpow
        # ax.errorbar(data_wedge[0],coef*data_wedge[1],yerr=coef*np.sqrt(np.diag(data_wedge[2])),fmt='+', label=lab)
        ax.errorbar(data_wedge[0],coef*data_wedge[1],yerr=coef*np.sqrt(np.diag(data_wedge[2])), label=lab, color=col, marker=mrk)

    ax.grid()
    ax.legend()
    ax.set_title(title, fontsize=20)
    ax.set_xlabel(r"$r \, [h^{-1}\mathrm{Mpc}]$",fontsize=20)
    if rpow == 2:
        ax.set_ylabel(r"$r^{2}\xi(r) \, [(h^{-1}\mathrm{Mpc})^2]$",fontsize=20)
    if rpow == 1:
        ax.set_ylabel(r"$r\xi(r) \, [h^{-1}\mathrm{Mpc}]$",fontsize=20)
    if rpow == 0:
        ax.set_ylabel(r"$\xi(r)$",fontsize=20)
    plt.show()
    return


def add_wedge(da, co, errorbar=True, mumin=0, mumax=1, rtmin=0, rtmax=200, rpmin=0, rpmax=200, nrt=50, nrp=50, absoluteMu=True, rpow=2, **kwargs):
    """
    da is a 1D array that contains correlation funtions
    co is a 2D array that contains covariance matrices
    """
    if not use_picca:
        print("Unable to use this function. Install picca first.")
        sys.exit(1)
    w = picca.wedgize.wedge(mumin=mumin,mumax=mumax, rtmax=rtmax, rpmax=rpmax, rtmin=rtmin, rpmin=rpmin, nrt=nrt, nrp=nrp,absoluteMu=absoluteMu)
    data_wedge = w.wedge(da,co)
    coef = data_wedge[0]**rpow
    if errorbar:
        plt.errorbar(data_wedge[0],coef*data_wedge[1],yerr=coef*np.sqrt(np.diag(data_wedge[2])), **kwargs)
    else:
        plt.plot(data_wedge[0],coef*data_wedge[1], **kwargs)
    return

def plot_cf1d(filenames, labels=None, legend=True):
    ymin = 1.e6
    ymax = -1.e6
    if labels==None:
        labels = 1 + np.arange(len(filenames))
    f, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
    for f, label in zip(filenames, labels):
        h = fitsio.FITS(f)
        y = h[1].read()['c1d'] #- in the future replace '1' by '1DCOR'
        y_err = np.sqrt(h[1].read()['v1d'])
        binsize = h[1].read_header()['DLL']
        bins = sp.arange(y.size)
        x = sp.power(10,bins*binsize)
        w = h[1]['nc1d'][:]>0.
        x = x[w]
        y = y[w]
        y_err = y_err[w]
        # ax1.errorbar(x, y, yerr=y_err, fmt='+', label=label)
        ax1.plot(x, y, marker='o', label=label)
        ymin = min(ymin,y.min())
        ymax = max(ymax,y[y!=1.].max())
        # ax2.errorbar(x, y, yerr=y_err, fmt='+', label=label)
        ax2.plot(x, y, marker='o', label=label)
        h.close()
    ax1.set_xlabel(r'$\lambda_{1}/\lambda_{2}$')
    ax1.set_ylabel(r'$\xi^{ff,1D}_{normed}$')
    ax1.legend()
    ax1.grid()
    ax2.set_xlim([0.999,1.1])
    ax2.set_ylim([-0.035,0.025])
    ax2.set_xlabel(r'$\lambda_{1}/\lambda_{2}$')
    ax2.set_ylabel(r'$\xi^{ff,1D}_{normed}$')
    if legend:
        ax2.legend()
    ax2.grid()
    plt.subplots_adjust(hspace=0.4)
    plt.tight_layout()
    plt.show()
    return

def growthRateStructure(z, omega_M_0=constant.omega_M_0):
    omega_m = omega_M_0*(1.+z)**3 / ( omega_M_0*(1.+z)**3+(1.-omega_M_0))
    f = sp.power(omega_m,0.55)
    return f
