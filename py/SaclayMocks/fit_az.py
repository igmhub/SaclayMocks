import os
import fitsio
import numpy as np
from SaclayMocks import util, powerspectrum, constant
from iminuit import Minuit
import time
import matplotlib.pyplot as plt
import numpy.ma as ma
from scipy import interpolate


"""
Fitter runs on raw spectra (spectra_merged), skewer transmissions.
You need a directory with at least the boxes, the QSO
Start by difining chi2 with indir and redshift
then read_data()
then read_mock()
and then minimize()
finally, export results of minimisation
can also check plots with check_plot()
"""
class Fitter(object):
    def __init__(self, indir, z, a_ref, cc, bb=1.58, Nreg=1, xx=100., pixel=0.2):
        self.data = {}
        self.mock = {}
        self.mock['a_ref'] = a_ref
        self.mock['indir'] = indir
        self.fit = {}
        self.z = z
        self.Nreg = Nreg
        self.pixel = pixel
        self.xx = xx  # parameter to go from km/s to Mpc/h
        self.mock['bb'] = bb
        self.mock['cc'] = cc
        # self.mock['dd'] = dd

    @classmethod
    def mock_init(Cls, indir, z, beta, cc, dd, a_ref=1e-3, Nhdu=8, bb=1.58, Nreg=1, xx=100., pixel=0.2, kmax=20., dk=0.001):
        '''
        This alternative init produces p1dmiss.fits and
        runs minimize_az which runs MakeSpectra
        and MergeSpectra $Nhdu times with the given parameters.
        '''
        # Running comp_p1dmiss.py:
        print("Running compute_p1dmiss.py...")
        os.system("python compute_p1dmiss.py -beta {} -dd {} -pixel {} -kmax {} -dk {} -outdir {}"
                  .format(beta, dd, pixel, kmax, dk, indir))
        # Producing spectra
        print("Computing spectra...")
        fit_p1d = True
        p1dfile = indir+"/p1dmiss.fits"
        os.system("bash minimize_az.sh {} {} {} {} {} {} {}".format(a_ref, z, indir, Nhdu, cc, fit_p1d, p1dfile))
        print("Done.")
        new = Cls(indir, z, a_ref, cc, bb=bb, xx=xx, pixel=pixel, Nreg=Nreg)
        new.beta = beta
        new.dd = dd
        new.Nhdu = Nhdu
        new.kmax = kmax
        new.dk = dk
        return new

    def read_data(self, filename=None):
        """
        Read Pk file from Nathalie. Si III oscillations have been removed.
        Format is z, k, pk, pkerr, 0, 0, 0
        k in km/s and pk and pkerr in (km/s)**-1 translated to Mpc/h
        """
        if filename is None:
            filename = os.path.expandvars("$SACLAYMOCKS_BASE/etc/pk_fft35bins_noSi.out")
        print("Reading data from: {}".format(filename))
        data = np.loadtxt(filename)
        z = data[:, 0]
        msk = np.where(z == self.z)[0]
        if len(msk) == 0:
            # z from 2.2 to 4.4
            raise ValueError("ERROR -- You entered a wrong redshift. Here is the list of redshift : {}".format(np.unique(z)))

        self.data['z'] = self.z
        self.data['k'] = data[:, 1][msk]*self.xx
        self.data['Pk'] = data[:, 2][msk]/self.xx
        self.data['Pkerr'] = data[:, 3][msk]/self.xx
        dk = self.data['k'][1] - self.data['k'][0]
        bins = np.append(self.data['k'] - dk/2, self.data['k'][-1]+dk/2)
        self.data['dk'] = dk
        self.data['bins'] = bins

    def read_mock(self, nfiles=None, debug=False):
        # print("Reading sigma_l...")
        # sigma_l = fitsio.read_header(self.mock['indir']+"/boxes/box-0.fits", ext=0)['sigma']
        # self.mock['sigma_l'] = sigma_l
        # print("sigma_l is {}".format(sigma_l))
        # print("Computing sigma_s...")
        # p1d_data = fitsio.read(self.mock['indir']+"/"+self.mock['p1dfilename'], ext=1)
        # P1Dmissing = interpolate.InterpolatedUnivariateSpline(p1d_data['k'], self.mock['dd']*p1d_data['P1D'])
        # sigma_s = util.sigma_p1d(P1Dmissing, pixel=self.pixel)
        # self.mock['sigma_s'] = sigma_s
        # print("sigma_s is {}".format(sigma_s))
        print("Reading mock spectra...")
        files = os.listdir(self.mock['indir']+'/spectra_merged/')
        if nfiles:
            files = np.sort(files)[:nfiles]
        first = True
        metadata = []
        spectra = []
        delta_l_list = []
        delta_s_list = []
        eta_par = []
        dt = [('RA', '>f4'), ('DEC', '>f4'), ('Z_noRSD', '>f4'), ('Z', '>f4'),
              ('HDU', '>i8'), ('THING_ID', '>i8'), ('PLATE', '>i8'),
              ('MJD', '>i4'), ('FIBERID', '>i4'), ('PMF', '<U22')]
        for f in files:
            if first:
                wav = fitsio.read(self.mock['indir']+'/spectra_merged/'+f, ext='LAMBDA')
                if not debug:
                    growthf = fitsio.read(self.mock['indir']+'/spectra_merged/'+f, ext='GROWTHF')
                    redshift = fitsio.read(self.mock['indir']+'/spectra_merged/'+f, ext='Z')
                first = False
            data = fitsio.read(self.mock['indir']+'/spectra_merged/'+f, ext='METADATA').astype(dt)
            spec = fitsio.read(self.mock['indir']+'/spectra_merged/'+f, ext='FLUX')
            msk = wav/(1+data['Z']).reshape(-1,1)
            msk = ((msk <= constant.lylimit) | (msk >= constant.lya))
            metadata.append(data)
            spectra.append(ma.array(spec, mask=msk))
            if not debug:
                delta_l = fitsio.read(self.mock['indir']+'/spectra_merged/'+f, ext='DELTA_L')
                delta_s = fitsio.read(self.mock['indir']+'/spectra_merged/'+f, ext='DELTA_S')
                eta = fitsio.read(self.mock['indir']+'/spectra_merged/'+f, ext='ETA_PAR')
                delta_l_list.append(ma.array(delta_l, mask=msk))
                delta_s_list.append(ma.array(delta_s, mask=msk))
                eta_par.append(ma.array(eta, mask=msk))

        self.mock['data'] = np.concatenate(metadata)
        self.mock['wav'] = np.float64(wav)
        self.mock['spectra'] = np.float64(ma.concatenate(spectra))
        if not debug:
            self.mock['growthf'] = np.float64(growthf)
            self.mock['redshift'] =np.float64(redshift)
            self.mock['delta_l'] = np.float64(ma.concatenate(delta_l_list))
            self.mock['delta_s'] = np.float64(ma.concatenate(delta_s_list))
            self.mock['delta'] = self.mock['delta_l'] + self.mock['delta_s']
            self.mock['eta_par'] = np.float64(ma.concatenate(eta_par))
            self.mock['g_field'] = self.mock['delta'] + self.mock['cc']*self.mock['eta_par']
            sigma_l = self.mock['delta_l'].std()
            sigma_s = self.mock['delta_s'].std()
            sigma = np.sqrt(sigma_l**2 + sigma_s**2)
            self.mock['sigma_s'] = sigma_s
            self.mock['sigma_l'] = sigma_l
            self.mock['sigma'] = sigma
            print("Sigma_l = {} ; sigma_s = {} --> sigma = {}".format(sigma_l, sigma_s, sigma))
        self.mock['mask_size'] = self.mock['spectra'].mask.size
        print("Done.")

    def comp_p1d(self, a, bins=None, debug=False):
        if bins is None:
            bins = self.data['bins']
        if debug:
            spec = self.mock['spectra']
        else:
            spec = self.comp_spectra(a)
        Fmean = ma.mean(spec)
        P1D = powerspectrum.ComputeP1D(self.pixel*self.Nreg)
        for s in spec:
            if len(s[s.mask==False]) < self.Nreg: continue
            if self.Nreg > 1:
                flux = util.regroup(s[s.mask==False], self.Nreg)
            else:
                flux = s[s.mask==False]
            delta = flux / Fmean - 1
            P1D.add_spectrum(delta)
        return P1D.P1D(bins)

    def comp_spectra(self, a):
        # spec = ma.exp(-a*self.mock['tau_a'])
        spec = util.fgpa(self.mock['delta'], self.mock['eta_par'], self.mock['growthf'],
                         a, self.mock['bb'], self.mock['cc'])
        if spec.mask.size != self.mock['mask_size']:
            print("WARNING: There is a different number of masked value:\n{} initially and {} now"
                  .format(self.mock['mask_size'], spec.mask.size))
        return spec

    def chi2(self, a):
        k, Pk, Pkerr = self.comp_p1d(a)
        msk = np.where(k > 0)
        chi2 = (((Pk[msk]-self.data['Pk'][msk]) / self.data['Pkerr'][msk])**2).sum()
        return chi2

    def minimize(self, a_init=0.01, a_err=0.1, a_min=0., a_max=3.,tol=1e4, print_level=1):
        t0 = time.time()
        print("Starting minimisation...")
        m=Minuit(self.chi2, a=a_init, error_a=a_err, limit_a=(a_min,a_max), print_level=print_level)
        m.tol = tol
        m.migrad()
        self.fit['a'] = m.values['a']
        self.fit['time'] = time.time()-t0
        self.fit['chi2'] = m.fval
        self.fit['tol'] = m.tol
        print("\nDone. Took {} s\n==> Optimal a is {} with chi2={}".format(time.time()-t0, m.values['a'], m.fval))

    def export(self, outdir):
        np.save(outdir+'/a.npy', self.fit['a'])
        np.save(outdir+'/time.npy', self.fit['time'])
        np.save(outdir+'/chi.npy', self.fit['chi2'])
        np.save(outdir+'/tol.npy', self.fit['tol'])

    def check_p1d(self, a=None, title='', debug=False):
        if debug:
            k, Pk, Pkerr = self.comp_p1d(1, debug=True)
        else:
            if not a:
                a = self.fit['a']
            k, Pk, Pkerr = self.comp_p1d(a)
        msk = np.where(k > 0)

        # Pk vs k [h.Mpc^-1]
        f1, ax1 = plt.subplots()
        # ax1.set_yscale('log')
        ax1.set_xlabel('k [h/Mpc]')
        ax1.set_ylabel('Pk')
        ax1.set_title(title)
        ax1.grid()
        ax1.errorbar(k[msk], Pk[msk], yerr=Pkerr[msk], fmt='.', label='mock')
        ax1.errorbar(self.data['k'][msk], self.data['Pk'][msk], yerr=self.data['Pkerr'][msk], fmt='+', label='data')
        # ax1.plot(k[msk], util.P1Dk(k[msk]/self.xx, z)/self.xx, '.', label='fit data')
        ax1.legend()

        # k*Pk vs k [s.km^-1]
        f2, ax2 = plt.subplots()
        ax2.set_yscale('log')
        ax2.set_title(title)
        ax2.set_xlabel('k [s/km]')
        ax2.set_ylabel('k * Pk / pi')
        ax2.errorbar(k[msk]/self.xx, k[msk]*Pk[msk]/np.pi, yerr=k[msk]*Pkerr[msk]/np.pi, fmt='.', label='mock')
        ax2.errorbar(self.data['k'][msk]/self.xx, self.data['k'][msk]*self.data['Pk'][msk]/np.pi, yerr=self.data['k'][msk]*self.data['Pkerr'][msk]/np.pi, fmt='+', label='data')
        # ax2.plot(k[msk]/self.xx, k[msk]*util.P1Dk(k[msk]/self.xx, z)/self.xx, '.', label='fit data')
        ax2.grid()
        ax2.legend()

        plt.show()

    def check_pdf(self, a=None, bins=100, title=''):
        if not a:
            a = self.fit['a']
        spec = self.comp_spectra(a)
        print("<F> = {}".format(spec.mean()))
        f, ax = plt.subplots()
        ax.set_xlabel('F')
        ax.set_title(title)
        ax.hist(spec[spec.mask==False], bins=bins)
        plt.show()
