import os
import fitsio
import numpy as np
import scipy as sp
from SaclayMocks import util, powerspectrum, constant
from iminuit import Minuit
import time
import matplotlib.pyplot as plt
import numpy.ma as ma
import pyfftw.interfaces.numpy_fft as fft


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
    def __init__(self, indir, z, a_ref, cc, bb=1.58, Nreg=1, pixel=0.2, cell_size=2.19, convergence_factor=1, convergence_criterium=None):
        self.data = {}
        self.mock = {}
        self.fit = {}
        self.mock['a_ref'] = a_ref
        self.mock['indir'] = indir
        self.z = z
        self.Nreg = Nreg
        self.mock['pixel'] = pixel
        self.mock['cell_size'] = cell_size
        self.mock['k_ny_spec'] = np.pi / pixel
        self.mock['k_ny_box'] = np.pi / cell_size
        self.mock['bb'] = bb
        self.mock['cc'] = cc
        self.niter = 0  # current iteration on p1d
        self.convergence_factor = convergence_factor
        self.convergence_criterium = convergence_criterium
        self.converged = False

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
        convert_factor = util.kms2mpc(self.z)
        self.data['k'] = data[:, 1][msk] * convert_factor
        self.data['Pk'] = data[:, 2][msk] / convert_factor
        self.data['Pkerr'] = data[:, 3][msk] / convert_factor
        dk = self.data['k'][1] - self.data['k'][0]
        bins = np.append(self.data['k'] - dk/2, self.data['k'][-1]+dk/2)
        self.data['dk'] = dk
        self.data['bins'] = bins

    def read_mock(self, nfiles=None, debug=False):
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

    def read_p1dmiss(self, filename=None):
        if filename is None:
            filename = self.mock['indir']+"/"+self.p1dmiss_filename()
        print("Reading P1Dmissing {}".format(filename))
        data = fitsio.read(filename, ext=1)
        kk = data['k']
        field = 'P1Dmiss'
        if self.mock['cc'] > 0:
            field += 'RSD'
        pk = data[field]
        # msk = kk > 0
        # msk = kk > 0.12  # prov
        if not hasattr(self, "data['k_fit']"):
            self.fit_data()
        msk = (kk > self.data['k_fit'].min()) & (kk < self.data['k_fit'].max())  # avoid interpolation issue in computation of rr in iterate()
        self.mock['kmiss'] = kk[msk]
        self.mock['p1dmiss'] = pk[msk]

    def fit_data(self, klim=1.6, filename=None):
        k_fit, p1d_fit = util.read_P1D_fit(self.z)
        convert_factor = util.kms2mpc(self.z)
        k_fit *= convert_factor
        p1d_fit /= convert_factor
        ### following line is to correct the fact that the fit on P1D data don't go to 0 at large k
        # msk1 = k_fit < klim
        # if filename is None:
        #     os.path.expandvars("$SACLAYMOCKS_BASE/etc/fit_p1d_mock_z3.4.fits")
        # data = fitsio.read(filename, ext=1)
        # k_mock = data['k']
        # p1d_mock = data['P1D']
        # msk2 = k_mock >= klim
        # x = p1d_fit[np.logical_not(msk1)][0] / p1d_mock[msk2][0]
        # k = np.append(k_fit[msk1], k_mock[msk2])
        # p1d = np.append(p1d_fit[msk1], x*p1d_mock[msk2])
        # p1d_interp = sp.interpolate.interp1d(k, p1d)
        k = k_fit
        p1d = p1d_fit
        p1d_interp = sp.interpolate.interp1d(k, p1d)
        self.data['k_fit'] = k
        self.data['p1d_fit'] = p1d
        self.data['p1d_fit_interp'] = p1d_interp

    def p1dmiss_filename(self, niter=None):
        if niter is None:
            niter=self.niter
        if niter == 0:
            filename = "p1dmiss.fits"
        else:
            filename = "p1dmiss-{}.fits".format(niter)
        return filename

    def p1d_filename(self, niter=None):
        if niter is None:
            niter=self.niter
        filename = "p1d-{}.fits".format(niter)
        return filename

    def gen_small_scales(self, filename=None):
        # if not filename:
        #     filename = self.mock['indir'] + "/" + self.p1dmiss_filename()
        # p1d_data = fitsio.read(filename, ext=1)
        # field = 'P1Dmiss'
        # if self.mock['cc'] > 0:
        #     field += 'RSD'
        # P1Dmissing = sp.interpolate.InterpolatedUnivariateSpline(p1d_data['k'], p1d_data[field])
        self.read_p1dmiss()
        P1Dmissing = sp.interpolate.InterpolatedUnivariateSpline(self.mock['kmiss'], self.mock['p1dmiss'])
        delta_s = []
        nz = len(self.mock['wav'])
        delta_s = np.random.normal(size=self.mock['spectra'].shape)
        delta_sk = fft.rfft(delta_s, axis=1)
        k = np.fft.rfftfreq(nz)*2*self.mock['k_ny_spec']
        # print("k interp: {} to {}".format(p1d_data['k'].min(), p1d_data['k'].max()))
        # print("interp applied from {} to {}".format(k.min(), k.max()))
        if self.mock['kmiss'].max() < k.max() or self.mock['kmiss'].min() > k.min():
            print("WARNING P1D is extrapolated with a spline outside interpolation range /!\ ")
        Pmis = np.maximum(P1Dmissing(k), 0)
        # # Correct the amplitude:
        # Pmis *= (self.sigma_eta_l / self.sigma_l)**2
        delta_sk *= np.sqrt(Pmis/self.mock['pixel'])
        delta_s = fft.irfft(delta_sk, axis=1)
        self.mock['delta_s'] = delta_s

    def compute_p1d(self, a, bins=None, debug=False):
        if bins is None:
            bins = self.data['bins']
        if debug:
            spec = self.mock['spectra']
        else:
            spec = self.comp_spectra(a)
        Fmean = ma.mean(spec)
        p1d = powerspectrum.ComputeP1D(self.mock['pixel']*self.Nreg)
        for s in spec:
            if len(s[s.mask==False]) < self.Nreg: continue
            if self.Nreg > 1:
                flux = util.regroup(s[s.mask==False], self.Nreg)
            else:
                flux = s[s.mask==False]
            delta = flux / Fmean - 1
            p1d.add_spectrum(delta)
        k, pk, pkerr = p1d.P1D(bins)
        msk = k > 0
        self.mock['k'] = k[msk]
        self.mock['p1d'] = pk[msk]
        self.mock['err_p1d'] = pkerr[msk]

    def smooth_p1d(self, k=None, pk=None, pkerr=None):
        if pk is None:
            k = self.mock['k']
            pk = self.mock['p1d']
            pkerr = self.mock['err_p1d']
            ret = False
        else:
            ret = True

        # s parameter for spline:
        s = len(k)*pkerr.mean()**2 * 3  # 1.5 is fine tuning, not to get fluctuation of P1D in fit
        # # Add a point at k=20 to avoid spline to diverge
        # kf = np.concatenate((kf, [self.k_fit[-1]]))
        # pkf = np.concatenate((pkf, [self.p1d_fit[-1]]))
        # pkferr = np.concatenate((pkferr, [pkferr[-1]]))
        # continuation with a straight line

        pk_interp = sp.interpolate.UnivariateSpline(k, pk, s=s)
        # Chi2 check:
        chi2_red = (((pk_interp(k)-pk)/pkerr)**2).sum() / len(pk)
        print("Reduced chi2 for smoothing is {}.".format(chi2_red))
        if ret:
            return pk_interp
        else:
            self.mock['p1d_interp'] = pk_interp

    def comp_spectra(self, a):
        '''
        Compute the spectra using delta_l, delta_s and eta_par
        and given the parameter a
        '''
        delta = self.mock['delta_l'] + self.mock['delta_s']
        spec = util.fgpa(delta, self.mock['eta_par'], self.mock['growthf'],
                         a, self.mock['bb'], self.mock['cc'])
        if spec.mask.size != self.mock['mask_size']:
            print("WARNING: There is a different number of masked value:\n{} initially and {} now"
                  .format(self.mock['mask_size'], spec.mask.size))
        return spec

    def chi2(self, a):
        '''
        Compute the chi2 between the P1D of mock and P1D of data
        '''
        self.compute_p1d(a)
        chi2 = (((self.mock['p1d']-self.data['Pk']) / self.data['Pkerr'])**2).sum()
        return chi2

    def minimize(self, a_init=0.01, a_err=0.1, a_min=0., a_max=3.,tol=1e4, print_level=1):
        '''
        Minimizing procedure to determine the optimal a parameter (for a given redshift)
        '''
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

    def iterate(self, a=None, bins=None, plot=False):
        '''
        Iterative procedure to tune the shape of the 1D powerspectrum
        '''
        # print(self.data['k_fit'].min())
        # print(self.data['k_fit'].max())
        # print(self.mock['kmiss'].min())
        # print(self.mock['kmiss'].max())
        rr = self.data['p1d_fit_interp'](self.mock['kmiss']) / self.mock['p1d_interp'](self.mock['kmiss'])
        p1dmiss = self.mock['p1dmiss'] * (1 + self.convergence_factor*(rr - 1))
        if plot:
            plt.plot(self.mock['kmiss'], self.data['p1d_fit_interp'](self.mock['kmiss']), label='fit data')
            plt.plot(self.mock['kmiss'], self.mock['p1d_interp'](self.mock['kmiss']), label='mock smoothed')
            plt.errorbar(self.mock['k'], self.mock['p1d'], yerr=self.mock['err_p1d'], fmt='.', label='mock')
            plt.plot(self.mock['kmiss'], self.mock['p1dmiss'] / 50, label='p1dmiss_{}'.format(self.niter))
            plt.plot(self.mock['kmiss'], rr, label='ratio rr')
            plt.plot(self.mock['kmiss'], p1dmiss / 50, label='p1dmiss_{}'.format(self.niter+1))
            plt.grid()
            plt.legend()
            plt.xlabel('k [h/Mpc/]')
            plt.title('iteration : {}'.format(self.niter))
            plt.show()

        if self.convergence_criterium:
            print("convergence is {}".format(np.max(np.abs(rr-1))))
            if np.max(np.abs(rr - 1)) < self.convergence_criterium:
                self.converged = True
                print("Iterative procedure has converged in {} iterations.".format(self.niter))
                return

        self.mock['p1dmiss'] = p1dmiss
        self.niter += 1
        self.export_p1dmiss()
        self.gen_small_scales()
        if a is None:
            a = self.fit['a']
        self.compute_p1d(a, bins=bins)
        self.smooth_p1d()
        self.export_p1d()
        print("Iteration {} done.\n".format(self.niter))

    def check_p1d(self, a=None, title='', debug=False):
        if debug:
            self.compute_p1d(1, debug=True)
        else:
            if not a:
                a = self.fit['a']
            self.compute_p1d(a)

        # Pk vs k [h.Mpc^-1]
        f1, ax1 = plt.subplots()
        # ax1.set_yscale('log')
        ax1.set_xlabel('k [h/Mpc]')
        ax1.set_ylabel('Pk')
        ax1.set_title(title)
        ax1.grid()
        ax1.errorbar(self.mock['k'], self.mock['p1d'], yerr=self.mock['err_p1d'], fmt='.', label='mock')
        ax1.errorbar(self.data['k'], self.data['Pk'], yerr=self.data['Pkerr'], fmt='+', label='data')
        # convert_factor = util.kms2mpc(self.z)
        # ax1.plot(k[msk], util.P1Dk(k[msk]/convert_factor, z)/convert_factor, '.', label='fit data')
        ax1.legend()

        # k*Pk vs k [s.km^-1]
        f2, ax2 = plt.subplots()
        ax2.set_yscale('log')
        ax2.set_title(title)
        ax2.set_xlabel('k [s/km]')
        ax2.set_ylabel('k * Pk / pi')
        convert_factor = util.kms2mpc(self.z)
        ax2.errorbar(self.mock['k']/convert_factor, self.mock['k']*self.mock['p1d']/np.pi, yerr=self.mock['k']*self.mock['err_p1d']/np.pi, fmt='.', label='mock')
        ax2.errorbar(self.data['k']/convert_factor, self.data['k']*self.data['Pk']/np.pi, yerr=self.data['k']*self.data['Pkerr']/np.pi, fmt='+', label='data')
        # ax2.errorbar(self.data['k'][msk]/convert_factor, self.data['k'][msk]*self.data['Pk'][msk]/np.pi, yerr=self.data['k'][msk]*self.data['Pkerr'][msk]/np.pi, fmt='+', label='data')
        # ax2.plot(k[msk]/convert_factor, k[msk]*util.P1Dk(k[msk]/convert_factor, z)/convert_factor, '.', label='fit data')
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

    def export(self, outdir):
        np.save(outdir+'/a.npy', self.fit['a'])
        np.save(outdir+'/time.npy', self.fit['time'])
        np.save(outdir+'/chi.npy', self.fit['chi2'])
        np.save(outdir+'/tol.npy', self.fit['tol'])

    def export_p1d(self):
        outfits = fitsio.FITS(self.mock['indir']+self.p1d_filename(), 'rw', clobber=True)
        table = [self.mock['k'], self.mock['p1d'], self.mock['err_p1d']]
        outfits.write(table, names=['k', 'P1D', 'P1Derr'])
        outfits[-1].write_key('pixel', self.mock['pixel'])
        outfits[-1].write_key('c', self.mock['cc'])
        outfits.close()
        print("Wrote {}".format(self.p1d_filename()))

    def export_p1dmiss(self):
        outfits = fitsio.FITS(self.mock['indir']+self.p1dmiss_filename(), 'rw', clobber=True)
        table = [self.mock['kmiss'], self.mock['p1dmiss']]
        field = 'P1Dmiss'
        if self.mock['cc'] > 0:
            field += 'RSD'
        names = ['k', field]
        outfits.write(table, names=names)
        outfits[-1].write_key('pixel', self.mock['pixel'])
        outfits[-1].write_key('c', self.mock['cc'])
        outfits.close()
        print("Wrote {}".format(self.p1dmiss_filename()))
