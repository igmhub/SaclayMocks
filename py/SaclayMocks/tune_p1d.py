import os
import fitsio
import numpy as np
import numpy.ma as ma
import scipy as sp
from SaclayMocks import util, powerspectrum, fit_az, constant
import pyfftw.interfaces.numpy_fft as fft
import matplotlib.pyplot as plt


'''
Needs in indir: at least boxes and QSO (and p1dmiss_eta_s.fits if dont use mock_init
P1dmiss compute the 1D power spectrum missing, used in MakeSpectra, iteratively
according to the following formulae:
P_n+1(eta_s) = P_n(eta_s)[1+lambda(P_F_data/P_F_mock - 1)]
First initialize P1dmiss: it produces spectra and reads it with read_spectra(),
then it generates eta for the small scales and compute P1D.
Then run iterate() as many times as you want
'''
class P1dmiss(object):
    def __init__(self, indir, z, aa, bb=1.58, cc=1, cell_size=2.19, pixel=0.2, kmax=20., dk=0.001, ll=1, Nhdu=8, xx=100., Nreg=1, comp_spec=False):
        self.z = z
        self.aa = aa
        self.bb = bb
        self.cc = cc
        # self.dd = dd
        self.Nhdu = Nhdu
        self.cell_size = cell_size
        self.k_ny_box = np.pi / cell_size
        self.pixel = pixel
        self.k_ny_spec = np.pi / pixel
        self.Nreg = Nreg
        self.kmax = kmax
        self.dk = dk
        self.ll = ll
        self.xx = xx
        self.indir = indir
        self.niter = 0

        # Read data fit
        self.fit_data()

        #Compute dD/dz:
        dgrowthfile = os.path.expandvars("$SACLAYMOCKS_BASE/etc/dgrowth.fits")
        dgrowth = util.InterpFitsTable(dgrowthfile, 'Z', 'dD/dz')
        self.dD_dz = dgrowth.interp

        # Read iteration 0 P1Dmissing :
        filename = indir+"/"+self.pkmiss_filename()
        print("Reading P1Dmissing {}".format(filename))
        data = fitsio.read(filename, ext=1)
        kk = data['k']
        pk = data['P1D']
        msk = kk >= 0.1  # prov
        self.k = kk[msk]
        self.pkmiss = pk[msk]
        print("Done.")

        if comp_spec:
            # Compute mock spectra
            fit_p1d = True
            print("Producing spectra...")
            os.system("bash minimize_az.sh {} {} {} {} {} {}".format(aa, z, indir, Nhdu, cc, fit_p1d))
            print("Done.")

        # Read mock spectra
        print("Reading spectra...")
        self.read_spectra()
        print("Done.")

        # Compute sigmas
        print("Computing sigma...")
        # sigma_l = 1.186
        sigma_l = self.delta_l.std()
        print("sigma_l = {}".format(sigma_l))
        data = fitsio.read(indir+"/p1dmiss.fits", ext=1)
        p1dmissing = sp.interpolate.InterpolatedUnivariateSpline(data['k'], data['P1D'])
        sigma_s = util.sigma_p1d(p1dmissing, pixel)
        print("sigma_s = {}".format(sigma_s))
        sigma = np.sqrt(sigma_l**2 + sigma_s**2)
        self.sigma_l = sigma_l
        self.sigma = sigma
        self.sigma_s = sigma_s
        print("sigma = {}".format(sigma))
        print("sigma_g = {}".format(self.g_field.std()))
        # sigma_eta_l = ma.std(self.eta_l)
        # self.sigma_eta_l = sigma_eta_l
        # print("sigma eta_l = {}".format(sigma_eta_l))

        # Generate small scales
        print("Generating small scales...")
        self.gen_small_scales()
        print("Done.")
        # sigma_eta_s = ma.std(self.eta_s)
        # self.sigma_eta_s = sigma_eta_s
        # print("sigma eta_s = {}".format(sigma_eta_s))
        tau = util.taubar_over_a(self.sigma, self.growthf, self.bb)
        self.taubar_over_a = tau
        print("Taubar_over_a = {}".format(tau))

        # Compute p1d for 0th iteration
        print("Computing P1D...")
        self.compute_pkf()
        self.smooth_pkf()
        self.export_pkf()
        print("Done.")

    @classmethod
    def mock_init(Cls, indir, z, aa, beta, cc, dd, bb=1.58, cell_size=2.19, pixel=0.2, kmax=20., dk=0.001, ll=1, Nhdu=8, xx=100., Nreg=1):
        '''
        This method is an alternative __init__ : it computes the p1dmiss file and the spectra files,
        according to the parameters.
        '''
        # Running comp_p1dmiss.py:
        print("Running compute_p1dmiss.py...")
        os.system("python compute_p1dmiss.py -beta {} -dd {} -pixel {} -kmax {} -dk {} -outdir {}"
                  .format(beta, dd, pixel, kmax, dk, indir))
        print("Producing spectra...")
        fit_p1d = True
        p1dfile = indir+"/p1dmiss.fits"
        os.system("bash minimize_az.sh {} {} {} {} {} {} {}".format(aa, z, indir, Nhdu, cc, fit_p1d, p1dfile))
        print("Done.")
        new = Cls(indir, z, aa, bb=bb, cc=cc, cell_size=cell_size, pixel=pixel, kmax=kmax, dk=dk, ll=ll, Nhdu=Nhdu, xx=xx, Nreg=Nreg)
        new.beta = beta
        new.dd = dd
        return new

    def pkmiss_filename(self, niter=None):
        if niter is None:
            niter=self.niter
        if niter == 0:
            filename = "p1dmiss.fits"
        else:
            filename = "p1dmiss-{}.fits".format(niter)
        return filename

    def pkf_filename(self, niter=None):
        if niter is None:
            niter=self.niter
        filename = "p1d-{}.fits".format(niter)
        return filename

    def fit_data(self,klim=1.6, filename=None):
        k_fit, p1d_fit = util.read_P1D_fit(self.z)
        k_fit *= self.xx
        p1d_fit /= self.xx
        msk1 = k_fit < klim
        if filename is None:
            os.path.expandvars("$SACLAYMOCKS_BASE/etc/fit_p1d_mock_z3.4.fits")
        data = fitsio.read(filename, ext=1)
        k_mock = data['k']
        p1d_mock = data['P1D']
        msk2 = k_mock >= klim
        x = p1d_fit[np.logical_not(msk1)][0] / p1d_mock[msk2][0]
        k = np.append(k_fit[msk1], k_mock[msk2])
        p1d = np.append(p1d_fit[msk1], x*p1d_mock[msk2])
        p1d_interp = sp.interpolate.interp1d(k, p1d)
        self.k_fit = k
        self.p1d_fit = p1d
        self.p1d_fit_interp = p1d_interp

    def read_spectra(self):
        print("Reading mock spectra...")
        files = os.listdir(self.indir+'/spectra_merged/')
        first = True
        spectra = []
        delta_l = []
        delta_s = []
        eta_par = []
        # eta_l = []
        # eta_s = []
        for f in files:
            if first:
                wav = fitsio.read(self.indir+'/spectra_merged/'+f, ext='LAMBDA')
                growthf = fitsio.read(self.indir+'/spectra_merged/'+f, ext='GROWTHF')
                redshift = fitsio.read(self.indir+'/spectra_merged/'+f, ext='Z')
                first = False
            data = fitsio.read(self.indir+'/spectra_merged/'+f, ext='METADATA')
            spec = fitsio.read(self.indir+'/spectra_merged/'+f, ext='FLUX')
            delta_l_tmp = fitsio.read(self.indir+'/spectra_merged/'+f, ext='DELTA_L')
            delta_s_tmp = fitsio.read(self.indir+'/spectra_merged/'+f, ext='DELTA_S')
            eta_par_tmp = fitsio.read(self.indir+'/spectra_merged/'+f, ext='ETA_PAR')
            # eta_l_tmp = fitsio.read(self.indir+'/spectra_merged/'+f, ext='ETA_L')
            # eta_s_tmp = fitsio.read(self.indir+'/spectra_merged/'+f, ext='ETA_S')
            msk = wav/(1+data['Z']).reshape(-1,1)
            msk = ((msk <= constant.lylimit) | (msk >= constant.lya))
            spectra.append(ma.array(spec, mask=msk))
            delta_l.append(ma.array(delta_l_tmp, mask=msk))
            delta_s.append(ma.array(delta_s_tmp, mask=msk))
            eta_par.append(ma.array(eta_par_tmp, mask=msk))
            # eta_l.append(ma.array(eta_l_tmp, mask=msk))
            # eta_s.append(ma.array(eta_s_tmp, mask=msk))

        self.wav = wav
        self.growthf = growthf
        self.redshift =redshift
        self.spectra_init = ma.concatenate(spectra)
        self.spectra = ma.concatenate(spectra)
        self.delta_l = ma.concatenate(delta_l)
        self.delta_s = ma.concatenate(delta_s)
        self.eta_par = ma.concatenate(eta_par)
        self.g_field = self.delta_l + self.delta_s + self.cc * self.eta_par
        # self.eta_l = ma.concatenate(eta_l)
        # self.eta_s = ma.concatenate(eta_s)
        self.mask_size = self.spectra.mask.size
        print("Done.")

    def iterate(self, plot=False):
        rr = self.p1d_fit_interp(self.k) / self.pkf_interp(self.k)
        pkmiss = self.pkmiss * (1 + self.ll*(rr - 1))
        if plot:
            plt.plot(self.k, self.p1d_fit_interp(self.k))
            plt.plot(self.k, self.pkf_interp(self.k))
            plt.plot(self.kf, self.pkf, '+')
            plt.plot(self.k, self.pkmiss)
            plt.plot(self.k, rr)
            plt.plot(self.k, pkmiss)
            plt.show()
        self.pkmiss = pkmiss
        self.niter += 1
        self.export_pkmiss()
        self.gen_small_scales()
        self.compute_pkf()
        self.smooth_pkf()
        self.export_pkf()

    def gen_small_scales(self, filename=None):
        if not filename:
            filename = self.indir + self.pkmiss_filename()
        p1d_data = fitsio.read(filename, ext=1)
        P1Dmissing = sp.interpolate.InterpolatedUnivariateSpline(p1d_data['k'], p1d_data['P1D'])
        if not hasattr(self, 'spectra'):
            self.read_spectra()
        delta_s = []
        nz = len(self.wav)
        delta_s = np.random.normal(size=self.spectra.shape)
        delta_sk = fft.rfft(delta_s, axis=1)
        k = np.fft.rfftfreq(nz)*2*self.k_ny_spec
        # print("k interp: {} to {}".format(p1d_data['k'].min(), p1d_data['k'].max()))
        # print("interp applied from {} to {}".format(k.min(), k.max()))
        if p1d_data['k'].max() < k.max() or p1d_data['k'].min() > k.min():
            print("WARNING P1D is extrapolated with a spline outside interpolation range /!\ ")
        Pmis = np.maximum(P1Dmissing(k), 0)
        # # Correct the amplitude:
        # Pmis *= (self.sigma_eta_l / self.sigma_l)**2
        delta_sk *= np.sqrt(Pmis/self.pixel)
        delta_s = fft.irfft(delta_sk, axis=1)
        self.delta_s = delta_s

    def comp_spectra(self):
        delta = self.delta_l + self.delta_s
        # eta_par = self.eta_l + self.eta_s
        # eta_par *= self.dD_dz(self.redshift) / self.dD_dz(0)
        # spectra = np.exp(-self.aa * (np.exp(self.bb*self.growthf*delta) -
        #                              self.cc*(1+self.redshift)*self.taubar_over_a*self.eta_par))
        spectra = util.fgpa(delta, self.eta_par, self.growthf, self.aa, self.bb, self.cc)
        if spectra.mask.size != self.mask_size:
            print("WARNING: There is a different number of masked value:\n{} initially and {} now"
                  .format(self.mask_size, spectra.mask.size))

        self.spectra = spectra
        self.g_field = self.delta_l + self.delta_s + self.cc*self.eta_par
        return spectra

    def compute_pkf(self, bins=300):
        # compute P1D of the transmitted flux fraction F
        spectra = self.comp_spectra()
        f_mean = ma.mean(spectra)
        p1d = powerspectrum.ComputeP1D(self.pixel*self.Nreg)
        for spec in spectra:
            if len(spec[spec.mask==False]) < self.Nreg: continue
            if self.Nreg > 1:
                flux = util.regroup(spec[spec.mask==False], self.Nreg)
            else:
                flux = spec[spec.mask==False]
            delta = flux / f_mean - 1
            p1d.add_spectrum(delta)
        kf, pkf, pkferr = p1d.P1D(bins//self.Nreg)
        self.kf = kf
        self.pkf = pkf
        self.pkferr =pkferr

    def smooth_pkf(self, kf=None, pkf=None, pkferr=None):
        if pkf is None:
            kf = self.kf
            pkf = self.pkf
            pkferr = self.pkferr
            ret = False
        else:
            ret = True

        # s parameter for spline:
        s = len(kf)*pkferr.mean()**2 * 1.5  # 1.5 is fine tuning, not to get fluctuation of P1D in fit
        # # Add a point at k=20 to avoid spline to diverge
        # kf = np.concatenate((kf, [self.k_fit[-1]]))
        # pkf = np.concatenate((pkf, [self.p1d_fit[-1]]))
        # pkferr = np.concatenate((pkferr, [pkferr[-1]]))
        # continuation with a straight line

        pkf_interp = sp.interpolate.UnivariateSpline(kf, pkf, s=s)
        # Chi2 check:
        chi2_red = (((pkf_interp(kf)-pkf)/pkferr)**2).sum() / len(pkf)
        print("Reduced chi2 for smoothing is {}.".format(chi2_red))
        if ret:
            return pkf_interp
        else:
            self.pkf_interp = pkf_interp

    def export_pkmiss(self):
        outfits = fitsio.FITS(self.indir+self.pkmiss_filename(), 'rw', clobber=True)
        table = [self.k, self.pkmiss]
        outfits.write(table, names=['k', 'P1D'])
        outfits[-1].write_key('kmax', self.kmax)
        outfits[-1].write_key('dk', self.dk)
        outfits[-1].write_key('pixel', self.pixel)
        outfits[-1].write_key('c', self.cc)
        outfits.close()
        print("Wrote {}".format(self.pkmiss_filename()))

    def export_pkf(self):
        outfits = fitsio.FITS(self.indir+self.pkf_filename(), 'rw', clobber=True)
        table = [self.kf, self.pkf, self.pkferr]
        outfits.write(table, names=['k', 'P1D', 'P1Derr'])
        outfits[-1].write_key('kmax', self.kmax)
        outfits[-1].write_key('dk', self.dk)
        outfits[-1].write_key('pixel', self.pixel)
        outfits[-1].write_key('c', self.cc)
        outfits.close()
        print("Wrote {}".format(self.pkf_filename()))

    def check_pkf(self, niter=None):
        if niter is None:
            niter = self.niter
        print("Reading {}".format(self.pkf_filename(niter)))
        data_pkf = fitsio.read(self.indir+self.pkf_filename(niter), ext=1)
        k = np.linspace(data_pkf['k'].min(), data_pkf['k'].max(), 1000)
        # Read data
        k_data, pk_data, pkerr_data = util.read_P1D(self.z)
        k_data *= self.xx
        pk_data /= self.xx
        pkerr_data /= self.xx
        ## Read data fit
        # k_fit, pk_fit = util.read_P1D_fit(self.z)
        # k_fit *= self.xx
        # pk_fit /= self.xx
        f, ax = plt.subplots()
        ax.set_xlabel('k [h/Mpc]')
        ax.set_ylabel('P1D')
        ax.errorbar(data_pkf['k'], data_pkf['P1D'], yerr=data_pkf['P1Derr'], fmt='+', label='mock-{}'.format(niter))
        ax.plot(k, self.smooth_pkf(data_pkf['k'], data_pkf['P1D'], data_pkf['P1Derr'])(k), label='mock fit-{}'.format(niter))
        ax.errorbar(k_data, pk_data, yerr=pkerr_data, fmt='+', label='data')
        ax.plot(self.k_fit, self.p1d_fit, label='data fit')
        try:
            print("Reading {}".format(self.pkf_filename(niter-1)))
            data_pkf = fitsio.read(self.indir+self.pkf_filename(niter-1), ext=1)
            ax.errorbar(data_pkf['k'], data_pkf['P1D'], yerr=data_pkf['P1Derr'], fmt='+', label='mock-{}'.format(niter-1))
            ax.plot(k, self.smooth_pkf(data_pkf['k'], data_pkf['P1D'], data_pkf['P1Derr'])(k), label='mock fit-{}'.format(niter-1))
        except:
            print("Can't read P1D_n-1")
        ax.grid()
        ax.legend()
        plt.show()

    def check_pkmiss(self, niter=None):
        if niter is None:
            niter = self.niter
        print("Reading {}".format(self.pkmiss_filename(niter)))
        data_pkmiss = fitsio.read(self.indir+self.pkmiss_filename(niter), ext=1)
        pkmiss_interp = sp.interpolate.InterpolatedUnivariateSpline(data_pkmiss['k'], data_pkmiss['P1D'])
        k_interp = np.linspace(0, data_pkmiss['k'].max(), 1e5)
        f, ax = plt.subplots()
        ax.set_xlabel('k [h/Mpc]')
        ax.set_ylabel('P1Dmissing')
        ax.plot(data_pkmiss['k'], data_pkmiss['P1D'], '.', label='n={}'.format(niter))
        ax.plot(k_interp, pkmiss_interp(k_interp))
        try:
            print("Reading {}".format(self.pkmiss_filename(niter-1)))
            data_pkmiss = fitsio.read(self.indir+self.pkmiss_filename(niter-1), ext=1)
            pkmiss_interp = sp.interpolate.InterpolatedUnivariateSpline(data_pkmiss['k'], data_pkmiss['P1D'])
            k_interp = np.linspace(0, data_pkmiss['k'].max(), 1e5)
            ax.plot(data_pkmiss['k'], data_pkmiss['P1D'], '.', label='n={}'.format(niter-1))
            ax.plot(k_interp, pkmiss_interp(k_interp))
        except:
            print("Cannot read P1Dmissing_n-1")
        ax.grid()
        ax.legend()
        plt.show()

    def check_pdf(self, bins=100, title=''):
        spec = self.comp_spectra()
        print("<F> = {}".format(spec.mean()))
        f, ax = plt.subplots()
        ax.set_xlabel('F')
        ax.set_title(title)
        ax.hist(spec[spec.mask==False], bins=bins)
        plt.show()
