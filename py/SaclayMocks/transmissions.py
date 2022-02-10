import fitsio
import os
import glob
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from scipy import stats as stats
from SaclayMocks import powerspectrum, util, constant


SMALL_SIZE = 16
MEDIUM_SIZE = 20
BIGGER_SIZE = 22
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.rc('figure', figsize=(9,7))

class ReadTransmission(object):
    '''
    Create hdu list of transmission files and extract informations
    hdu are stored opened and can be read with the methods read_*
    If the list of transmission files is too long, use alt_init()
    to instance the class : the method will read directly the info,
    store them in attributes and close de fits files.
    '''
    def __init__(self, inDir, read_metadata=True, read_dla=True, nfiles=None, lmin=constant.lylimit, lmax=constant.lya):
        metadata = []
        transmission = []
        dla = []
        #print(constant.lylimit,lmin,lmax)
        files = glob.glob(inDir+"/*/*/transmission*")
        if nfiles is not None and nfiles < len(files):
            files = files[:np.int32(nfiles)]
        wav = fitsio.read(files[0], ext='WAVELENGTH')
        for f in files:
            print("Reading : {}".format(f))
            data = fitsio.read(f, ext='METADATA')
            spec = fitsio.read(f, ext='TRANSMISSION')
            msk = wav/(1+data['Z']).reshape(-1,1)
            msk = ((msk <= lmin) | (msk >= lmax))
            transmission.append(ma.array(spec, mask=msk))
            if read_metadata:
                metadata.append(data)
            if read_dla:
                data_dla = fitsio.read(f, ext='DLA')
                if len(data_dla) > 0:
                    dla.append(data_dla)

        self.wavelength = wav
        transmission = ma.concatenate(transmission)
        self.transmission = transmission
        self.nspec = len(self.transmission)
        self.npix = len(self.wavelength)
        if read_metadata:
            self.metadata = np.concatenate(metadata)
        if read_dla:
            self.dla = np.concatenate(dla)

    def comp_wav_rf(self):
        print("Computing wavelength in rest frame")
        wav_rf = self.wavelength / (self.metadata['Z'].reshape(-1,1) + 1)  # (nspec, npix)
        msk = ((wav_rf >= lmax) | (wav_rf <= lmin))
        self.wav_rf = ma.array(wav_rf, mask=msk)
        print("Done.")

    def trans_of_wav(self, bins=200, plot=True, title='', err=True):
        if not hasattr(self, 'wav_rf'):
            self.comp_wav_rf()
        if err:
            wav_rf, trans, err_trans = util.MakeProfileHisto(self.wav_rf[self.wav_rf.mask==False], self.transmission[self.transmission.mask==False], bins=bins)
        else:
            wav_rf, trans = stats.binned_statistic(self.wav_rf[self.wav_rf.mask==False],
                    [self.wav_rf[self.wav_rf.mask==False], self.transmission[self.transmission.mask==False]],
                                                   bins=bins, statistic='mean').statistic
        if plot:
            f, ax = plt.subplots()
            if err:
                ax.errorbar(wav_rf, trans, yerr=err_trans)
            else:
                ax.errorbar(wav_rf, trans)
            ax.set_ylim(0,1)
            ax.set_xlabel(r'$\lambda_{RF}$', fontsize=20)
            ax.set_ylabel('T', fontsize=20)
            ax.set_title(title, fontsize=20)
            ax.grid()
            plt.show()
        else:
            if err:
                return wav_rf, trans, err_trans
            else:
                return wav_rf, trans

    def trans_of_z(self, plot=True, title=''):
        z = self.wavelength / constant.lya - 1  # (npix)
        t = ma.mean(self.transmission, axis=0)
        if plot:
            zz = np.linspace(1.8, 3.6, 10000)
            # tt = np.exp(-0.0051*(1+zz)**3.2)
            # tt2 = np.exp(-0.0028*(1+zz)**3.45)
            tt2 = np.exp(-0.00211*(1+zz)**3.7)
            f, ax = plt.subplots()
            ax.set_ylim(0, 1)
            ax.plot(z, t, label='mock')
            # ax.plot(zz, tt, label='data', linestyle='--')
            ax.plot(zz, tt2, label='data', linestyle='--')
            ax.set_xlabel(r'$z$')
            ax.set_ylabel(r'$<F>(z)$')
            ax.set_title(title)
            ax.grid()
            # ax.legend()
            plt.show()
        else:
            return z, t

    def pdf_of_F(self, Nreg=1, bins=400, title=''):
        if Nreg > 1:
            loss = self.npix % Nreg
            spec = rebin(self.transmission[:,:-loss], shape=(self.nspec, self.npix // Nreg))
        else:
            spec = self.transmission
        f, ax = plt.subplots()
        ax.set_xlabel('F')
        ax.set_title(title+" - Nreg = {}".format(Nreg))
        ax.hist(spec[spec.mask==False], bins=bins, density=True)
        plt.show()

    def p1d(self, add_data=False, redshift=None, Nreg=1, bins=300, title='', filename=None, pixel=0.2, plot=False):
        if add_data:
            # P1D of data
            if filename is None:
                filename = os.path.expandvars("$SACLAYMOCKS_BASE/etc/pk_fft35bins_noSi.out")
            print("Reading P1D data from: {}".format(filename))
            data = np.loadtxt(filename)
            z = data[:, 0]
            msk = np.where(z == redshift)[0]
            if len(msk) == 0:
                # z from 2.2 to 4.4
                raise ValueError("ERROR -- You entered a wrong redshift. Here is the list of redshift : {}".format(np.unique(z)))

            convert_factor = util.kms2mpc(redshift)
            k_data = data[:, 1][msk] * convert_factor
            pk_data = data[:, 2][msk] / convert_factor
            pkerr_data = data[:, 3][msk] / convert_factor

        # P1D of mock
        print("Computing P1D of mocks")
        Fmean = ma.mean(self.transmission)
        P1D = powerspectrum.ComputeP1D(pixel*Nreg)
        for s in self.transmission:
            if len(s[s.mask==False]) < Nreg: continue
            if Nreg > 1:
                flux = regroup(s[s.mask==False], Nreg)
            else:
                flux = s[s.mask==False]
            delta = flux / Fmean - 1
            P1D.add_spectrum(delta)
        k, pk, pkerr = P1D.P1D(bins//Nreg)
        print("Done")
        if plot:
            print(". Plotting...")
            f1, ax1 = plt.subplots()
            ax1.set_xlabel('k [h/Mpc]')
            ax1.set_ylabel('Pk')
            ax1.set_title(title+" - Nreg = {}".format(Nreg))
            ax1.grid()
            ax1.errorbar(k, pk, yerr=pkerr, fmt='.', label='mock')
            if add_data:
                ax1.errorbar(k_data,pk_data, yerr=pkerr_data, fmt='+', label='data')
            ax1.legend()
            plt.show()
        else:
            return k, pk, pkerr

    def footprint(self, bins=200, title=''):
        # if not hasattr(self, 'ra'):
        #     self.read_meta('RA')
        # if not hasattr(self, 'dec'):
        #     self.read_meta('DEC')
        f, ax = plt.subplots()
        ax.set_xlabel('R.A. [deg]', fontsize=20)
        ax.set_ylabel('sin(Dec)', fontsize=20)
        ax.set_title(title, fontsize=20)
        h = ax.hist2d(self.metadata['RA'], np.sin(np.radians(self.metadata['DEC'])), bins=bins)
        cbar = f.colorbar(h[3], ax=ax)
        cbar.set_label('#QSO', fontsize=20, rotation=90)
        f.tight_layout()
        plt.show()

    def z_footprint(self, bins=200, title='', zmin=1.3, zmax=3.6):
        # if not hasattr(self, 'ra'):
        #     self.read_meta('RA')
        # if not hasattr(self, 'dec'):
        #     self.read_meta('DEC')
        # if not hasattr(self, 'z'):
        #     self.read_meta('Z')
        res = stats.binned_statistic_2d(self.metadata['RA'], np.sin(np.radians(self.metadata['DEC'])),
                    self.metadata['Z'], bins=bins)
        X, Y = np.meshgrid(res.x_edge, res.y_edge)
        f, ax = plt.subplots()
        ax.set_xlabel('R.A. [deg]', fontsize=20)
        ax.set_ylabel('sin(Dec)', fontsize=20)
        ax.set_title(title, fontsize=20)
        cm = ax.pcolormesh(X, Y, np.nan_to_num(res.statistic.T), vmin=zmin, vmax=zmax)
        cbar = f.colorbar(cm)
        cbar.set_label('<z>', fontsize=20, rotation=90)
        f.tight_layout()
        plt.show()
