#!/usr/bin/env python
import fitsio
import os, sys
import numpy as np
import scipy as sp
import argparse
import time
from SaclayMocks import util, constant
import pyfftw
import pyfftw.interfaces.numpy_fft as fft
import glob
# import matplotlib.pyplot as plt


# @profile
def main():
# if True:
    t_init = time.time()

    # ........... Read arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-inDir", help="dir for spectra fits file")
    parser.add_argument("-outDir", help="dir for merged spectra fits file")
    parser.add_argument("-i", type=int, help="Treated slice")
    parser.add_argument("-aa", type=float, help="a param in FGPA. Default is to read a(z) from etc/params.fits", default=-1)
    parser.add_argument("-bb", type=float, help="b param in FGPA. Default 1.58", default=1.58)
    parser.add_argument("-cc", type=float, help="c param in FGPA. Default is to read c(z) from etc/params.fits.", default=-1)
    parser.add_argument("-paramfile", help="fits file for parameters, default is etc/params.fits", default=None)
    parser.add_argument("-p1dfile", help="P1Dmissing fits file, default is etc/pkmiss_interp.fits.gz", default=None)
    parser.add_argument("-pixsize", type=float, help="spectrum pixsize in Mpc/h, default 0.2", default=0.2)
    parser.add_argument("-nside", type=int, help="nside for healpix. Default 16", default=16)
    parser.add_argument("-nest", help="If True, healpix scheme is nest. Default True", default='True')
    parser.add_argument("-rsd", help="If True, rsd are added, default True", default='True')
    parser.add_argument("-addnoise", help="If True, small scales are added, default True", default='True')
    parser.add_argument("-dla", help="If True, store delta and growth skewers, default False", default='False')
    parser.add_argument("-zfix", type=float, help="Specify a redshift to evaluate parameters in FGPA, if None parameter are taken at the absorber redshift", default=None)
    parser.add_argument("--fit-p1d", help="If True, do the fitting procedure. Default False", default='False')
    parser.add_argument("--store-g", help="If True, store delta_l, delta_s and eta_par. Default False", default='False')
    parser.add_argument("-seed", type=int, help="specify a seed", default=None)
    parser.add_argument("--check-id", help="If True, check if the spectra ID matches the QSO ID by looking at (ra,dec), default True", default='True')
    parser.add_argument("-ncpu", type=int, help="number of cpu, default = 2", default=2)
    args = parser.parse_args()

    inpath = args.inDir
    outpath = args.outDir
    islice = args.i
    nside = args.nside
    ncpu = args.ncpu
    nest_option = util.str2bool(args.nest)
    rsd = util.str2bool(args.rsd)
    add_noise = util.str2bool(args.addnoise)
    dla = util.str2bool(args.dla)
    fit_p1d = util.str2bool(args.fit_p1d)
    store_g = util.str2bool(args.store_g)
    check_id = util.str2bool(args.check_id)
    pixsize = args.pixsize
    k_ny = np.pi / pixsize
    if args.paramfile is None:
        filename = os.path.expandvars("$SACLAYMOCKS_BASE/etc/params.fits")

    if args.aa > 0:
        aa = args.aa
        print("a  parameter has been fixed to {}".format(aa))
    else:
        a_of_z = util.InterpFitsTable(filename, 'z', 'a')
    if args.bb > 0:
        bb = args.bb
        print("b  parameter has been fixed to {}".format(bb))
    else:
        b_of_z = util.InterpFitsTable(filename, 'z', 'b')
    if args.cc > 0:
        cc = args.cc
        print("c  parameter has been fixed to {}".format(cc))
    else:
        c_of_z = util.InterpFitsTable(filename, 'z', 'c')

    seed = args.seed
    if seed is None:
        seed = np.random.randint(2**31 -1, size=1)[0]
        np.random.seed(seed)
        print("Seed has not been specified. Seed is set to {}".format(seed))
    else:
        seed = seed + islice
        np.random.seed(seed)
        print("Specified seed is {}".format(seed))

    print("Arguments read.\nFits file read from {}\nOut file will be written in {}\nTreated slice : {}".format(inpath, outpath, islice))

    # .......... Compute growth factor
    Om = constant.omega_M_0
    growthf_24 = util.fgrowth(2.4, Om)

    # .......... Load QSO files
    if check_id:
        qso_data = []
        files = glob.glob(inpath+"/../qso/*")
        for f in files:
            try:
                d = fitsio.read(f, ext=1)
                qso_data.append(d)
            except IOError:
                print("Can't read file: {}".format(f))
        qso_data = np.concatenate(qso_data)

    # .......... Load P1D missing
    filename = args.p1dfile
    if filename is None:
        filename = os.path.expandvars("$SACLAYMOCKS_BASE/etc/pkmiss_interp.fits.gz")
    print("Reading P1D file {}".format(filename))
    if fit_p1d:
        store_g = True
        p1d_data = fitsio.read(filename, ext=1)
        field = 'P1Dmiss'
        if rsd: field += 'RSD'
        p1dmiss = sp.interpolate.InterpolatedUnivariateSpline(p1d_data['k'], p1d_data[field])
    else:
        p1dmiss = util.InterpP1Dmissing(filename)
        sigma_s_interp = sp.interpolate.interp1d(fitsio.read(filename, ext='z'), fitsio.read(filename, ext='sigma'))

    # ........... Open fits files
    print("Opening fits files...")
    fits = []

    for f in os.listdir(inpath):
        index1 = f.rfind('-')+1
        index2 = f.find('.')
        if f[index1:index2] == str(islice):
            try:
                fits.append(fitsio.FITS(inpath+'/'+f))
                print("{} opened".format(f))
            except:
                print("*WARNING* Fits file {} cannot be read".format(f))

    if len(fits) == 0:
        print("No fits file opened. Exit.")
        sys.exit()
    else:
        print("Fits files opened - {} s".format(time.time() - t_init))

    # ........... Read ID
    print("Reading IDs...")
    IDs = []
    RA = []
    DEC = []
    Z_QSO_NO_RSD = []
    Z_QSO_RSD = []
    PLATE = []
    MJD = []
    FIBERID = []
    PMF = []
    LAMBDA = []
    delta_l = []
    eta_par = []
    velo_par = []
    redshift = []
    cpt1 = 0

    first = True
    for f in fits:
        data = f[1].read()
        if len(data) == 0:
            print("\n*WARNING* fits table empty:")
            print(f[1])
            continue
        if first:
            first = False
            header = f[1].read_header()
            z0 = header['z0']
            pixel = header['pixel']
            NX = header['NX']
            dmax = header['dmax']
            ra0 = header['ra0']
            dec0 = header['dec0']
            npixeltot = header['Npixel']

        IDs.append(data['THING_ID'])
        RA.append(data['RA'])
        DEC.append(data['DEC'])
        Z_QSO_NO_RSD.append(data['Z_noRSD'])
        Z_QSO_RSD.append(data['Z'])
        PLATE.append(data['PLATE'])
        MJD.append(data['MJD'])
        FIBERID.append(data['FIBERID'])
        PMF.append(data['PMF'])
        wav = f['LAMBDA'].read()
        delta = f['DELTA_L'].read()
        if rsd:
            eta = f['ETA_PAR'].read()
            if dla:
                velo = f['VELO_PAR'].read()
        zzz = f['REDSHIFT'].read()
        for i in range(len(wav)):
            msk = wav[i] > 0
            LAMBDA.append(wav[i][msk])
            delta_l.append(delta[i][msk])
            if rsd:
                eta_par.append(eta[i][msk])
                if dla:
                    velo_par.append(velo[i][msk])
            redshift.append(zzz[i][msk])
        cpt1 += len(data)
        f.close()

    IDs = np.concatenate(IDs)
    RA = np.concatenate(RA)
    DEC = np.concatenate(DEC)
    Z_QSO_NO_RSD = np.concatenate(Z_QSO_NO_RSD)
    Z_QSO_RSD = np.concatenate(Z_QSO_RSD)
    PLATE = np.concatenate(PLATE)
    MJD = np.concatenate(MJD)
    FIBERID = np.concatenate(FIBERID)
    PMF = np.concatenate(PMF)
    LAMBDA = np.array(LAMBDA)
    delta_l = np.array(delta_l)
    if rsd:
        eta_par = np.array(eta_par)
        if dla:
            velo_par = np.array(velo_par)
    redshift = np.array(redshift)
    healpix = util.radec2pix(nside, RA, DEC, nest=nest_option)
    print("IDs read - {} s".format(time.time()-t_init))

    #...............................    get wisdom to save time on FFT
    wisdom_path = os.path.expandvars("$SACLAYMOCKS_BASE/etc/")
    wisdomFile = wisdom_path+"wisdom1D_"+str(ncpu)+".npy"

    if os.path.isfile(wisdomFile) :
        pyfftw.import_wisdom(sp.load(wisdomFile))
        save_wisdom = False
    else :
        print("{f} file not found. Saving wisdom file to {f}".format(f=wisdomFile))
        save_wisdom = True

    # .......... Merge spectra
    print("Merging spectra...")
    names = ['RA', 'DEC', 'Z_noRSD', 'Z', 'HDU', 'THING_ID', 'PLATE', 'MJD', 'FIBERID', 'PMF']
    hlist = [{'name':"z0", 'value':z0, 'comment':"redshift of box center"},
             {'name':"NX", 'value':NX},
             {'name':"dmax", 'value':dmax},
             {'name':"ra0", 'value':ra0, 'comment':"right ascension of box center"},
             {'name':"dec0", 'value':dec0, 'comment':"declination of box center"}]
    cpt2 = 0
    cpt3 = 0
    t2 = time.time()
    timer = 0
    for pix in np.unique(healpix):
        cut = (healpix == pix)  # select only ID in pix healpix
        spectra_list = []
        ra_list = []
        dec_list = []
        Z_QSO_NO_RSD_list = []
        Z_QSO_RSD_list = []
        HDU_list = []
        ID_list = []
        PLATE_list = []
        MJD_list = []
        FIBERID_list = []
        PMF_list = []
        if dla:
            delta_list = []
            velo_list = []
        if store_g:
            delta_list = []
            eta_list = []
            noise_list = []
        for ID in np.unique(IDs[cut]):
            # if ID != 3134001249: continue
            msk = np.where((IDs[cut] == ID))[0]
            if len(msk) > 0:
                if check_id:
                    if not np.isin(ID, qso_data['THING_ID']):
                        print("WARNING ID: {} didn't match any QSO ID".format(ID))
                        continue
                    # mm = np.where((qso_data['RA'] == RA[cut][msk][0]) & (qso_data['DEC'] == DEC[cut][msk][0]))[0]
                    # if len(mm) == 0:
                    #     print("No matching QSO for RA: {}; DEC: {}; ID: {}".format(RA[cut][msk][0], DEC[cut][msk][0], ID))
                    #     continue
                    # if qso_data['THING_ID'][mm] != ID:
                    #     print("QSO matched (ra,dec) but has different ID: {}; spectra_id: {}".format(qso_data['THING_ID'][mm], ID))
                    #     continue
                wav_tmp = np.concatenate(LAMBDA[cut][msk])
                arg_wav_sorted = np.argsort(wav_tmp)
                wav_tmp = wav_tmp[arg_wav_sorted]
                if args.zfix:
                    z = args.zfix * np.ones_like(wav_tmp)
                else:
                    z = np.concatenate(redshift[cut][msk])[arg_wav_sorted]
                    z_0 = np.ones_like(z) * 2.2466318099484273  # prov
                if args.aa <= 0:
                    aa = a_of_z.interp(z)
                if args.bb <= 0:
                    bb = b_of_z.interp(z)
                if args.cc <= 0:
                    cc = c_of_z.interp(z)
                growthf_tmp = growthf_24*(1+2.4) / (1+z)
                # growthf_tmp = growthf_24*(1+2.4) / (1+z_0)  # prov
                if rsd:
                    eta_par_tmp = np.concatenate(eta_par[cut][msk])[arg_wav_sorted]
                delta_l_tmp = np.concatenate(delta_l[cut][msk])[arg_wav_sorted]

                # Generate small scales
                if add_noise:
                    timer_init = time.time()
                    wav_rf = wav_tmp / (1+Z_QSO_RSD[cut][msk][0])
                    mmm = np.where((wav_rf<constant.lya) & (wav_rf>constant.lylimit))[0]
                    if len(mmm) > 0:
                        nz = 256
                        while (nz < len(wav_tmp)) : nz *= 2
                        nz *= 2  # make sure that the size of delta_s is larger than the size of the forest, to avoid spurious correlations
                        delta_s = np.random.normal(size=nz)   # latter, produce directly in k space
                        delta_sk = fft.rfftn(delta_s, threads=ncpu)
                        k = np.fft.rfftfreq(nz) * 2 * k_ny
                        zeff = z[mmm].mean()
                        # zeff = z_0[mmm].mean()  # prov
                        if fit_p1d:
                            pmis = p1dmiss(k)
                        else:
                            pmis = p1dmiss(zeff, k)
                            pmis[pmis<0] = 0
                        delta_sk *= np.sqrt(pmis/pixsize)
                        delta_s = fft.irfftn(delta_sk, threads=ncpu)
                        delta_s = delta_s[0:len(wav_tmp)]
                        if not fit_p1d:  # correct the z dependence
                            delta_s *= sigma_s_interp(z) / sigma_s_interp(zeff)
                        delta = delta_l_tmp + delta_s
                        # delta = delta_l_tmp  # prov
                    else:
                        delta = delta_l_tmp
                        if store_g:
                            delta_s = np.zeros_like(delta_l_tmp)
                    timer += time.time() - timer_init
                else:
                    delta = delta_l_tmp

                # Apply FGPA:
                if rsd:
                    spec = util.fgpa(delta, eta_par_tmp, growthf_tmp, aa, bb, cc)
                else:
                    spec = np.exp(-aa * np.exp(bb * growthf_tmp * delta))

                if len(spec) != npixeltot:
                    print("/!\ WARNING Spectrum hasn't the nominal lenght /!\ \n  ID: {} in healpix: {} and z: {}\n    Lenght {} != {}".format(ID, pix, Z_QSO_RSD[cut][msk][0], len(spec), npixeltot))
                    continue
                if len(wav_tmp) != npixeltot:
                    print("/!\ WARNING Wavelength hasn't the nominal lenght /!\ \n  ID: {} in healpix: {} and z: {}\n    Lenght {} != {}".format(ID, pix, Z_QSO_RSD[cut][msk][0], len(wav_tmp), npixeltot))
                    continue
                wav = wav_tmp
                if len(growthf_tmp) != npixeltot:
                    print("/!\ WARNING Growth factor hasn't the nominal lenght /!\ \n  ID: {} in healpix: {} and z: {}\n    Lenght {} != {}".format(ID, pix, Z_QSO_RSD[cut][msk][0], len(growthf_tmp), npixeltot))
                    continue
                growthf = growthf_tmp

                spectra_list.append(spec)
                if dla:
                    delta_list.append(delta_l_tmp)
                    if rsd:
                        vpar = np.concatenate(velo_par[cut][msk])[arg_wav_sorted]
                    else:
                        vpar = np.zeros_like(delta_l_tmp)
                    velo_list.append(vpar)
                if store_g:
                    delta_list.append(delta_l_tmp)
                    eta_list.append(eta_par_tmp)
                    noise_list.append(delta_s)

                ra_list.append(RA[cut][msk][0])
                dec_list.append(DEC[cut][msk][0])
                Z_QSO_NO_RSD_list.append(Z_QSO_NO_RSD[cut][msk][0])
                Z_QSO_RSD_list.append(Z_QSO_RSD[cut][msk][0])
                HDU_list.append(islice)
                ID_list.append(ID)
                PLATE_list.append(PLATE[cut][msk][0])
                MJD_list.append(MJD[cut][msk][0])
                FIBERID_list.append(FIBERID[cut][msk][0])
                PMF_list.append(PMF[cut][msk][0])
                if len(msk) > 1:
                    cpt2 += 1
        if len(ra_list) == 0: continue
        outfits = fitsio.FITS(outpath+'/spectra_merged-{}-{}.fits.gz'.format(pix, islice), 'rw', clobber=True)
        table = [np.array(ra_list), np.array(dec_list),
                 np.array(Z_QSO_NO_RSD_list), np.array(Z_QSO_RSD_list),
                 np.array(HDU_list), np.array(ID_list),
                 np.array(PLATE_list), np.array(MJD_list),
                 np.array(FIBERID_list), np.array(PMF_list)]
        outfits.write(table, names=names, header=hlist, extname='METADATA')
        outfits.write(np.float32(wav), extname='LAMBDA')
        outfits.write(np.float32(spectra_list), extname='FLUX')
        if dla:
            outfits.write(np.float32(delta_list), extname='DELTA')
            outfits.write(np.float32(growthf), extname='GROWTHF')
            outfits.write(np.float32(velo_list), extname='VELO_PAR')
        if store_g:
            outfits.write(np.float32(delta_list), extname='DELTA_L')
            outfits.write(np.float32(eta_list), extname='ETA_PAR')
            outfits.write(np.float32(growthf), extname='GROWTHF')
            outfits.write(np.float32(z), extname='Z')
            outfits.write(np.float32(noise_list), extname='DELTA_S')
        outfits.close()
        cpt3 += len(ra_list)

    print("Merging and writting done. {} s".format(time.time() - t2))
    # Save wisdom
    if save_wisdom:
        print("Saving wisdom file in {}".format(wisdomFile))
        sp.save(wisdomFile, pyfftw.export_wisdom())
        save_wisdom = False
    print("Spectra merged and fits file saved.")
    print("{} initial forests.".format(cpt1))
    print("{} forest mergers.".format(cpt2))
    print("{} total forest written.".format(cpt3))
    print("Slice {} done. Took {}s".format(islice, time.time()-t_init))
    print("FFT timer = {} s".format(timer))

if __name__ == "__main__":
    main()
