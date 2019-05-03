#!/usr/bin/env python
import numpy as np
import healpy as hp
import fitsio
from fitsio import FITS
import argparse
import os
import sys
from time import time
from SaclayMocks import util
from SaclayMocks import constant
from memory_profiler import profile


def iterfiles(root, prefix):
    for d in os.listdir(root):
        if d.startswith('chunk'):
            for f in os.listdir(root+'/'+d+'/spectra_merged/'):
                if f.startswith(prefix):
                    yield os.path.join(root, d, 'spectra_merged', f)


# @profile
def main():
# if True:
    t_init = time()

    parser = argparse.ArgumentParser()
    parser.add_argument("-inDir", help="spectra file, default DesiMocks/spectra/", default="DesiMocks/spectra/")
    parser.add_argument("-outDir", help="out directory. Default DesiMocks/output/", default="DesiMocks/output/")
    parser.add_argument("-job", type=int, help="index of current job", default=0)
    parser.add_argument("-ncpu", type=int, help="total number of cpu", default=64)
    parser.add_argument("-nside", type=int, help="nside for healpix. Default 16", default=16)
    parser.add_argument("-nest", help="If True, healpix scheme is nest. Default True", default='True')
    parser.add_argument("-prod", help="True: use mock prod architecture ; False: use a special directory. Default is True", default='True')
    parser.add_argument("-dla", help="If True, store delta and growth skewers, default False", default='False')
    args = parser.parse_args()

    overwrite = True
    inDir = args.inDir
    nside = args.nside
    npixel = hp.nside2npix(nside)
    nest_option = util.str2bool(args.nest)
    dla_cond = util.str2bool(args.dla)
    if nest_option:
        nest_val = 'T'
    outDir = args.outDir
    prod = util.str2bool(args.prod)
    job = args.job
    ncpu = args.ncpu
    if job == ncpu-1:
        pixels = range(int(job*npixel/ncpu), npixel)
    else:
        pixels = range(int(job*npixel/ncpu), int((job+1)*npixel/ncpu))

    print("Treated healpix pixels: {}".format(pixels))

    # Read catalogs
    if dla_cond:
        t0 = time()
        print("Reading DLA catalog...")
        dla_cat = fitsio.read(outDir+"/master_DLA.fits", ext=1)
        print("Done. {} s".format(time()-t0))

    # Select files to read
    print("Getting the treated files...")
    t0 = time()
    fits = []
    fits_pixels = []
    files = iterfiles(inDir, 'spectra_merged')
    for f in files:
        i = int(f[:].find('spectra_merged-')) + 15
        j = i + int(f[i:].find('-'))
        if int(f[i:j]) in pixels:
            fits_pixels.append(int(f[i:j]))
    print("Done. {} s".format(time() - t0))
    if len(fits_pixels) == 0:
        print("No file found for these healpix pixels")
        sys.exit()

    names = ['RA', 'DEC', 'Z_noRSD', 'Z', 'MOCKID']
    hlist = [{'name':"HPXNSIDE", 'value':nside, 'comment':"healpix nside parameter"},
             {'name':"HPXNEST", 'value':nest_val, 'comment':"healpix scheme"},
             {'name':"OL", 'value':constant.omega_lambda_0, 'comment':"Omega_lambda_0"},
             {'name':"OM", 'value':constant.omega_M_0, 'comment':"Omega_M_0"},
             {'name':"OK", 'value':constant.omega_k_0, 'comment':"Omega_k_0"},
             {'name':"H0", 'value':constant.h*100, 'comment':"H0"}]
    dla_names = ['Z_DLA_NO_RSD', 'Z_DLA_RSD', 'N_HI_DLA', 'MOCKID', 'DLAID']

    t0 = time()
    cpt = 0
    for pix in np.unique(fits_pixels):
        fname = outDir+'/{}/{}/transmission-{}-{}.fits.gz'.format(pix//100, pix, nside, pix)
        if not overwrite:
            if os.path.isfile(fname): continue
        outfits = FITS(fname, 'rw', clobber=True)
        ra = []
        dec = []
        fluxes = []
        Z_QSO_NO_RSD = []
        Z_QSO_RSD = []
        ID = []
        z_dla = []
        z_dla_rsd = []
        n_hi_dla = []
        mock_id = []
        dla_id = []
        files = iterfiles(inDir, 'spectra_merged-{}-'.format(pix))
        first = True
        for f in files:
            try :
                fits = fitsio.FITS(f)
            except:
                print("*WARNING* Fits file {} cannot be read".format(f))
                continue
            data = fits['METADATA'].read()
            ra.append(data['RA'])
            dec.append(data['DEC'])
            Z_QSO_NO_RSD.append(data['Z_noRSD'])
            Z_QSO_RSD.append(data['Z'])
            ID.append(data['THING_ID'])
            if first:
                wavelength = fits['LAMBDA'].read()
                first = False
            fluxes.append(fits['FLUX'].read())
            fits.close()
        if len(ra) == 0:
            continue
        cpt += len(ra)
        ra = np.concatenate(ra)
        dec = np.concatenate(dec)
        Z_QSO_NO_RSD = np.concatenate(Z_QSO_NO_RSD)
        Z_QSO_RSD = np.concatenate(Z_QSO_RSD)
        ID = np.concatenate(ID)
        fluxes = np.concatenate(fluxes)
        table = [np.float32(np.array(ra)), np.float32(np.array(dec)),
                 np.float32(np.array(Z_QSO_NO_RSD)), np.float32(np.array(Z_QSO_RSD)), np.array(ID)]

        # get dla from catalog
        if dla_cond:
            cut = dla_cat['PIXNUM'] == pix
            dla_cat_tmp = dla_cat[cut]
            for i in ID:
                msk = np.where(dla_cat_tmp['MOCKID'] == i)[0]
                if len(msk) < 1: continue  # no dla for this id
                z_dla.append(dla_cat_tmp['Z_DLA_NO_RSD'][msk])
                z_dla_rsd.append(dla_cat_tmp['Z_DLA_RSD'][msk])
                n_hi_dla.append(dla_cat_tmp['N_HI_DLA'][msk])
                mock_id.append(dla_cat_tmp['MOCKID'][msk])
                dla_id.append(dla_cat_tmp['DLAID'][msk])
            if len(z_dla) > 0:
                z_dla = np.concatenate(z_dla)
                z_dla_rsd = np.concatenate(z_dla_rsd)
                n_hi_dla = np.concatenate(n_hi_dla)
                mock_id = np.concatenate(mock_id)
                dla_id = np.concatenate(dla_id)
                dla_table = [np.float32(z_dla), np.float32(z_dla_rsd),
                             np.float32(n_hi_dla), mock_id, dla_id]
            else:
                x = np.array([], dtype=np.float32)
                dla_table = [x, x, x, x, x]

        # Write to outfile
        outfits.write(table, names=names, header=hlist, extname='METADATA')  # METADATA HDU
        outfits[-1].write_key("HPXPIXEL", pix, comment='Healpix pixel')
        outfits[-1].write_key("NSIDE", nside, comment='Healpix parameter')
        outfits.write(np.float32(wavelength), extname='WAVELENGTH')  # WAVELENGTH HDU
        outfits.write(np.float32(np.array(fluxes)), extname='TRANSMISSION')  # TRANSMISSION HDU
        if dla_cond:
            outfits.write(dla_table, names=dla_names, extname='DLA')  # DLA HDU
        outfits.close()

    print("{} Sorted spectra written. {} s".format(cpt, time() - t0))
    print("Job {} done. Took {}s".format(job, time()-t_init))

if __name__ == "__main__":
    main()
