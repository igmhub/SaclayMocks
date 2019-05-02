import fitsio
import numpy as np
import argparse
import time
import os
from SaclayMocks import util, constant


def iterfiles(root, random_cond):
    for d in os.listdir(root):
        if d.startswith('chunk'):
            print(d)
            if random_cond:
                for f in os.listdir(root+'/'+d+'/randoms/'):
                    if f.startswith('randoms'):
                        yield os.path.join(root, d, 'randoms', f)
            else:
                for f in os.listdir(root+'/'+d+'/qso/'):
                    if f.startswith('QSO'):
                        yield os.path.join(root, d, 'qso', f)


if True:
# def main():
    t_init = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument("-inDir", help="Directory where the QSO files are located")
    parser.add_argument("-outDir", help="Directory where the out file is saved")
    parser.add_argument("-nside", type=int, help="nside for healpix. Default 16", default=16)
    parser.add_argument("-nest", help="If True, healpix scheme is nest. Default True", default='True')
    parser.add_argument("-random", type=str, help="If True, generate randoms. Default: False", default='False')
    parser.add_argument("-prod", help="True: use mock prod architecture ; False: use a special directory. Default is True", default='True')

    args = parser.parse_args()
    nside = args.nside
    nest_option = util.str2bool(args.nest)
    random_cond = util.str2bool(args.random)
    prod = util.str2bool(args.prod)
    inDir = args.inDir
    outDir = args.outDir
    print("Merging QSO fits files...")

    if not random_cond:
        names = ["Z_QSO_NO_RSD", "Z_QSO_RSD", "RA", "DEC", "HDU", "THING_ID", "PLATE", "MJD", "FIBERID", "PMF", "PIXNUM", "MOCKID"]
        outfits = fitsio.FITS(outDir+'/master.fits', 'rw', clobber=True)
    else:
        names = ["Z", "RA", "DEC", "HDU", "THING_ID", "PLATE", "MJD", "FIBERID", "PMF", "PIXNUM", "MOCKID"]
        outfits = fitsio.FITS(outDir+'/master_randoms.fits', 'rw', clobber=True)

    if prod:
        files = iterfiles(inDir, random_cond)
    else:
        files = os.listdir(inDir)

    first = True
    for f in files:
        data = fitsio.read(f, ext=1)
        if not random_cond:
            zzz_norsd = data['Z_QSO_NO_RSD']
            zzz_rsd = data['Z_QSO_RSD']
        else:
            zzz = data['Z']

        ra = data['RA']
        dec = data['DEC']
        hdu = data['HDU']
        THING_ID = data['THING_ID']
        plate = data['PLATE']
        mjd = data['MJD']
        fiber = data['FIBERID']
        pmf = data['PMF']
        pix = util.radec2pix(nside, ra, dec, nest=nest_option)
        mockid = THING_ID
        if not random_cond:
            array_list = [np.float32(zzz_norsd), np.float32(zzz_rsd), np.float32(ra), np.float32(dec), np.int32(hdu), THING_ID, plate, np.int32(mjd), np.int32(fiber), pmf, np.int32(pix), mockid]
        else:
            array_list = [np.float32(zzz), np.float32(ra), np.float32(dec), np.int32(hdu), THING_ID, plate, np.int32(mjd), np.int32(fiber), pmf, np.int32(pix), mockid]

        if first:
            first = False
            outfits.write(array_list, names=names, extname='CATALOG')
            outfits[1].write_key("RA", None, comment="right ascension in degrees")
            outfits[1].write_key("DEC", None, comment="declination in degrees")
            if random_cond:
                outfits[1].write_key("Z_QSO_NO_RSD", None, comment="redshift without rsd")
                outfits[1].write_key("Z_QSO_RSD", None, comment="redshift with rsd")
            else:
                outfits[1].write_key("Z", None, comment="redshift")

            outfits[1].write_key("HDU", None, comment="Slice index of the QSO")
            outfits[1].write_key("PIXNUM", None, comment="healpix pixel number")
            outfits[1].write_key("Om", constant.omega_M_0, comment="Omega_M_0")
            outfits[1].write_key("Ol", constant.omega_lambda_0, comment="Omega_lambda_0")
            outfits[1].write_key("Ob", constant.omega_b_0, comment="Omega_baryon_0")
            outfits[1].write_key("Ok", constant.omega_k_0, comment="Omega_k_0")
            outfits[1].write_key("h", constant.h, comment="h")
            outfits[1].write_key("z0", constant.z0, comment="box center redshift")

        else:
            outfits[1].append(array_list, names=names)

    if random_cond:
        print("Merged fits file written in :"+outDir+'/master_randoms.fits')
    else:
        print("Merged fits file written in :"+outDir+'/master.fits')
    print("Took {}s".format(time.time()-t_init))


# if __name__ == "__main__":
#     main()
