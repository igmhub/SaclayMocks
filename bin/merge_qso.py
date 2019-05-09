#!/usr/bin/env python
import fitsio
import numpy as np
import argparse
import time
import os
from SaclayMocks import util, constant
import glob


# if True:
def main():
    t_init = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument("-inDir", help="Directory where the QSO files are located")
    parser.add_argument("-outDir", help="Directory where the out file is saved")
    parser.add_argument("-nside", type=int, help="nside for healpix. Default 16", default=16)
    parser.add_argument("-nest", help="If True, healpix scheme is nest. Default True", default='True')
    parser.add_argument("-random", type=str, help="If True, generate randoms. Default: False", default='False')
    parser.add_argument("-prod", help="True: use mock prod architecture ; False: use a special directory. Default is True", default='True')
    parser.add_argument("-zmin", type=float, help="minimal redshift for drawing QSO, default is 0", default=0)
    parser.add_argument("-zmax", type=float, help="maximal redshift for drawing QSO, default is 10", default=10)
    parser.add_argument("-dgrowthfile", help="dD/dz file, default etc/dgrowth.fits", default=None)

    args = parser.parse_args()
    nside = args.nside
    nest_option = util.str2bool(args.nest)
    random_cond = util.str2bool(args.random)
    prod = util.str2bool(args.prod)
    outDir = args.outDir
    print("Reading QSO fits files...")

    if not random_cond:
        names = ["Z_QSO_NO_RSD", "Z_QSO_RSD", "RA", "DEC", "HDU", "THING_ID", "PLATE", "MJD", "FIBERID", "PMF", "PIXNUM", "MOCKID"]
        outfits = fitsio.FITS(outDir+'/master.fits', 'rw', clobber=True)
    else:
        names = ["Z", "RA", "DEC", "HDU", "THING_ID", "PLATE", "MJD", "FIBERID", "PMF", "PIXNUM", "MOCKID"]
        outfits = fitsio.FITS(outDir+'/master_randoms.fits', 'rw', clobber=True)

    if prod:
        if random_cond:
            files = glob.glob(args.inDir+"/*/randoms/randoms*.fits")
        else:
            files = glob.glob(args.inDir+"/*/qso/QSO*.fits")
    else:
        files = os.listdir(args.inDir)

    ra = []
    dec = []
    hdu = []
    THING_ID =[]
    plate = []
    mjd = []
    fiber = []
    pmf = []
    if not random_cond:
        zzz_norsd = []
        zzz_rsd = []
    else:
        zzz = []
    for f in files:
        data = fitsio.read(f, ext=1)
        if not random_cond:
            zzz_norsd.append(data['Z_QSO_NO_RSD'])
            zzz_rsd.append(data['Z_QSO_RSD'])
        else:
            zzz.append(data['Z'])

        ra.append(data['RA'])
        dec.append(data['DEC'])
        hdu.append(data['HDU'])
        THING_ID.append(data['THING_ID'])
        plate.append(data['PLATE'])
        mjd.append(data['MJD'])
        fiber.append(data['FIBERID'])
        pmf.append(data['PMF'])

    ra = np.concatenate(ra)
    dec = np.concatenate(dec)
    if not random_cond:
        zzz_norsd = np.concatenate(zzz_norsd)
        zzz_rsd = np.concatenate(zzz_rsd)
    else:
        zzz = np.concatenate(zzz)
    hdu = np.concatenate(hdu)
    THING_ID = np.concatenate(THING_ID)
    plate = np.concatenate(plate)
    mjd = np.concatenate(mjd)
    fiber = np.concatenate(fiber)
    pmf = np.concatenate(pmf)
    mockid = THING_ID
    pix = util.radec2pix(nside, ra, dec, nest=nest_option)
    print("Reading done. {} s".format(time.time()-t_init))
    if not random_cond:
        array_list = [np.float32(zzz_norsd), np.float32(zzz_rsd), np.float32(ra), np.float32(dec), np.int32(hdu), THING_ID, plate, np.int32(mjd), np.int32(fiber), pmf, np.int32(pix), mockid]
    else:
        array_list = [np.float32(zzz), np.float32(ra), np.float32(dec), np.int32(hdu), THING_ID, plate, np.int32(mjd), np.int32(fiber), pmf, np.int32(pix), mockid]

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

    # Write distance relations
    cosmo_fid = {'omega_M_0':constant.omega_M_0,
                 'omega_lambda_0':constant.omega_lambda_0,
                 'omega_k_0':constant.omega_k_0, 'h':constant.h}
    R_of_z = dist.quick_distance_function(dist.comoving_distance, return_inverse=False, **cosmo_fid)
    growthf_24 = util.fgrowth(2.4, Om)
    redshift = np.linspace(args.zmin, args.zmax, 10000)
    growthf = growthf_24*(1+2.4) / (1+redshift)
    rz = R_of_z(redshift)*constant.h  # in comobile coordinates
    if args.dgrowthfile is None:
        filename = os.path.expandvars("$SACLAYMOCKS_BASE/etc/dgrowth.fits")
    if constant.omega_M_0 == fitsio.read_header(filename, ext=1)['OM']:
        Dgrowth = util.InterpFitsTable(filename, 'Z', 'dD/dz')
    else:
        raise ValueError("Omega_M_0 in SaclayMocks.constant ({}) != OM in {}".format(
            constant.omega_M_0, fitsio.read_header(filename, ext=1)['OM']))
    fz = Dgrowth.interp(redshift)
    table = [redshift, rz, growthf, fz]
    names = ['Z', 'DC', 'G', 'F']
    outfits.write(table, names=names, extname='COSMODIS')

    # Write input power spectrum
    filename = os.path.expandvars("$SACLAYMOCKS_BASE/etc/PlanckDR12.fits")
    data = fitsio.read(filename, ext=1)
    head = fitsio.read_header(filename, ext=1)
    outfits.write(data, header=head)

    outfits.close()
    if random_cond:
        print("Merged fits file written in :"+outDir+'/master_randoms.fits')
    else:
        print("Merged fits file written in :"+outDir+'/master.fits')
    print("{} QSO written".format(len(ra)))
    print("Took {}s".format(time.time()-t_init))


if __name__ == "__main__":
    main()
