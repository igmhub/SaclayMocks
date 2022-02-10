#!/usr/bin/env python
import fitsio
import time
import argparse
import os
import numpy as np
from SaclayMocks import util
import glob


# if True:
def main():
    t_init = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument("-indir", help="Root directory of the various Chunks")
    parser.add_argument("-outdir", help="Path to out directory")
    parser.add_argument("-nside", type=int, help="nside for healpix. Default 16", default=16)
    parser.add_argument("-nest", help="If True, healpix scheme is nest. Default True", default='True')
    parser.add_argument("-random",type = str, default='False', help="If True, generate randoms")
    parser.add_argument('--nhi-low-cut', type = float, default=None,
                        help='cut HCDs with log(n_HI) < nhi-low_cut')
    parser.add_argument('--nhi-high-cut', type = float, default=None,
                        help='cut HCDs with log(n_HI) > nhi-high_cut')

    args = parser.parse_args()
    nside = args.nside
    nest_option = util.str2bool(args.nest)
    random_cond = util.str2bool(args.random)
    if not random_cond:
        filenames = args.indir+"/*/dla"
    else:
        filenames = args.indir+"/*/dla_randoms"
    if args.nhi_low_cut is not None and args.nhi_high_cut is not None:
        filenames += "_nhi_{}_{}".format(args.nhi_low_cut, args.nhi_high_cut)
    filenames += ".fits"
    print("Merging DLA files: {}".format(filenames))
    files = glob.glob(filenames)

    # Read files
    print("Reading files...")
    mockid = []
    z = []
    dz = []
    nhi = []
    z_qso = []
    z_qso_rsd = []
    ra = []
    dec = []
    for f in files:
        data = fitsio.read(f, ext=1)
        mockid.append(data['MOCKID'])
        z.append(data['Z_DLA'])
        dz.append(data['DZ_DLA'])
        nhi.append(data['N_HI_DLA'])
        z_qso_rsd.append(data['Z_QSO'])
        z_qso.append(data['Z_QSO_NO_RSD'])
        ra.append(data['RA'])
        dec.append(data['DEC'])

    mockid = np.concatenate(mockid)
    z = np.concatenate(z)
    dz = np.concatenate(dz)
    nhi = np.concatenate(nhi)
    z_qso = np.concatenate(z_qso)
    z_qso_rsd = np.concatenate(z_qso_rsd)
    ra = np.concatenate(ra)
    dec = np.concatenate(dec)
    pixnum = util.radec2pix(nside, ra, dec, nest=nest_option)
    dla_id = np.arange(len(ra)) + 1
    z_rsd = z + dz
    t1 = time.time()
    print("Reading done. {} s".format(t1-t_init))

    # # Read master file
    # print("Getting QSO infos...")
    # master = fitsio.read(args.indir+"/output/master.fits", ext=1)

    # qso_ids = master['THING_ID']
    # id_idx = util.find_A_in_B(mockid, qso_ids)
    # ra = master['RA'][id_idx]
    # dec = master['DEC'][id_idx]
    # z_qso = master['Z_QSO_NO_RSD'][id_idx]
    # z_qso_rsd = master['Z_QSO_RSD'][id_idx]
    # pixnum = master['PIXNUM'][id_idx]
    # dla_id = np.arange(len(ra)) + 1
    # z_rsd = z + dz
    # t2 = time.time()
    # print("Done. {} s".format(t2 - t1))


    # Write DLA catalog
    t2 = time.time()
    print("Writting merged file...")
    if not random_cond:
        filename = args.outdir + "/master_DLA.fits"
        table = [ra, dec, z_qso, z_qso_rsd, z, z_rsd, nhi, mockid, dla_id, pixnum]
        names = ['RA', 'DEC', 'Z_QSO_NO_RSD', 'Z_QSO_RSD', 'Z_DLA_NO_RSD', 'Z_DLA_RSD', 'N_HI_DLA', 'MOCKID', 'DLAID', 'PIXNUM']
    else:
        filename = args.outdir + "/master_DLA_randoms.fits"
        table = [ra, dec, z_qso, z_qso_rsd, z, mockid]
        names = ['RA', 'DEC', 'Z_QSO_NO_RSD', 'Z_QSO_RSD', 'Z', 'MOCKID']
    outfits = fitsio.FITS(filename, 'rw', clobber=True)
    outfits.write(table, names=names, extname='DLACAT')
    print("Done. {} s. Written {} DLAs in {}".format(time.time()-t2, len(ra), filename))
    print("Took {} s".format(time.time()-t_init))


if __name__ == "__main__":
    main()
