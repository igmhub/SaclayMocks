#!/usr/bin/env python
import fitsio
import time
import argparse
import os
import numpy as np
import glob

from SaclayMocks import util

# if True:
def main():
    t_init = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument("-indir", help="Root directory of the various Chunks")
    parser.add_argument("-outfile", help="Path to out file")
    parser.add_argument("-nside", type=int, help="nside for healpix. Default 16", default=16)
    parser.add_argument("-nest", help="If True, healpix scheme is nest. Default True", default='True')

    args = parser.parse_args()
    nside = args.nside
    nest_option = util.str2bool(args.nest)
    files = glob.glob(args.indir+"/*/dla.fits")

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
    outfits = fitsio.FITS(args.outfile, 'rw', clobber=True)
    table = [ra, dec, z_qso, z_qso_rsd, z, z_rsd, nhi, mockid, dla_id, pixnum]
    names = ['RA', 'DEC', 'Z_QSO_NO_RSD', 'Z_QSO_RSD', 'Z_DLA_NO_RSD', 'Z_DLA_RSD', 'N_HI_DLA', 'MOCKID', 'DLAID', 'PIXNUM']
    outfits.write(table, names=names, extname='DLACAT')
    print("Done. {} s. Written {} DLAs".format(time.time()-t2, len(ra)))
    print("Took {} s".format(time.time()-t_init))


if __name__ == "__main__":
    main()
