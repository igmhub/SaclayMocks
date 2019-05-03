import time
import argparse
import numpy as np
import fitsio
import os

'''
This code is taken from LyaCoLoRe:
https://github.com/igmhub/LyaCoLoRe/blob/master/desi/generate_rnd_dla.py
'''

def generate_rnd(factor=3, infile=None, outfile=None):
    """
    Routine to generate a random catalog in 3D following
    certain N(z) distribution
    Args:
    ----
    factor: Size of the generated catalog (before masking)
    out_path: Output path
    """
    #Creating random that follows N(z)
    data = fitsio.read(infile)
    zvec = data['Z_DLA_RSD']
    ntot = int(len(data)*factor)

    # z_rnd = np.random.choice(zvec,size=ntot)+0.0025*np.random.normal(size=ntot)
    z_rnd = np.random.choice(zvec,size=ntot)
    qso_rnd = np.random.choice(np.arange(len(data)), size=ntot)

    #Get the RA, DEC, MOCKID and Z_QSO of the quasars
    ra_rnd = data['RA'][qso_rnd]
    dec_rnd = data['DEC'][qso_rnd]
    MOCKID_rnd = data['MOCKID'][qso_rnd]
    Z_QSO_RSD_rnd = data['Z_QSO_RSD'][qso_rnd]
    Z_QSO_NO_RSD_rnd = data['Z_QSO_NO_RSD'][qso_rnd]

    if outfile is not None:
        outfits = fitsio.FITS(outfile, 'rw', clobber=True)
        table = [ra_rnd,dec_rnd,z_rnd,Z_QSO_NO_RSD_rnd,Z_QSO_RSD_rnd,MOCKID_rnd]
        names=['RA','DEC','Z','Z_QSO_NO_RSD','Z_QSO_RSD','MOCKID']
        outfits.write(table, names=names, extname='DLACAT')
        outfits[-1].write_key("x", factor, comment="Nrand is x times Ndata")

    return None

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-infile', type = str, default = None, required = True,
                    help='path to DLA catalog fits file')
    parser.add_argument('-outfile', type = str, default = None, required = True,
                    help='out path to randoms DLA catalog fits file')
    parser.add_argument('-factor', type = float, default = 3, required = False,
                    help='factor to dermine the number of randoms based on number of data')
    args = parser.parse_args()
    t0 = time.time()
    print("Generating randoms DLA from:\n{}".format(args.infile))
    print("Will be written in {}".format(args.outfile))
    generate_rnd(factor=args.factor, infile=args.infile, outfile=args.outfile)
    print("Done. Took {} s".format(time.time() - t0))


if __name__ == "__main__":
    main()
