#!/usr/bin/env python
from fitsio import FITS
import numpy as np
import argparse
import time
from memory_profiler import profile


# @profile
def main():
    t0 = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument("-inDir", help="Directory where files are merged")
    parser.add_argument("-outDir", help="Directory where the out file is saved")
    parser.add_argument("-N", type=int, help="total number of Pk files, default 64", default=64)
    parser.add_argument("-NX", type=int, help="number of pixels along x, default 256", default=256)
    parser.add_argument("-NY", type=int, help="number of pixels along y, default 256", default=256)
    parser.add_argument("-NZ", type=int, help="number of pixels along z, default 256", default=256)

    args = parser.parse_args()
    N_HDU = args.N
    inDir = args.inDir
    outDir = args.outDir
    NX = args.NX
    NY = args.NY
    NZ = args.NZ
    print("Merging {} Pk fits files".format(N_HDU))

    if (NX == NY and NX == NZ):
        outfits = FITS(outDir+'/P{}.fits'.format(NX), 'rw', clobber=True)
        infits = [FITS(inDir+'/P{}_{}_{}.fits'.format(NX, i, N_HDU)) for i in range(N_HDU)]
        head = infits[0][0].read_header()
        hdict = {'Dcell': head['Dcell'], 'NX': head['NX'], 'NY': head['NY'], 'NZ': head['NZ']}

        Pln = [f[0].read() for f in infits]
        Pln = np.concatenate(Pln, axis=0)
        outfits.write(Pln, header=hdict)
        del Pln

        P0 = [f[1].read() for f in infits]
        P0 = np.concatenate(P0, axis=0)
        outfits.write(P0, header=hdict)
        del P0

    else:
        outfits = FITS(outDir+'/P{}-{}-{}.fits'.format(NX, NY, NZ), 'rw', clobber=True)
        infits = [FITS(inDir+'/P{}-{}-{}_{}_{}.fits'.format(NX, NY, NZ, i, N_HDU)) for i in range(N_HDU)]
        head = infits[0][0].read_header()
        hdict = {'Dcell': head['Dcell'], 'NX': head['NX'], 'NY': head['NY'], 'NZ': head['NZ']}

        Pln = [f[0].read() for f in infits]
        Pln = np.concatenate(Pln, axis=0)
        outfits.write(Pln, header=hdict)
        del Pln

        P0 = [f[1].read() for f in infits]
        P0 = np.concatenate(P0, axis=0)
        outfits.write(P0, header=hdict)
        del P0

    print("Merged fits file written in {}".format(outDir))
    print("Took {}s".format(time.time()-t0))

if __name__ == "__main__":
    main()
