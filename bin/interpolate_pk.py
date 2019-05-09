#!/usr/bin/env python
# interpolates a table of P(k) in a file
# to compute sqrt(P(K)/Vcell) for each mode
# store result in fits file
import os
import numpy as np
from fitsio import FITS
import time
import argparse
from SaclayMocks import powerspectrum
from SaclayMocks import constant
from SaclayMocks import util
from multiprocessing import Pool
from memory_profiler import profile


#********************************************************************
def Power_Spectrum_ln(k, growth, bias, Vcell):
    filename = os.path.expandvars("$SACLAYMOCKS_BASE/etc/PlanckDR12.fits")
    myPln = np.float32(powerspectrum.P_ln(filename,G_times_bias=growth*bias).P(k))
    return np.float32( np.sqrt(myPln/Vcell) )

#********************************************************************
def Power_Spectrum(k, Vcell):
    filename = os.path.expandvars("$SACLAYMOCKS_BASE/etc/PlanckDR12.fits")
    myP = np.float32(powerspectrum.P_0(filename).P(k))
    return np.float32( np.sqrt(myP/Vcell) )

# @profile
def main():
    t0 = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument("-NX", type=int, help="number of pixels along x, default 256", default=256)
    parser.add_argument("-NY", type=int, help="number of pixels along y, default 256", default=256)
    parser.add_argument("-NZ", type=int, help="number of pixels along z, default 256", default=256)
    parser.add_argument("-pixel", type=float, help="pixel size (Mpc/h), default 2.19", default=2.19)
    parser.add_argument("-i", type=int, help="index of treated slice, default 0", default=0)
    parser.add_argument("-N", type=int, help="total number of slices, default 64", default=64)
    parser.add_argument("-outDir", help="Directory where the out file is saved")
    args = parser.parse_args()

    NX = args.NX
    NY = args.NY
    NZ = args.NZ
    Dcell = args.pixel
    iSlice = args.i
    NSlice = args.N
    outDir = args.outDir

    print("InterpolatePk : {}th slice over {}.".format(iSlice, NSlice))
    Vcell = np.float32(Dcell**3)
    PI = np.pi
    k_ny = PI / Dcell
    # ncpu = 2

    # .............................  make (kx,ky,kz) array
    #
    # extrapolation from 1024^3: 3072 x 3072 x 1152 requires 122 Go and 703 s
    # this costs 703 x 80 / 3600 = 15.6 h
    # should run 64 jobs, each computing 1/64 of the box
    # then merge it

    kx = np.fft.fftfreq(NX) * 2 * k_ny  # rfftfreq 0 -> 0.5   (NX)
    ky = np.fft.fftfreq(NY) * 2 * k_ny   # (NY)
    kz = np.fft.rfftfreq(NZ) * 2 * k_ny   # (NZ/2+1)        fftw has no rfftfreq function
        # is that correct ?   <=
    # print(kx.shape,ky.shape,kz.shape)

    # Select treated slice of kx
    kx = kx[iSlice*len(kx)//NSlice:(iSlice+1)*len(kx)//NSlice]  # NX/NSlice


    kz = np.float32(kz)
    # print(kx.shape,type(kx[0]))
    ky = np.float32(ky.reshape(-1, 1))        # (NY,1)
    kx = np.float32(kx.reshape(-1, 1, 1))        # (NX,1,1)
    k = np.sqrt(kx*kx + ky*ky + kz*kz)  # shape = (NX, NY, NZ/2+1)

    #print(k.nbytes/1024/1024," Mb for k")
    #print(type(kx[0,0,0]),type(k[0,0,0]))

    if (NY == NX and NZ == NX):
        Pfilename = outDir + "/P{}_{}_{}.fits".format(NX, iSlice, NSlice)
    else:
        Pfilename = outDir + "/P{}-{}-{}_{}_{}.fits".format(NX, NY, NZ, iSlice, NSlice)

    # constants
    omega_M_0 = constant.omega_M_0
    bias = constant.QSO_bias
    z_QSO_bias = constant.z_QSO_bias    # V1 produced with z_0 =2.5 !!!
    growth = util.fgrowth(z_QSO_bias, omega_M_0)

    # Open fits file
    fits = FITS(Pfilename, 'rw', clobber=True)

    # Check number of modes
    if len(k.ravel()) > 2**31:
        print("Number of modes too large : {} > 2**31 !\nExit.".format(len(k.ravel())))
        exit()

    # .............................   P_ln
    Pln = Power_Spectrum_ln(k, growth, bias, Vcell)
    hdict = {'Dcell': Dcell, 'NX': NX, 'NY': NY, 'NZ': NZ}
    fits.write(Pln, header=hdict)
    fits[-1].write_key("QSO_bias", bias, comment="QSO bias at z={}".format(z_QSO_bias))
    fits[-1].write_key("G", growth, comment="growth factor at z={}".format(z_QSO_bias))
    fits[-1].write_key("Om", omega_M_0, comment="Omega matter today")
    del Pln

    # .............................   P_0
    P_0 = Power_Spectrum(k, Vcell)
    hdict = {'Dcell': Dcell, 'NX': NX, 'NY': NY, 'NZ': NZ}
    fits.write(P_0, header=hdict)
    del P_0
    fits.close()
    print("produced ", Pfilename)
    print("Took {}s".format(time.time() - t0))


if __name__ == "__main__":
    main()
