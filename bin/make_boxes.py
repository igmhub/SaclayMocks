#!/usr/bin/env python
# python DrawGRF -NX= -pixel=
# Draw GRF in k-space
# multiply by sqrt{ P(k) }
# iFFT to get box in r space
# save to fits file

import scipy as sp
import numpy as np
import fitsio
from fitsio import FITS,FITSHDR
import sys
import time
#from scipy import interpolate
#import matplotlib.pyplot as plt
from memory_profiler import profile
#import scipy.fftpack as fft
#import numpy.fft as fft
# import pyfftw.interfaces.numpy_fft as fft
#from subprocess import call
#from time import time
import pyfftw
import os
import argparse
import subprocess
from SaclayMocks import powerspectrum
from SaclayMocks import constant
from SaclayMocks import util
import gc


# z0 = constant.z0    # V1 produced with z_0 =2.5 !!!
# bias = constant.QSO_bias
# ncpu = 1
# use_pool = True


#********************************************************************
#@profile
def DrawGRF_boxk(NX,NY,NZ, ncpu, wisdomFile, box_null=False):
#        Draw GRF box in k space in numpy.fft format
# later we will directly draw boxk, with appropriate symetries
# and with var(delta_k)=NX*NY*NZ(see cahier simu FFT normalization)
  t0 = time.time()
  # box = np.float32(np.random.normal(size=[NX, NY, NZ]))
  box = np.zeros((NX,NY,NZ),dtype=float32)
  for iz in range(NZ):
    box[:,:,iz] =  np.float32(np.random.normal(size=[NX, NY]))
  t1 = time.time()
  print(box.nbytes/1024/1024, " Mbytes box drawn",  t1-t0, " s")
  print(box.dtype)
  boxk = np.zeros([NX,NY,NZ//2+1],dtype=np.complex64)
  myfft = pyfftw.FFTW(box,boxk,axes=(0,1,2),threads=ncpu, flags=('FFTW_DESTROY_INPUT',))
  myfft.execute()
  t2 =time.time()
  print("FFT", t2-t1, " s")

  # Next lines is a test of the FFT:
  # if the produced box is null, then try to save wisdom and redo FFT once again
  # using flag box_null = True : this will call DrawGRF_boxk a second time,
  # with the wisdom file previously saved, a try another FFT. If the box is still
  # null, then exit the code.
  if (not np.any(boxk)) or (np.any(np.isnan(boxk))):
    sp.save(wisdomFile, pyfftw.export_wisdom())
    if box_null:
      raise ValueError("/!\ boxk is null /!\ \n    Wisdom saved before exiting.")
    print("/!\ Boxk was null. Wisdom has been saved, trying again FFTW...")
    pyfftw.import_wisdom(sp.load(wisdomFile))
    del boxk
    boxk = DrawGRF_boxk(NX, NY, NZ, ncpu, wisdomFile, True)

  return boxk


#********************************************************************
# @profile
def FFTandStore(Dcell, nHDU, boxfilename, ncpu, wisdomFile, box_null=False):
#.............................  FFT
    global boxk   # use global variable boxk
    t0=time.time()
    t2 = time.time()
    NX = boxk.shape[0]
    NY = boxk.shape[1]
    NZ = boxk.shape[2]
    # box = np.zeros([NX,NY,2*(NZ-1)],dtype=np.float32)
    box = pyfftw.empty_aligned((NX, NY, 2*(NZ-1)), dtype='float32')
    pyfftw.FFTW(boxk,box,axes=(0,1,2),direction='FFTW_BACKWARD',threads=ncpu, flags=('FFTW_DESTROY_INPUT',)).execute()
    del boxk
    box /= NX*NY*2*(NZ-1)
    t3 = time.time()
    print("FFT done", t3-t2, "s")
    sigma = np.std(box)
    print("sigma = {}".format(sigma))

    # Next lines is a test of the FFT:
    # if the produced box is null, then try to save wisdom and redo FFT once again
    # using flag box_null = True : this will call FFTandStore a second time,
    # with the wisdom file previously saved, a try another FFT. If the box is still
    # null, then exit the code.
    if (not np.any(box)) or (np.any(np.isnan(box))):
      sp.save(wisdomFile, pyfftw.export_wisdom())
      if box_null:
        raise ValueError("/!\ box is null /!\ \n    Box name: {}\nWisdom saved before exiting.".format(boxfilename))
      print("/!\ Box was null. Wisdom has been saved, trying again FFTW...")
      box_null = True
    else:
      #...............................      write to fits file
      if box_null:
        print("Box is not null this time. Continuing...")
      for i in np.arange(0, nHDU):
        fits = FITS(boxfilename+'-{}.fits'.format(i),'rw',clobber=True)
        hdict = {'DX': Dcell, 'DY': Dcell, 'DZ':Dcell, 'NX':NX, 'NY':NY, 'NZ':(NZ-1)*2}
        fits.write(box[i*NX//nHDU:(i+1)*NX//nHDU],header=hdict)
        if i == 0:
          fits[0].write_key("sigma", np.float32(sigma), comment="std of the box")
          fits[0].write_key("seed", np.int32(seed), comment="seed used to generate randoms")
        fits.close()
      t4 = time.time()
      print(boxfilename, "written", t4-t3,"s")

    del box
    return


#********************************************************************
# if True:
# @profile
def main() :
  print("Starting DrawGRF...")
  h = constant.h
  omega_M_0 = constant.omega_M_0
  t_init = time.time()

  parser = argparse.ArgumentParser()
  parser.add_argument("-pixel", type=float, help="pixel size (Mpc/h), default 2.19", default=2.19)
  parser.add_argument("-NX", type=int, help="number of pixels along x, default 256", default=256)
  parser.add_argument("-NY", type=int, help="number of pixels along y, default = NX", default=-1)  # -1 means equal to NX
  parser.add_argument("-NZ", type=int, help="number of pixels along z, default = NX", default=-1)
  parser.add_argument("-nHDU", type=int, help="number of HDU box.fits, default 1", default=1)
  parser.add_argument("-ncpu", type=int, default=2)
  parser.add_argument("-PkDir", help="directory of Pk fits file")
  parser.add_argument("-seed", type=int, help="specify a seed", default=None)
  parser.add_argument("-rsd", type=str, help="If True, rsd are added, default True", default='True')
  parser.add_argument("-dgrowthfile", help="dD/dz file, default etc/dgrowth.fits", default=None)
  parser.add_argument("-outDir", help="directory where the box are saved")

  args = parser.parse_args()
  rsd = util.str2bool(args.rsd)
  Dcell = args.pixel
  NX = args.NX
  if (args.NY < 0):
    NY = NX
  else:
    NY = args.NY
  if (args.NZ < 0):
    NZ = NX
  else:
    NZ = args.NZ

  global seed
  seed = args.seed
  if seed is None:
    seed = np.random.randint(2**31 -1, size=1)[0]
    np.random.seed(seed)
    print("Seed has not been specified. Seed is set to {}".format(seed))
  else:
    np.random.seed(seed)
    print("Specified seed is {}".format(seed))

  nHDU = args.nHDU
  ncpu = args.ncpu
  PkDir = args.PkDir
  outDir = args.outDir

  PI = np.pi
  k_ny = PI / Dcell
  nCell  = NX * NY * NZ
  Vcell = np.float32(Dcell**3)
  volume = nCell * Vcell
  print("volume = ",volume)

  # Reading power spectra
  if (NY==NX and NZ==NX):
    Pfilename = PkDir+"/P"+str(NX)+".fits"
  else :
    Pfilename = PkDir+"/P"+str(NX)+"-"+str(NY)+"-"+str(NZ)+".fits"

  #...............................    get wisdom to save time on FFT
  wisdom_path = os.path.expandvars("$SACLAYMOCKS_BASE/etc/")
  if (NY==NX and NZ==NX):
    wisdomFile = wisdom_path+"wisdom."+str(NX)+"."+str(ncpu)+".npy"
  else :
    wisdomFile = wisdom_path+"wisdom."+str(NX)+"-"+str(NY)+"-"+str(NZ)+"."+str(ncpu)+".npy"

  if os.path.isfile(wisdomFile) :
    pyfftw.import_wisdom(sp.load(wisdomFile))
    save_wisdom = False
  else :
    print(wisdomFile," wisdomfile not found !!! set save_wisdom = true")
    save_wisdom = True

  #............................. Draw GRF in k space
  boxkfile = outDir + "/boxk.npy"
  global boxk
  if os.path.isfile(boxkfile):
    print("{} already exists ! Reading boxk.npy file to compute density and velocity boxes...".format(boxkfile))
    boxk = np.load(boxkfile)
    print("Done.")
    sigma = boxk.std()
    if sigma > 180000:
      print("Sigma of boxk is {} > 180000:".format(sigma))
      print("dividing boxk by P0(k) and saving...")
      t0 = time.time()
      boxk /= fitsio.read(Pfilename, ext='P0')
      boxk[0,0,0] = 0j
      np.save(boxkfile, boxk)
      print("Done. {} s".format(time.time()-t0))
  else:
    t0 = time.time()
    print(NX,NY,NZ)
    boxk = DrawGRF_boxk(NX,NY,NZ, ncpu, wisdomFile)
    t1 = time.time()
    np.save(boxkfile,boxk)
    np.save(outDir+"/seed_boxk.npy", seed)
    t2 = time.time()
    print("boxk produced and saved:",t1-t0,t2-t1," s ")

  #............................. multiply by sqrt(P/Vcell), FFT and store
  print("Computing delta boxes...")
  # First lognormal
  boxfile = outDir+'/boxln_1'
  command = "ls -l {}* | wc -l".format(boxfile)
  if nHDU == int(subprocess.run(command, shell=True, capture_output=True).stdout.decode('UTF-8')[:-1]):
    print("{} files already exist ! Skiping this step.".format(boxfile))
  else:
    boxk *= fitsio.read(Pfilename, ext='Pln1')
    box_null = FFTandStore(Dcell, nHDU, boxfile, ncpu, wisdomFile)
    # This box_null thing is in case the FFT went bad with the wisdom file
    if box_null:
      print("Starting again FFTandStore...")
      print("Loading wisdom {}".format(wisdomFile))
      pyfftw.import_wisdom(sp.load(wisdomFile))
      boxk = np.load(boxkfile)
      boxk *= fitsio.read(Pfilename, ext=0)
      FFTandStore(Dcell, nHDU, boxfile, ncpu, wisdomFile, box_null)

  # Second lognormal
  boxfile = outDir+'/boxln_2'
  command = "ls -l {}* | wc -l".format(boxfile)
  if nHDU == int(subprocess.run(command, shell=True, capture_output=True).stdout.decode('UTF-8')[:-1]):
    print("{} files already exist ! Skiping this step.".format(boxfile))
  else:
    boxk = np.load(boxkfile)
    boxk *= fitsio.read(Pfilename, ext='Pln2')
    FFTandStore(Dcell, nHDU, boxfile, ncpu, wisdomFile)
  # Third lognormal
  boxfile = outDir+'/boxln_3'
  command = "ls -l {}* | wc -l".format(boxfile)
  if nHDU == int(subprocess.run(command, shell=True, capture_output=True).stdout.decode('UTF-8')[:-1]):
    print("{} files already exist ! Skiping this step.".format(boxfile))
  else:
    boxk = np.load(boxkfile)
    boxk *= fitsio.read(Pfilename, ext='Pln3')
    FFTandStore(Dcell, nHDU, boxfile, ncpu, wisdomFile)
  # Density field
  nHDU_bis = NX   # we want 1 HDU per ix
  boxfile = outDir+'/box'
  command = "ls -l {}-* | wc -l".format(boxfile)
  if nHDU_bis == int(subprocess.run(command, shell=True, capture_output=True).stdout.decode('UTF-8')[:-1]):
    print("{} files already exist ! Skiping this step.".format(boxfile))
    print("Reading and multiplying boxk by P0(k), and saving...")
    t0 = time.time()
    boxk = np.load(boxkfile)
    boxk *= fitsio.read(Pfilename, ext='P0')
    np.save(boxkfile, boxk)
    print("Done. {} s".format(time.time() - t0))
  else:
    boxk=np.load(boxkfile)
    boxk *= fitsio.read(Pfilename, ext='P0')
    np.save(boxkfile, boxk)
    FFTandStore(Dcell, nHDU_bis, boxfile, ncpu, wisdomFile)

  if rsd:
    # ............................ Compute eta
    H0 = constant.H0
    Om = constant.omega_M_0
    print("Computing eta boxes:")
    kx = np.fft.fftfreq(NX) * 2 * k_ny  # rfftfreq 0 -> 0.5   (NX)
    ky = np.fft.fftfreq(NY) * 2 * k_ny   # (NY)
    kz = np.fft.rfftfreq(NZ) * 2 * k_ny   # (NZ/2+1)        fftw has no rfftfreq function
    kz = np.float32(kz)
    ky = np.float32(ky.reshape(-1, 1))        # (NY,1)
    kx = np.float32(kx.reshape(-1, 1, 1))        # (NX,1,1)
    kk = kx*kx + ky*ky + kz*kz  # shape = (NX, NY, NZ/2+1)
    kk[0, 0, 0] = 1  # avoid dividing by 0
    if args.dgrowthfile is None:
        filename = os.path.expandvars("$SACLAYMOCKS_BASE/etc/dgrowth.fits")
    if Om != fitsio.read_header(filename, ext=1)['OM']:
      raise ValueError("Omega_M_0 in SaclayMocks.constant ({}) != OM in {}".format(Om,
                            fitsio.read_header(filename, ext=1)['OM']))

    dgrowth0 = fitsio.read(filename, ext=1)['dD/dz'][0]  # value for z=0

    # etak_xx
    boxfile = outDir+'/eta_xx'
    command = "ls -l {}* | wc -l".format(boxfile)
    if nHDU_bis == int(subprocess.run(command, shell=True, capture_output=True).stdout.decode('UTF-8')[:-1]):
      print("{} files already exist ! Skiping this step.".format(boxfile))
    else:
      print("eta_xx...")
      t0 = time.time()
      boxk = np.load(boxkfile)
      boxk *= kx*kx / kk
      FFTandStore(Dcell, nHDU_bis, boxfile, ncpu, wisdomFile)
      print("Done. {} s".format(time.time() -t0))

    # etak_yy
    boxfile = outDir+'/eta_yy'
    command = "ls -l {}* | wc -l".format(boxfile)
    if nHDU_bis == int(subprocess.run(command, shell=True, capture_output=True).stdout.decode('UTF-8')[:-1]):
      print("{} files already exist ! Skiping this step.".format(boxfile))
    else:
      print("eta_yy...")
      t0 = time.time()
      boxk = np.load(boxkfile)
      boxk *= ky*ky / kk
      FFTandStore(Dcell, nHDU_bis, boxfile, ncpu, wisdomFile)
      print("Done. {} s".format(time.time() -t0))

    # etak_zz
    boxfile = outDir+'/eta_zz'
    command = "ls -l {}* | wc -l".format(boxfile)
    if nHDU_bis == int(subprocess.run(command, shell=True, capture_output=True).stdout.decode('UTF-8')[:-1]):
      print("{} files already exist ! Skiping this step.".format(boxfile))
    else:
      print("eta_zz...")
      t0 = time.time()
      boxk = np.load(boxkfile)
      boxk *= kz*kz / kk
      FFTandStore(Dcell, nHDU_bis, boxfile, ncpu, wisdomFile)
      print("Done. {} s".format(time.time() -t0))

    # etak_xy
    boxfile = outDir+'/eta_xy'
    command = "ls -l {}* | wc -l".format(boxfile)
    if nHDU_bis == int(subprocess.run(command, shell=True, capture_output=True).stdout.decode('UTF-8')[:-1]):
      print("{} files already exist ! Skiping this step.".format(boxfile))
    else:
      print("eta_xy...")
      t0 = time.time()
      boxk = np.load(boxkfile)
      boxk *= kx*ky / kk
      FFTandStore(Dcell, nHDU_bis, boxfile, ncpu, wisdomFile)
      print("Done. {} s".format(time.time() -t0))

    # etak_xz
    boxfile = outDir+'/eta_xz'
    command = "ls -l {}* | wc -l".format(boxfile)
    if nHDU_bis == int(subprocess.run(command, shell=True, capture_output=True).stdout.decode('UTF-8')[:-1]):
      print("{} files already exist ! Skiping this step.".format(boxfile))
    else:
      print("eta_xz...")
      t0 = time.time()
      boxk = np.load(boxkfile)
      boxk *= kx*kz / kk
      FFTandStore(Dcell, nHDU_bis, boxfile, ncpu, wisdomFile)
      print("Done. {} s".format(time.time() -t0))

    # etak_yz
    boxfile = outDir+'/eta_yz'
    command = "ls -l {}* | wc -l".format(boxfile)
    if nHDU_bis == int(subprocess.run(command, shell=True, capture_output=True).stdout.decode('UTF-8')[:-1]):
      print("{} files already exist ! Skiping this step.".format(boxfile))
    else:
      print("eta_yz...")
      t0 = time.time()
      boxk = np.load(boxkfile)
      boxk *= ky*kz / kk
      FFTandStore(Dcell, nHDU_bis, boxfile, ncpu, wisdomFile)
      print("Done. {} s".format(time.time() -t0))

    print("Computing velocity boxes:")
    # vx
    boxfile = outDir+'/vx'
    command = "ls -l {}* | wc -l".format(boxfile)
    if nHDU == int(subprocess.run(command, shell=True, capture_output=True).stdout.decode('UTF-8')[:-1]):
      print("{} files already exist ! Skiping this step.".format(boxfile))
    else:
      print("vx...")
      t0 = time.time()
      boxk = np.load(boxkfile)
      boxk *= -1j*kx / kk * H0 * dgrowth0
      FFTandStore(Dcell, nHDU, boxfile, ncpu, wisdomFile)
      print("Done. {} s".format(time.time() -t0))

    # vy
    boxfile = outDir+'/vy'
    command = "ls -l {}* | wc -l".format(boxfile)
    if nHDU == int(subprocess.run(command, shell=True, capture_output=True).stdout.decode('UTF-8')[:-1]):
      print("{} files already exist ! Skiping this step.".format(boxfile))
    else:
      print("vy...")
      t0 = time.time()
      boxk = np.load(boxkfile)
      boxk *= -1j*ky / kk * H0 * dgrowth0
      FFTandStore(Dcell, nHDU, boxfile, ncpu, wisdomFile)
      print("Done. {} s".format(time.time() -t0))

    # vz
    boxfile = outDir+'/vz'
    command = "ls -l {}* | wc -l".format(boxfile)
    if nHDU == int(subprocess.run(command, shell=True, capture_output=True).stdout.decode('UTF-8')[:-1]):
      print("{} files already exist ! Skiping this step.".format(boxfile))
    else:
      print("vz...")
      t0 = time.time()
      boxk = np.load(boxkfile)
      boxk *= -1j*kz / kk * H0 * dgrowth0
      FFTandStore(Dcell, nHDU, boxfile, ncpu, wisdomFile)
      print("Done. {} s".format(time.time() -t0))

  print("NX=", NX,"nCPU=", ncpu)  #, "use_pool=",  use_pool
  if (save_wisdom):
    sp.save(wisdomFile, pyfftw.export_wisdom())

  print("Took {}s".format(time.time()-t_init))


if __name__ == "__main__":
    main()
