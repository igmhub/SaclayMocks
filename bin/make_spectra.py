#!/usr/bin/env python
#  MakeSpectra.py 
#	reads QSO list
# 	reads non parallel l.o.s. through boxes
#   add small scale power
#	applies Gunn Peterson to make Lya forest
from __future__ import division, print_function
from LyaMocks import box
from LyaMocks import constant 
from LyaMocks import powerspectrum
from LyaMocks import util
import fitsio
from fitsio import FITS
import sys
import scipy as sp
import numpy as np
import cosmolopy.distance as dist
import os
import time
from numba import jit
import pyfftw.interfaces.numpy_fft as fft
import pyfftw
import json
import argparse
import scipy.stats as stats
import matplotlib.pyplot as plt 
import cosmolopy.perturbation as pert
from scipy import interpolate, integrate
from memory_profiler import profile


# @profile
def main():
# if True:
    #  .................... hardcoded param
    PI = np.pi
    plotPkMis = False
    #*************************************************************
    @jit       #    @jit improves from 73 to 32 ms
    def ComputeWeight(X,Y,Z,grid,cells,onesline, sig2) : # LX,LY,LZ,DX,DY,DZ
        ix = int((X +LX/2)/DX)  # -LX/2 < X < -LX/2 + LX/Nslice so 0 < ix < NX/Nslice
        iy = int((Y +LY/2)/DY)  # -LY/2 < Y < LY/2 so 0 < iy < NY
        iz = int((Z +LZ/2 -R0)/DZ)  
        ixyz = sp.array([ix,iy,iz])  	# (3,)
        cell_center = sp.array([ix*DX-LX/2,iy*DY-LY/2,iz*DZ-LZ/2+R0])
        XYZ = sp.array([X,Y,Z]) 
        
        lgrid = grid + XYZ	# grid around XYZ
        lcells = cells + ixyz
        cell_center_line = cell_center * onesline	
        # (3,) * (343,1) => (343,3) 	343= (2*dmax+1)**3 
        #    if (icell==0): print  
        
        xx = (cell_center_line-lgrid)**2
        Delta_r2 = xx[:,0]+xx[:,1]+xx[:,2]    # Delta_r2 = ((cell_center_line-lgrid)**2).sum(axis=1)
        weight = sp.exp(-Delta_r2 / sig2)
        return lcells, weight
    
    #*************************************************************
    #@jit   # @jit degrades from 22 to 27 ms
    def computeRho(myrho,weight) :
        sumweight = weight.sum(axis=0)
        sumrho = (weight*myrho).sum(axis=0)
        return sumrho / sumweight
    
    #*************************************************************
    #   @jit + python -m cProfile fails 
    #   @jit degrades from 80 t0 94 for the full treatment of a QSO
    #@jit      
    def ReadSpec(Xvec, XvecSlice, Yvec, Zvec, grid, cells, onesline, imin=0, imax=sys.maxint):
        # LX,LY,LZ,DX,DY,DZ        
        # (Xvec, Yvec, Zvec) is the vector along line of sight,
        # XvecSlice is in [-LX/2, -LX/2 + LX/NSlice]
        # cells is the list of indices used for G.S. around (0,0,0),
        # and grid its value in Mpc/h,
        # imin imax are the indices delimiting the lya forest 
        if rsd:
            eta_par = sp.zeros_like(XvecSlice)  # so that exp(-a(exp(b*g) + taubar_a*eta_par)) = 1
            if dla:
                vpar = sp.zeros_like(XvecSlice)

        spectrum = -1000000 * sp.ones_like(XvecSlice)
        imax = np.minimum(imax,XvecSlice.size)
        sig2=2*DX*DX
        for icell in range(imin,imax):
            X = XvecSlice[icell]
            Xtrue = Xvec[icell]
            Y = Yvec[icell]
            Z = Zvec[icell]
            lcells, weight = ComputeWeight(X,Y,Z,grid,cells,onesline, sig2)
            myrho = fullrho[lcells[:,0],lcells[:,1],lcells[:,2]]
            spectrum[icell] = computeRho(myrho,weight)
            if rsd:
                RR = Xtrue**2+Y**2+Z**2
                myeta_xx = eta_xx[lcells[:,0],lcells[:,1],lcells[:,2]]
                myeta_yy = eta_yy[lcells[:,0],lcells[:,1],lcells[:,2]]
                myeta_zz = eta_zz[lcells[:,0],lcells[:,1],lcells[:,2]]
                myeta_xy = eta_xy[lcells[:,0],lcells[:,1],lcells[:,2]]
                myeta_xz = eta_xz[lcells[:,0],lcells[:,1],lcells[:,2]]
                myeta_yz = eta_yz[lcells[:,0],lcells[:,1],lcells[:,2]]
                myeta_xx = computeRho(myeta_xx, weight)
                myeta_yy = computeRho(myeta_yy, weight)
                myeta_zz = computeRho(myeta_zz, weight)
                myeta_xy = computeRho(myeta_xy, weight)
                myeta_xz = computeRho(myeta_xz, weight)
                myeta_yz = computeRho(myeta_yz, weight)
                eta_par[icell] = (Xtrue*myeta_xx*Xtrue + Y*myeta_yy*Y
                                  + Z*myeta_zz*Z + 2*Xtrue*myeta_xy*Y
                                 + 2*Xtrue*myeta_xz*Z + 2*Y*myeta_yz*Z) / RR
                if dla:
                    vx = velo_x[lcells[:,0],lcells[:,1], lcells[:,2]]
                    vy = velo_y[lcells[:,0],lcells[:,1], lcells[:,2]]
                    vz = velo_z[lcells[:,0],lcells[:,1], lcells[:,2]]
                    vx = computeRho(vx, weight)
                    vy = computeRho(vy, weight)
                    vz = computeRho(vz, weight)
                    vpar[icell] = (vx*Xtrue + vy*Y + vz*Z)/np.sqrt(RR)

        if rsd:
            if dla:
                return spectrum, eta_par, vpar
            else:
                return spectrum, eta_par
        else:
            return spectrum


    #************************************************************* main
    #..................  PARAMETERS 
    t_init = time.time()
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-dmax", type=int, help="compute rho over +- dmax cells, default 3", default=3)
    parser.add_argument("-pixel", type=float, help="spectrum pixel in Mpc/h, default 0.2", default=0.2)
    parser.add_argument("-zmin", type=float, help="min redshift. Default is 1.3", default=1.3)
    parser.add_argument("-zmax", type=float, help="max redshift. Default is 3.6", default=3.6)
    parser.add_argument("-QSOfile", help="QSO fits file")
    parser.add_argument("-boxdir", help="path to box fits files")
    parser.add_argument("-outDir", help="dir for spectra fits file")
    parser.add_argument("-i", type=int, help="index of treated slice")
    parser.add_argument("-N", type=int, help="total number of slices")
    parser.add_argument("-NQSOfile", type=int, help="number of QSO files, default = N", default=-1)
    parser.add_argument("-NQSO", type=int, help="cut at QNSO, default = -1, no cut", default=-1)
    parser.add_argument("-rsd", help="If True, rsd are added, default True", default='True')
    parser.add_argument("-dla", help="If True, store delta and growth skewers, default False", default='False')
    parser.add_argument("-dgrowthfile", help="dD/dz file, default data/dgrowth.fits", default="data/dgrowth.fits")
    args = parser.parse_args()

    iSlice = args.i
    NSlice = args.N
    dmax = args.dmax
    nqsomax = args.NQSO
    DeltaR = args.pixel
    boxdir = args.boxdir
    rsd = util.str2bool(args.rsd)
    dla = util.str2bool(args.dla)

    Om = constant.omega_M_0
    if Om == fitsio.read_header(args.dgrowthfile, ext=1)['OM']:
        Dgrowth = util.InterpFitsTable(args.dgrowthfile, 'Z', 'dD/dz')
    else:
        raise ValueError("Omega_M_0 in LyaMocks.constant ({}) != OM in {}".format(Om,
                            fitsio.read_header(args.dgrowthfile, ext=1)['OM']))
    dgrowth0 = Dgrowth.interp(0)

    print("Begining of MakeSpectra - {}".format(iSlice))
    h = constant.h
    Om = constant.omega_M_0
    lya = constant.lya
    lylimit = constant.lylimit
    lambda_min = constant.lambda_min
    lyb = constant.lyb
    Om = constant.omega_M_0
    OL = constant.omega_lambda_0
    Ok = constant.omega_k_0
    z0 = constant.z0

    #...............................   read box file
    print("Reading delta box...")
    t0 = time.time()
    head = fitsio.read_header(boxdir+"/box-0.fits", ext=0)
    DX = head["DX"]
    DY = head["DY"]
    DZ = head["DZ"]
    NZ = head["NAXIS1"]
    NY = head["NAXIS2"]
    NX = head["NAXIS3"]
    nHDU = head["NX"]
    LX = DX * nHDU
    LY = DY * NY
    LZ = DZ * NZ
    
    if (NX != 1):
        print ("NX=", NX, " != 1  => abort !")
        exit(1)
        
    if (iSlice >= NSlice):
        print ('iSlice=',iSlice,">= NSlice =" , NSlice , "=> abort !")
        exit(1)
        
    if (nHDU%NSlice != 0):
        print (NSlice, " slices not divider of nHDU:", nHDU, "=> abort !")
        exit(1)

    iXmin = np.maximum((iSlice * nHDU)//NSlice -dmax,0)
    iXmax = np.minimum(((iSlice+1) * nHDU)//NSlice +dmax,nHDU)
    fullrho = []
    for iX in np.arange(iXmin,iXmax):
        fullrho.append(fitsio.read(boxdir+"/box-{}.fits".format(iX), ext=0))

    fullrho = np.concatenate(fullrho)
    print("Done. {} s".format(time.time() - t0))

    if rsd:
        print("Reading eta boxes...")
        t1 = time.time()
        eta_xx = []
        eta_yy = []
        eta_zz = []
        eta_xy = []
        eta_xz = []
        eta_yz = []
        for iX in np.arange(iXmin, iXmax):
            eta_xx.append(fitsio.read(boxdir+"/eta_xx-{}.fits".format(iX), ext=0))
            eta_yy.append(fitsio.read(boxdir+"/eta_yy-{}.fits".format(iX), ext=0))
            eta_zz.append(fitsio.read(boxdir+"/eta_zz-{}.fits".format(iX), ext=0))
            eta_xy.append(fitsio.read(boxdir+"/eta_xy-{}.fits".format(iX), ext=0))
            eta_xz.append(fitsio.read(boxdir+"/eta_xz-{}.fits".format(iX), ext=0))
            eta_yz.append(fitsio.read(boxdir+"/eta_yz-{}.fits".format(iX), ext=0))

        print("Concatenating... {} s".format(time.time()-t1))
        eta_xx = np.concatenate(eta_xx)
        eta_yy = np.concatenate(eta_yy)
        eta_zz = np.concatenate(eta_zz)
        eta_xy = np.concatenate(eta_xy)
        eta_xz = np.concatenate(eta_xz)
        eta_yz = np.concatenate(eta_yz)
        print("Done. {} s".format(time.time()-t1))

        if dla:
            print("Reading velocity boxes...")
            t1 = time.time()
            velo_x = []
            velo_y = []
            velo_z = []
            i = int(iXmin*NSlice/nHDU)
            j = np.minimum(int(iXmax*NSlice/nHDU + 1), NSlice)
            for iX in np.arange(i , j):
                velo_x.append(fitsio.read(boxdir+"/vx-{}.fits".format(iX), ext=0))
                velo_y.append(fitsio.read(boxdir+"/vy-{}.fits".format(iX), ext=0))
                velo_z.append(fitsio.read(boxdir+"/vz-{}.fits".format(iX), ext=0))

            print("Concatenating... {} s".format(time.time()-t1))
            ii = int(iXmin % (nHDU / NSlice))
            jj = int(((j-i-1)*nHDU/NSlice + dmax))
            if iSlice == NSlice-1: jj=int((j-i)*nHDU/NSlice)
            velo_x = np.concatenate(velo_x)[ii:jj]
            velo_y = np.concatenate(velo_y)[ii:jj]
            velo_z = np.concatenate(velo_z)[ii:jj]
            print("Done. {} s".format(time.time()-t1))

    # print rho cells to check <==
    xSlicemin = LX * iSlice / NSlice - LX/2
    xSlicemax = LX * (iSlice+1) / NSlice - LX/2
    print("Box {} - {} - {} with LX = {}, LY = {}, LZ = {}".format(nHDU, NY, NZ, LX, LY, LZ))
    print ("slice #",iSlice,"of box: ", xSlicemin," < x < ",xSlicemax,)
    print (",  due to dmax=",dmax,"requires", iXmin,"<= ix <",iXmax,fullrho.shape,"   ")

    #.................................................  	set the box at z0
    # http://roban.github.io/CosmoloPy/docAPI/cosmolopy.distance-module.html
    # cosmolopy returns distance in Mpc not Mpc/h
    cosmo_fid = {'omega_M_0':Om, 'omega_lambda_0':OL, 'omega_k_0':Ok, 'h':h}
    R_of_z, z_of_R = dist.quick_distance_function(dist.comoving_distance, return_inverse=True, **cosmo_fid)
    h_of_z = interpolate.interp1d(np.arange(0,5,0.01), dist.hubble_z(np.arange(0,5,0.01), **cosmo_fid))
    R0 = h * R_of_z(z0)
    Rmin,Rmax,tanx_max,tany_max = box.box_limit(LX,LY,LZ,R0,dmax*DX)
    z_low = args.zmin
    Rmin = h * R_of_z(z_low)
    z_high = args.zmax
    Rmax = h * R_of_z(z_high)
    
    npixeltot = int((Rmax - Rmin) /DeltaR +0.5) 
    R_vec = Rmin + np.arange(npixeltot) * DeltaR
    lambda_vec = lya * (1 + z_of_R(R_vec/h))
    cut = (lambda_vec > lambda_min)  # cut pixel bellow 3530 A
    R_vec = R_vec[cut]
    lambda_vec = lambda_vec[cut]
    print ("z0 =",z0, "=> R0=",R0,"and", Rmin,"<R<", Rmax,"thetax,y <",tanx_max, tany_max)
    print ("   ",z_low,"< z <",z_high,";   ",(1+z_low)*lya,"< lambda <",(1+z_high)*lya,"  =>",npixeltot,"pixels")
    print("+ cut spectro : {} < lambda < {} => {} pixels".format(
        np.min(lambda_vec), np.max(lambda_vec), len(lambda_vec)))
    npixeltot = len(lambda_vec)

    #...................................................	read QSO file
    if args.NQSOfile < 0:
        NQSOfile = NSlice
    else:
        NQSOfile = args.NQSOfile
    QSOfilename = args.QSOfile
    if (iSlice >= NSlice//2):
        # ifile0 = iSlice * NQSOfile // NSlice  # tbc when ca tombe pas rond <==
        ifile0 = NSlice//2
        # not optimal, we are trying irrelevant QSOs, but conservative
        ifile1 = NQSOfile   # could be reduced <==
        #Xmin = DX * iSlice -LX/2 # left edge of the slice
        #ZoverXmin = (R0-LZ/2) / (DX * (iSlice+1) -LX/2)   
        # Z/X > Z/X of lower right edge of the slice
        tanx_slice_max = (DX * (nHDU/NSlice) * (iSlice+1) -LX/2) / (R0-LZ/2)  
        # lower right edge of the slice
        #print ((DX * nHDU/NSlice) * (iSlice+1), LX/2, R0-LZ/2, tanx_slice_max)
        #ra_max = np.atan(tanx_slice_max)
        #   ra and dec properly computed, here and in DrawQSO ?   <==
    else :
        ifile0 = 0 # could be increased <==
        # ifile1 = (iSlice+1) * NQSOfile // NSlice # why +1 ?  <==
        ifile1 = NSlice//2   
        # not optimal, we are trying irrelevant QSOs, but conservative
        #Xmax = DX * (iSlice+1) * nHDU/NSlice - LX/2 # right edge of the slice
        #ZoverXmax = (R0-LZ/2) / (DX * iSlice -LX/2)
        # |Z/X| > |Z/X| of lower left edge of the slice, but X < 0
        tanx_slice_min = (DX * (nHDU/NSlice) * iSlice -LX/2) / (R0-LZ/2)  
        # lower right edge of the slice
        tanx_slice_max = np.abs(tanx_slice_min)
        #print ((DX * nHDU/NSlice) * (iSlice), LX/2, R0-LZ/2, tanx_slice_max)
        #ra_min = np.atan(tan_min)
    print ("use QSO files:",ifile0, "to", ifile1-1)
    
    qsos = []
    first = True
    for ifile in np.arange(ifile0,ifile1) :
        QSOfilename = args.QSOfile + str(ifile)+'-'+str(NQSOfile)+'.fits'
        try:
            fits = FITS(QSOfilename,'r')
        except IOError:
            print("*Warning* Fits file {} cannot be read.".format(QSOfilename))
            continue

        if first: 
            qsos.append(fits[1].read())
            head = fits[1].read_header()
            ra0 = head["RA0"]
            dec0 = head["DEC0"]
            first = False
        else : 
            qsos.append(fits[1].read())

        fits.close()
    qsos = np.concatenate(qsos)
    if len(qsos) == 0:
        print("No QSO read. ==> Exit.")
        sys.exit(0)
    print (len(qsos),"QSO read")

    #............................. draw non parallel l.o.s. between Rmin and R_QSO    
    #   for dmax=3 produce array([[-3, -3, -3], [-3, -3, -2], [-3, -3, -1],
    #       ...,        [ 3,  3,  1], [ 3,  3,  2], [ 3,  3,  3]])
    xx = 2 * dmax + 1
    cells = sp.indices((xx,xx,xx)).reshape(3,-1) - dmax
    grid = sp.array([DX*cells[0],DY*cells[1],DZ*cells[2]])
    cells = cells.T
    grid = grid.T
    onesline = sp.ones((xx**3,1))	# (343,1)
    
    wisdomFile = "wisdom1D" # should depend on #CPU  <==
    if os.path.isfile(wisdomFile):
        g = open(wisdomFile, 'r')
        wisd=json.load(g)
        pyfftw.import_wisdom(wisd)
        g.close()
        save_wisdom = False
    else :
        print ("Warning ***** no wisdom file, set save_wisdom = true ******")
        save_wisdom = True

    iqso=0
    pixtot=0
    maxsize = 0
    lambda_list = []
    delta_l_list = []
    eta_list = []
    velo_list = []
    redshift_list = []
    ra_list = []
    dec_list = []
    zQSO_norsd_list = []
    zQSO_rsd_list = []
    QSOhdu_list = []
    QSOid_list = []
    plate_list = []
    mjd_list = []
    fiber_list = []
    pmf_list = []
    t0 = time.time()
    for qso in (qsos) :		#............................  loop over qso
        if ((nqsomax > 0) & (iqso >= nqsomax)) :break 
        ra = qso['RA']
        dec = qso['DEC']
        if rsd:
            zQSO = qso['Z_QSO_RSD']
        else:
            zQSO = qso['Z_QSO_NO_RSD']

        zQSO_norsd = qso['Z_QSO_NO_RSD']
        zQSO_rsd = qso['Z_QSO_RSD']
        QSOid = qso['THING_ID']
        QSOhdu = qso['HDU']
        plate = qso['PLATE']
        mjd = qso['MJD']
        fiber = qso['FIBERID']
        pmf = qso['PMF']
        R_QSO = h * R_of_z(zQSO)
        X_QSO, Y_QSO, Z_QSO = box.ComputeXYZ2(np.radians(ra),np.radians(dec),
                                    R_QSO,np.radians(ra0),np.radians(dec0))
        tanx = X_QSO / Z_QSO  
        tany = Y_QSO / Z_QSO     
        if (np.abs(tanx) > tanx_slice_max ): 
            continue 
        # this cut is not mandatory, if missed results in npixel = 0
        if ((zQSO<z_low) | (zQSO>z_high)):
            continue
    
        # Next lines can be optimized by regrouping all cuts into one
        sinx = np.sign(tanx) * np.sqrt(tanx*tanx/(1+tanx*tanx))
        siny = np.sign(tany) * np.sqrt(tany*tany/(1+tany*tany))
        cut = (R_vec * X_QSO/R_QSO > xSlicemin) # Xvec > xSlicemin
        Rvec=R_vec[cut]
        mylambda = lambda_vec[cut]
        cut = (Rvec * X_QSO/R_QSO <= xSlicemax) #  Xvec < xSlicemax
        Rvec=Rvec[cut]
        mylambda = mylambda[cut]
        redshift = z_of_R(Rvec/h)
        Xvec = Rvec * X_QSO/R_QSO
        Yvec = Rvec * Y_QSO/R_QSO
        Zvec = Rvec * Z_QSO/R_QSO
        npixel = len(Rvec)
        if (npixel < 1) :
            continue

        XvecSlice = Xvec-LX * iSlice / NSlice  # Xvec within the considered slice
        # -LX/2 < Xvec < LX/2, while the read box is only LX/Nslice wide, 
        # so XvecSlice must start at zero in order to get correct indices
        # (cf ComputeWeight())
        if (iSlice > 0) : XvecSlice += DX * dmax # Xvec within the considered slice
        # Actually the read box includes a DX*dmax patch before the slice itself
        # so xVecSlice must start at DX*dmax
        # this is not true for slice#0, which should really start at zero.

        # select indices to compute spectra
        xyz = np.where(mylambda > lylimit*(1+zQSO))[0]
        if (len(xyz) == 0) : imin=len(mylambda)
        else : imin = xyz[0]
        xyz = np.where(mylambda < lya*(1+zQSO))[0]
        if (len(xyz) == 0) : imax=0
        else : imax = xyz[-1] + 1
        iqso += 1

        # Read boxes along l.o.s and apply smoothing
        if rsd:
            if dla:
                try:
                    delta_l, eta_par, velo_par = ReadSpec(Xvec, XvecSlice, Yvec, Zvec, grid, cells, onesline, imin=imin, imax=imax)
                except:
                    print("***WARNING ReadSpec:\n    ID {}***".format(QSOid))
                    continue
            else:
                try:
                    delta_l, eta_par = ReadSpec(Xvec, XvecSlice, Yvec, Zvec, grid, cells, onesline, imin=imin, imax=imax)
                except:
                    print("***WARNING ReadSpec:\n    ID {}***".format(QSOid))
                    continue
        else:
            delta_l = ReadSpec(Xvec, XvecSlice, Yvec, Zvec, grid, cells, onesline, imin=imin, imax=imax)

        # LX,LY,LZ,DX,DY,DZ are hidden parameters,
        # as well as fullrho and eta_x, eta_y eta_z

        lrf =  mylambda/(1+zQSO)
        cut = ((lrf<lya) & (lrf>lylimit))
        pixtot += delta_l[cut].size
        if maxsize < delta_l.size:
            maxsize = delta_l.size

        # Append to list:
        lambda_list.append(np.float32(mylambda))
        delta_l_list.append(np.float32(delta_l))
        if rsd:
            eta_list.append(np.float32(eta_par))
            if dla:
                # extra (1+z) factor for dz = (1+z)*v/c
                # velo_par *= (1+redshift)**2 * Dgrowth.interp(redshift) / dgrowth0
                velo_par *= (1+redshift) * h_of_z(redshift) / h_of_z(0) * Dgrowth.interp(redshift) / dgrowth0
                velo_list.append(np.float32(velo_par))
        redshift_list.append(np.float32(redshift))
        ra_list.append(ra)
        dec_list.append(dec)
        zQSO_norsd_list.append(zQSO_norsd)
        zQSO_rsd_list.append(zQSO_rsd)
        QSOhdu_list.append(QSOhdu)
        QSOid_list.append(QSOid)
        plate_list.append(plate)
        mjd_list.append(mjd)
        fiber_list.append(fiber)
        pmf_list.append(pmf)

    # end of loop on QSO: write spectra to fits files
    print("End of loop: {}s".format(time.time()-t0))
    print("Writting...")
    t1 = time.time()
    names = ["RA", "DEC", "Z_noRSD", "Z", "HDU", "THING_ID", "PLATE", "MJD", "FIBERID", "PMF"]
    hlist = [{'name':"z0", 'value':z0, 'comment':"redshift of box center"},
                  {'name':"pixel", 'value':DeltaR},
                  {'name':"Npixel", 'value':npixeltot},
                  {'name':"NX", 'value':NX},
                  {'name':"dmax", 'value':dmax},
                  {'name':"ra0", 'value':ra0, 'comment':"right ascension of box center"},
                  {'name':"dec0", 'value':dec0, 'comment':"declination of box center"}]
    lambda_list = np.array(lambda_list)
    delta_l_list = np.array(delta_l_list)
    if rsd:
        eta_list = np.array(eta_list)
        if dla:
            velo_list = np.array(velo_list)
    redshift_list = np.array(redshift_list)

    # # Bug check
    # id_idx = util.find_A_in_B(np.array(QSOid_list), qsos['THING_ID'])
    # c1 = np.array_equal(np.array(ra_list), qsos['RA'][id_idx])
    # c2 = np.array_equal(np.array(dec_list), qsos['DEC'][id_idx])
    # c3 = np.array_equal(np.array(zQSO_norsd_list), qsos['Z_QSO_NO_RSD'][id_idx])
    # c4 = np.array_equal(np.array(zQSO_rsd_list), qsos['Z_QSO_RSD'][id_idx])
    # if not c1 & c2 & c3 & c4:
    #     print("WARNING: IDS DOES NOT CORRESPOND:")
    #     print("c1: {}; c2: {}; c3: {}; c4: {}".format(c1,c2,c3,c4))
    #     print("Saving in {}".format(args.outDir))
    #     np.save(args.outDir+"/ra.npy", np.array(ra_list))
    #     np.save(args.outDir+"/dec.npy", np.array(dec_list))
    #     np.save(args.outDir+"/z.npy", np.array(zQSO_norsd_list))
    #     np.save(args.outDir+"/zrsd.npy", np.array(zQSO_rsd_list))
    #     np.save(args.outDir+"/wav.npy", lambda_list)
    #     np.save(args.outDir+"/delta.npy", delta_l_list)
    #     np.save(args.outDir+"/redshift.npy", redshift_list)
    #     if rsd:
    #         np.save(args.outDir+"/eta.npy", eta_list)
    #         if dla:
    #             np.save(args.outDir+"/velo.npy", velo_list)

    for ID in np.unique(QSOhdu_list):
        outfits = fitsio.FITS(args.outDir+'/spectra-{}-{}.fits'.format(iSlice, ID), 'rw', clobber=True)
        msk = (QSOhdu_list == ID)
        table = [np.array(ra_list)[msk], np.array(dec_list)[msk],
                 np.array(zQSO_norsd_list)[msk], np.array(zQSO_rsd_list)[msk],
                 np.array(QSOhdu_list)[msk], np.array(QSOid_list)[msk],
                 np.array(plate_list)[msk], np.array(mjd_list)[msk],
                 np.array(fiber_list)[msk], np.array(pmf_list)[msk]]
        outfits.write(table, names=names, header=hlist, extname='METADATA')
        wav = list(lambda_list[msk])
        delta = list(delta_l_list[msk])
        if rsd:
            eta = list(eta_list[msk])
            if dla:
                velo = list(velo_list[msk])
        redshift = list(redshift_list[msk])
        for i in range(len(delta)):
            wav[i] = np.concatenate((wav[i], -1*np.ones(maxsize-len(wav[i]))))
            delta[i] = np.concatenate((delta[i], -2e6*np.ones(maxsize-len(delta[i]))))
            if rsd:
                eta[i] = np.concatenate((eta[i], -2e6*np.ones(maxsize-len(eta[i]))))
                if dla:
                    velo[i] = np.concatenate((velo[i], -2e6*np.ones(maxsize-len(velo[i]))))
            redshift[i] = np.concatenate((redshift[i], -1*np.ones(maxsize-len(redshift[i]))))
        outfits.write(np.float32(wav), extname='LAMBDA')
        outfits.write(np.float32(delta), extname='DELTA_L')
        if rsd:
            outfits.write(np.float32(eta), extname='ETA_PAR')
            if dla:
                outfits.write(np.float32(velo), extname='VELO_PAR')
        outfits.write(np.float32(redshift), extname='REDSHIFT')
        outfits.close()

    print("Done. {} s".format(time.time()-t1))
    
    if (save_wisdom) :
        wisd=pyfftw.export_wisdom()
        f = open(wisdomFile, 'w')
        json.dump(wisd,f)
        f.close()
    
    print (iqso, "QSO written")
    if (iqso != 0):
        print (1000*(t1-t0)/iqso, "ms per spectra")
        if (pixtot != 0) : print (1000*(t1-t0)/pixtot, "ms per pixel")
        print (pixtot," / ",iqso," = ", pixtot/iqso,"pixels / spectra")
    
    print("Slice {} done. Took {}s".format(iSlice, time.time()-t_init))

if __name__ == "__main__":
    main()
    sys.exit()

                                # anything bellow useful ?  <==
    
    #=====================================================================================
    
    #.............................................. QSO catalogue
    #   a la pyfits
    #col_ra = pyfits.Column(name='RA',  format='D', array=ralist, unit='deg')
    #col_de = pyfits.Column(name='DEC',  format='D', array=declist, unit='deg')
    #col_zz = pyfits.Column(name='Z',  format='D', array=zlist, unit='0 (redshift)')
    #tbhdu = pyfits.BinTableHDU.from_columns([col_ra, col_de, col_zz])
    #tbhdu.writeto('DesiMocks/delta/QSO.fits.gz', clobber=True)
    
    
    if (Flux):
        k,PF1D,errPF1D = MakeProfileHisto(kvec,PkFvec)
        plt.errorbar(k, PF1D, yerr=0, fmt='+')
        plt.plot(k,0*k)
    
        # McDonald z=2.4 
        k = np.array([0.00141, 0.00178, 0.00224, 0.00282, 0.00355, 0.00447, 0.00562, 0.00708, 0.00891, 0.01122, 0.01413, 0.01778])  
        P  = np.array([20.85, 22.94, 22.87, 21.63, 18.17, 15.40, 17.76, 13.00, 12.50, 10.42, 8.24, 6.18]) 
        errP  = np.array([1.85, 2.03, 1.35, 1.21, 0.87, 0.64, 0.67, 0.50, 0.41, 0.41, 0.36, 0.34]) 
        k *= 100
        P /= 100
        errP /= 100
        #plt.errorbar(k, P, yerr=errP, fmt='o',color='green')
        
        
        filename="data/PF1D-DR9.dat"
        data = np.loadtxt(filename, skiprows=1) 
        z = data[:,0]
        k = data[:,1]
        P = data[:,2]
        dP = data[:,3]
        
        beta=1.666
        cut = np.where(z==2.4)
        k=k[cut] * 100
        P=P[cut] / 100 #/(1+beta)**2
        dP=dP[cut] /100 #/(1+beta)**2
        plt.errorbar(k,P,yerr=dP,fmt='o',color='red')
        
        # from sinu1D.cpp (fit by C.Y>) but for k in km/s
        # 100 km/s ~ 1 Mpc/h
        kk = k/100
        b= -0.35;
        c=-0.075;
        P1Dexp = -3.224883 +b*np.log(kk)+c*np.log(kk)*np.log(kk);
        P1Dexp = PI * np.exp(P1Dexp) /kk
        P1Dexp /= 100
        #plt.plot(k,P1Dexp)
        
        P1Dexp = -3.039048 +b*np.log(kk)+c*np.log(kk)*np.log(kk);
        P1Dexp = PI * np.exp(P1Dexp) /kk
        P1Dexp /= 100
        #plt.plot(k,P1Dexp)
        
        #plt.show()
    
    k,P1D,errP1D = MakeProfileHisto(kvec,Pkvec)
    #plt.errorbar(k, P1D, yerr=errP1D, fmt='o')
    plt.errorbar(k, P1D, yerr=0, fmt='o')
    #print k
    
    kmax = 10
    dk=0.001
    kk=np.arange(kmax/dk)*dk
    Pcamb = powerspectrum.P_0("data/PlanckDR12.fits").P(kk)
    P1Dcamb = powerspectrum.P_1D(kk,Pcamb).P1D(k)
    plt.plot(k,P1Dcamb)
    plt.plot(k,0*k)
    
    kmax = k[-1]
    #print kmax
    #print kk
    cut = np.where(kk<=kmax)
    kk=kk[cut]
    Pcamb=Pcamb[cut]
    P1Dcutcamb = powerspectrum.P_1D(kk[cut],Pcamb[cut]).P1D(k)
    #plt.plot(k,P1Dcutcamb)
    
    
    DX = 1.6 
    W = (DX/np.sqrt(2*PI))* np.exp(- DX*DX*kk*kk/2)
    W = W/W[0]
    Pcamb *= W*W
    P1DWcutcamb = powerspectrum.P_1D(kk[cut],Pcamb[cut]).P1D(k)
    plt.plot(k,P1DWcutcamb)
    
    
    #dk=0.001
    #kk=np.arange(kmax/dk)*dk
    #Pcamb = powerspectrum.P_0("data/PkCamb.dat").P(kk)
    #P1Dcamb = powerspectrum.P_1D(kk,Pcamb).P1D(k)
    #plt.plot(k,P1Dcamb)
    
    
    #Pcmv = h*powerspectrum.P_0("data/Pk_cmv_z0.dat",colp=3).P(k*h)
    #P1Dcmv = powerspectrum.P_1D(k,Pcmv).P1D(k)
    #plt.plot(k,P1Dcmv)
    plt.show()
    
    
    '''
    XTENSION= 'BINTABLE'           / binary table extension
    BITPIX  =                    8 / 8-bit bytes
    NAXIS   =                    2 / 2-dimensional binary table
    NAXIS1  =                   32 / width of table in bytes
    NAXIS2  =                  133 / number of rows in table
    PCOUNT  =                    0 / size of special data area
    GCOUNT  =                    1 / one data group (required keyword)
    TFIELDS =                    4 / number of fields in each row
    TTYPE1  = 'LOGLAM  '           / label for field   1
    TFORM1  = 'D       '           / data format of field: 8-byte DOUBLE
    TTYPE2  = 'DELTA   '           / label for field   2
    TFORM2  = 'D       '           / data format of field: 8-byte DOUBLE
    TTYPE3  = 'WEIGHT  '           / label for field   3
    TFORM3  = 'D       '           / data format of field: 8-byte DOUBLE
    TTYPE4  = 'CONT    '           / label for field   4
    TFORM4  = 'D       '           / data format of field: 8-byte DOUBLE
    MJD     =                57162
    PLATE   =                 7431
    DEC     =    0.845550407798142
    THING_ID=            479651402
    FIBERID =                  770
    PMF     = '7431-57162-770'
    RA      =     3.59432233582233
    Z       =     2.30414414405823
    '''
    
    
    
        #        quasars = [data.quasar(x,y,z) for x,y,z in zip(tanx,tany,zzz)]
    
    #for qso in (qsos) :
    #  print qso.spectrum[0:5]
    
    #    print sumrho,sumweight,mean_rho
    #    print "weight min max:",weight.min(axis=0), weight.max(axis=0) 
      
    #  ic = 1
    #  print lcells[ic]
    #  print lcells[ic,0],lcells[ic,1],lcells[ic,2]
    #  print rho[lcells[ic,0],lcells[ic,1],lcells[ic,2]]
    #  print myrho[ic]
    #  sys.exit(0)
      
    
    
    
