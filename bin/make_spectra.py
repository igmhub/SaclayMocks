#!/usr/bin/env python
#!/usr/bin/env python
#  MakeSpectra.py
#	reads QSO list
# 	reads non parallel l.o.s. through boxes
#   add small scale power
#	applies Gunn Peterson to make Lya forest
from __future__ import division, print_function
from SaclayMocks import box
from SaclayMocks import constant
from SaclayMocks import powerspectrum
from SaclayMocks import util
import fitsio
from fitsio import FITS
import sys
import scipy as sp
import numpy as np
import cosmolopy.distance as dist
import os
import time
from numba import jit
import argparse
import scipy.stats as stats
import matplotlib.pyplot as plt
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
    def ComputeWeight(X,Y,Z, sig2) : 
        # returns local cells around (X,Y,Z) and Gaussian weights
        # also uses constants grid,cell,LX,LY,LZ,DX,DY,DZ
        # dmax=3 cell = array([[-3, -3, -3], [-3, -3, -2], [-3, -3, -1],
        #   ...,   [ 3,  3,  1], [ 3,  3,  2], [ 3,  3,  3]])   (343,3)
        # grid idem multiplied by DX,DY,DZ

        #  .......  cell that contains (X,Y,Z) and surrounding cells
        ix = int((X +LX/2)/DX)  # -LX/2 < X < -LX/2 + LX/Nslice so 0 < ix < NX/Nslice
        iy = int((Y +LY/2)/DY)  # -LY/2 < Y < LY/2 so 0 < iy < NY
        iz = int((Z +LZ/2 -R0)/DZ)
        ixyz = sp.array([ix,iy,iz])  	# (3,)
        lcells = cells + ixyz   # surrounding cells (343,3) 

        # weights
        cell_center = sp.array([(ix+0.5)*DX-LX/2,(iy+0.5)*DY-LY/2,
                (iz+0.5)*DZ-LZ/2+R0])  # (3,)
        xx = (grid+cell_center - sp.array([X,Y,Z]))**2  # (343,3)
        Delta_r2 = xx[:,0]+xx[:,1]+xx[:,2]    # 41-46 s (343,)
        #Delta_r2 = np.sum( ((XYZ-grid-cell_center)**2) ,axis=1 ) # 80-95 s
        #Delta_r2 = ((XYZ-grid-cell_center)**2).sum(axis=1) # 69 - 74 s
        # equivalent, but longer !!  (time for full MakeSpectra)
        weight = sp.exp(-Delta_r2 / sig2)
        return lcells, weight

    #*************************************************************
    #@jit   # @jit degrades from 22 to 27 ms
    def computeRho(myrho,weight) :
        return (weight*myrho).sum() / weight.sum()
        #sumweight = weight.sum()
        #sumrho = (weight*myrho).sum()
        #return sumrho / sumweight

    #*************************************************************
    # this seems marginally faster, significant ? 
    def computeRhob(rho,lcells,weight) :
        return (weight*rho[lcells[:,0],lcells[:,1],lcells[:,2]]).sum() / weight.sum()

    #*************************************************************
    #   @jit + python -m cProfile fails
    #   @jit degrades from 80 t0 94 for the full treatment of a QSO
    #@jit
    def ReadSpec(Xvec, XvecSlice, Yvec, Zvec, imin=0, imax=sys.maxint):
        # reads spectrum for (Xvec, Yvec, Zvec)
        # XvecSlice is in [-LX/2, -LX/2 + LX/NSlice]
        # cells is the list of indices used for G.S. around (0,0,0),
        # and grid its value in Mpc/h, both shapes are (343,3)
        # imin imax are the indices delimiting the lya forest
        # function also uses cosntants LX,LY,LZ,DX,DY,DZ

        spectrum = -1000000 * sp.ones_like(XvecSlice) # so that exp(-a(exp(b*g))) = 1
        if rsd:
            eta_par = sp.zeros_like(XvecSlice)
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
            lcells, weight = ComputeWeight(X,Y,Z, sig2)
            myrho = fullrho[lcells[:,0],lcells[:,1],lcells[:,2]]
            spectrum[icell] = computeRho(myrho,weight)
            #spectrum[icell] = computeRhob(fullrho,lcells,weight)
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
    parser.add_argument("-dgrowthfile", help="dD/dz file, default etc/dgrowth.fits", default=None)
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
    if args.dgrowthfile is None:
        filename = os.path.expandvars("$SACLAYMOCKS_BASE/etc/dgrowth.fits")
    if Om == fitsio.read_header(filename, ext=1)['OM']:
        Dgrowth = util.InterpFitsTable(filename, 'Z', 'dD/dz')
    else:
        raise ValueError("Omega_M_0 in SaclayMocks.constant ({}) != OM in {}".format(Om,
                            fitsio.read_header(filename, ext=1)['OM']))
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
    cells = sp.indices((xx,xx,xx)).reshape(3,-1) - dmax  # (3,343)
    grid = sp.array([DX*cells[0],DY*cells[1],DZ*cells[2]])
    cells = cells.T  # (343,3)
    grid = grid.T

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
                    delta_l, eta_par, velo_par = ReadSpec(Xvec, XvecSlice, Yvec, Zvec, imin=imin, imax=imax)
                except:
                    print("***WARNING ReadSpec:\n    ID {}***".format(QSOid))
                    continue
            else:
                try:
                    delta_l, eta_par = ReadSpec(Xvec, XvecSlice, Yvec, Zvec, imin=imin, imax=imax)
                except:
                    print("***WARNING ReadSpec:\n    ID {}***".format(QSOid))
                    continue
        else:
            delta_l = ReadSpec(Xvec, XvecSlice, Yvec, Zvec, imin=imin, imax=imax)

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
    print (iqso, "QSO written")
    if (iqso != 0):
        print (1000*(t1-t0)/iqso, "ms per spectra")
        if (pixtot != 0) : print (1000*(t1-t0)/pixtot, "ms per pixel")
        print (pixtot," / ",iqso," = ", pixtot/iqso,"pixels / spectra")

    print("Slice {} done. Took {}s".format(iSlice, time.time()-t_init))

if __name__ == "__main__":
    main()
