#!/usr/bin/env python
# python DrawQSO.py DesiMocks/box1024.fits
# draw a QSO in each pixel of the box with a proba proportional to exp(delta) and n(z)
# n(g,z) missing and no Gaussian smearing
import scipy as sp
import numpy as np
import healpy as hp
import sys
import matplotlib.pyplot as plt
import os
import cosmolopy.distance as dist
import fitsio
from fitsio import FITS,FITSHDR
from SaclayMocks import box
from SaclayMocks import constant
from SaclayMocks import util
import argparse
from time import time
from memory_profiler import profile
import gc
PI = np.pi

fullBox = False  # to fill the box completely without N(z) in order to check P(k)
drawPlot = False


#********************************************************************
# David Kirkby mail 29/6/16
# https://github.com/dkirkby/ArgonneLymanAlpha/blob/master/simulate/qso_grid.py
# also ... master/notebooks/SampleProperties.ipynb
def luminosity_function(data, z_max=6.0, area_sq_deg=10000.):
    """Transform a data array from Nathalie into a tuple gbin, zbin, nqso.
    """
    ng, nz = data.shape
    # g-band magnitude bin centers are in the first column.
    gbin = data[:, 0]
    nz = nz - 1
    # Check that g-band bins are equally spaced.
    assert np.allclose(np.diff(gbin),  gbin[1] - gbin[0])
    # redshift bins are equally spaced from 0 up to z_max.
    zbin = z_max * (0.5 + np.arange(nz)) / nz
    # The remaining columns give predicted numbers of QSO in a 10,000 sq.deg. sample.
    # Normalize to densisities per sq.deg.
    nqso = data[:, 1:].reshape((ng, nz)) / area_sq_deg
    return gbin, zbin, nqso

#********************************************************************
def bin_index(bin_centers, low_edge):
    """Find the index of the bin with the specified low edge, where bins is an array of equally-spaced bin centers.
    """
    delta = bin_centers[1] - bin_centers[0]
    min_value = bin_centers[0] - 0.5 * delta
    index = int(round((low_edge - min_value) / delta))
    if abs((low_edge - min_value) / delta - index) > 1e-5:
        raise ValueError('low_edge = {} is not aligned with specified bins.'.format(low_edge))
    return index

#********************************************************************
def sample(data, num_samples, g_min=19, g_max=23, z_min=0, z_max=6, seed=None):
    """Generate random samples of (g,z) within the specified cuts.
    """
    gbin, zbin, nqso = luminosity_function(data)
    z_min_cut = bin_index(zbin, z_min)
    z_max_cut = bin_index(zbin, z_max)
    g_min_cut = bin_index(gbin, g_min)
    g_max_cut = bin_index(gbin, g_max)
    nqso_cut = nqso[g_min_cut:g_max_cut, z_min_cut:z_max_cut]
    # Calculate the flattened CDF
    cdf = np.cumsum(nqso_cut.ravel())
    cdf /= np.float(cdf[-1])
    # Pick random flattened indices.
    generator = np.random.RandomState(seed)
    r = generator.rand(num_samples)
    indices = np.searchsorted(cdf, r)
    # Unravel to get g,z indices.
    g_indices, z_indices = np.unravel_index(indices, (len(gbin), len(zbin)))
    # Spread points out uniformly in each 2D bin.
    dg, dz = gbin[1] - gbin[0], zbin[1] - zbin[0]
    g = gbin[g_min_cut + g_indices] + dg * (generator.rand(num_samples) - 0.5)
    z = zbin[z_min_cut + z_indices] + dz * (generator.rand(num_samples) - 0.5)
    return g, z


#********************************************************************
#def zDistribution(gbin,zbin,nqso) :
#	np.linalg.matmut(nqso,)

#********************************************************************
# @profile
def main():
# if True:
#........................................  parameters and constants
    t_init = time()

    parser = argparse.ArgumentParser()
    parser.add_argument("-dmax", type=int, help="no quasar < dmax cells from edges, default 3", default=3)
    parser.add_argument("-indir", help="directory which contains boxes")
    parser.add_argument("-outpath", help="QSO fit path")
    parser.add_argument("-i", type=int, help="index of treated HDU")
    parser.add_argument("-Nslice", type=int, help="total number of slice")
    parser.add_argument("-chunk", type=int, help="index of treated chunk, default 0", default=0)
    parser.add_argument("-ra0", type=float, help="center box RA in degrees, default 0", default=0)
    parser.add_argument("-dec0", type=float, help="center box DEC in degrees, default 0", default=0)
    parser.add_argument("-dra", type=float, help="|ra-ra0|<dra in degrees, default -1 no cut", default=-1)
    parser.add_argument("-ddec", type=float, help="|dec-dec0|<ddec in degrees, default -1 no cut", default=-1)
    parser.add_argument("-random",type=str, help="If True, generate randoms. Default: False", default='False')
    parser.add_argument("-zmin", type=float, help="minimal redshift for drawing QSO, default is None", default=-1)
    parser.add_argument("-zmax", type=float, help="maximal redshift for drawing QSO, default is None", default=-1)
    parser.add_argument("-desi", type=str, help="select only objects in desi footprint, default is True", default="True")
    parser.add_argument("-seed", type=int, help="specify a seed", default=None)
    parser.add_argument("-rsd", help="If True, rsd are added, default True", default='True')
    parser.add_argument("-dgrowthfile", help="dD/dz file, default etc/dgrowth.fits", default=None)
    args = parser.parse_args()
    zmin = args.zmin
    zmax = args.zmax
    dmax = args.dmax    # margin of dmax cells
    i_slice = args.i
    chunk = args.chunk
    ra0 = args.ra0
    dec0 = args.dec0
    dra = args.dra
    ddec= args.ddec
    random_cond = util.str2bool(args.random)
    desi = util.str2bool(args.desi)
    rsd = util.str2bool(args.rsd)

    h = constant.h
    H0 = constant.H0
    Om = constant.omega_M_0
    OL = constant.omega_lambda_0
    Ok = constant.omega_k_0
    z0 = constant.z0
    if args.dgrowthfile is None:
        filename = os.path.expandvars("$SACLAYMOCKS_BASE/etc/dgrowth.fits")
    if Om == fitsio.read_header(filename, ext=1)['OM']:
        Dgrowth = util.InterpFitsTable(filename, 'Z', 'dD/dz')
    else:
        raise ValueError("Omega_M_0 in SaclayMocks.constant ({}) != OM in {}".format(Om,
                            fitsio.read_header(filename, ext=1)['OM']))
    dgrowth0 = Dgrowth.interp(0)

    print("\n\nBegining of DrawQSO.")
    if random_cond:
        print("***Random catalogs***")
        print(random_cond)

    global seed
    seed = args.seed
    if seed is None:
        seed = np.random.randint(2**31 -1, size=1)[0]
        np.random.seed(seed)
        print("Seed has not been specified. Seed is set to {}".format(seed))
    else:
        np.random.seed(seed)
        print("Specified seed is {}".format(seed))

    #........................................   read box.fits header
    boxfits = fitsio.FITS(args.indir+"/boxln-{}.fits".format(i_slice))
    head = boxfits[0].read_header()
    Nslice = args.Nslice
    DX = head["DX"]
    DY = head["DY"]
    DZ = head["DZ"]
    NZ = head["NAXIS1"]
    NY = head["NAXIS2"]
    NX = head["NAXIS3"]
    NX_fullbox = head["NX"]

    print("Treating: slice = {} ; Nslice = {}\n".format(i_slice, Nslice))
    # out_file = args.out
    if random_cond:
        out_file = args.outpath+'/randoms-{}-{}.fits'.format(i_slice, Nslice)
    else:
        out_file = args.outpath+'/QSO-{}-{}.fits'.format(i_slice, Nslice)

    if (not random_cond):
        t0=time()
        rho = boxfits[0].read()
        boxfits.close()
        t1=time()
        print "read box in ",t1-t0,"s ",rho.shape
        print "sigma(rho)=",rho.std()
        exprho = np.exp(rho)
        del rho  # unassign variable
        gc.collect()  # force memory clearing
        if rsd:
            print("Reading velocity boxes...")
            vx = fitsio.read(args.indir+"/vx-{}.fits".format(i_slice), ext=0)
            vy = fitsio.read(args.indir+"/vy-{}.fits".format(i_slice), ext=0)
            vz = fitsio.read(args.indir+"/vz-{}.fits".format(i_slice), ext=0)

    LX = NX * DX
    LX_fullbox = NX_fullbox * DX
    LY = NY * DY
    LZ = NZ * DZ

    #........................................     some geometry and cosmology
    cosmo_fid = {'omega_M_0':Om, 'omega_lambda_0':OL, 'omega_k_0':Ok, 'h':h}
    R_of_z, z_of_R = dist.quick_distance_function(dist.comoving_distance, return_inverse=True, **cosmo_fid)
    R0 = h * R_of_z(z0)
    x_axis = (sp.arange(NX)+0.5)*DX + (2*i_slice-Nslice) * LX_fullbox/(2*Nslice)
    y_axis = (sp.arange(NY)+0.5)*DY - LY/2
    z_axis = (sp.arange(NZ)+0.5)*DZ + R0 - LZ/2     # Z at cell center
    z_edges = sp.arange(NZ+1)*DZ + R0 - LZ/2        # Z at edges
    dz = z_of_R(z_edges[1:]/h) - z_of_R(z_edges[0:-1]/h)   #  z_of_R[m+1] - z_of_R[m]

    # margin of dmax cells
    print LX,LY,LZ,R0,dmax*DX # prov
    Rmin,Rmax,tanx_max,tany_max = box.box_limit(LX_fullbox,LY,LZ,R0,dmax*DX)
    # tanx_max and tany_max should be calculated again with new Rmin/Rmax ?
    if zmin > 0:
        z_min = max(z_of_R(Rmin/h), zmin)
        Rmin = R_of_z(z_min)*h
    else:
        z_min = z_of_R(Rmin/h)
    if zmax > 0:
        z_max = min(z_of_R(Rmax/h), zmax)
        Rmax = R_of_z(z_max)*h
    else:
        z_max = z_of_R(Rmax/h)
    print "zmin, z0, zmax:", z_min, z0, z_max
    print "Rmin, R0, Rmax:", Rmin, R0, Rmax
    print "theta_max:",tanx_max,tany_max
    tetax_max = np.arctan(tanx_max)
    tetay_max = np.arctan(tany_max)
    if ((dra>0) & (ddec>0)) :
        Ldec = np.sin(np.radians(dec0+ddec)) - np.sin(np.radians(dec0-ddec))
        Lra = 2*np.radians(dra)	
        surface = Ldec * Lra	# rad^2
        print ddec,dra,Ldec,Lra,surface
    else :
        surface = 4 * tetax_max * tetay_max     # rad^2
    surfaceDeg = surface * (180/PI)**2
    nQSOexp = constant.n_qso_exp * surfaceDeg
    nQSOexp /= Nslice
    volFrac= (surface/3)*(Rmax**3-Rmin**3)/LX_fullbox/LY/LZ
    print "Fraction of box volume used",volFrac
    print "surface: ", surfaceDeg, "deg^2  -> ",int(nQSOexp),"QSOs expected"
    xx = np.sqrt(Rmax*Rmax+LY*LY+LZ*LZ)
    print "far corner of the box",xx,"Mpc/h -> z=",z_of_R(xx/h)

    #...............................  read dN/dz assuming constant Delta z
    filename = os.path.expandvars("$SACLAYMOCKS_BASE/etc/nz_qso_desi.dat")
    d = np.loadtxt(filename)
    delta_z = d[1,0]-d[0,0]
    zlow = d[:,0]
    zhigh = d[:,1]
    dndz = d[:,2] #number of QSO per sq deg per dz=0.1
    if (z_of_R((R0+LZ/2)/h) > zhigh[-1]) :
        print "z_max = ",z_of_R((R0+LZ/2)/h)," larger than dNdz upper range ",zhigh[-1]
        exit(0)

    # smooth dndz
    f_dndz = sp.interpolate.interp1d((zlow+zhigh)/2, dndz)
    dz_interp = np.linspace((zlow[0]+zhigh[0])/2, (zlow[-1]+zhigh[-1])/2, 2*NZ)
    dndz_interp = f_dndz(dz_interp)

    # dn/dz/deg^2 => n/cell(z)
    dn_cell = dndz_interp * DZ / dist.hubble_distance_z(dz_interp,**cosmo_fid) 	
    # dZ = (c/H(z) dz = d_H dz   =>  dN/dZ = dN/dz Dz/dZ = dN/dz /d_H
    dn_cell = dn_cell * DX * DY / R_of_z(dz_interp)**2  # deg^-2 -> per cell

    density_max = dn_cell.max()
    density_mean = dn_cell[(dz_interp>z_min)*(dz_interp<z_max)].mean()
    if z_max > 2.1:  # Count only QSO for z > 2.1 to have the right N/deg^2
        N_zmin_zmax = dn_cell[(dz_interp>z_min)*(dz_interp<z_max)].sum()
        N_21_zmax = dn_cell[(dz_interp>2.1)*(dz_interp<z_max)].sum()
    print "dN per cell max:",density_max,", mean:", density_mean,density_max/density_mean
    dn_cell = np.append(dn_cell, np.zeros(10*NZ))  # artificially increasing dNdz range

    if (drawPlot) :
        hrho , hh = np.histogram(rho,100)
        plt.plot(hrho)
        plt.hist(exprho.reshape(np.size(exprho)),1000)
        plt.xscale('log')
        plt.yscale('log')
        plt.show()

    if (not random_cond):
        rho_max = exprho.max()
        rho_sum = exprho.sum()
        kk = nQSOexp / rho_sum
        norm = nQSOexp / rho_sum    # corresponds to cond1
        norm *= density_max / density_mean  # correction for cond2
        norm /= volFrac             # correction for cond3
        if z_max > 2.1:
            norm /= N_21_zmax / N_zmin_zmax  # correct for z > 2.1 only
        # gives ~4% more QSO than expected
        # this is due to the fact that the 3 cond are not independent

        print "exp(rho) max and sum = ",rho_max,rho_sum
        print "k exp(rho_max) =",kk*rho_max
        if (kk*rho_max>1):
            print "k exp(rho) > 1 in ", np.size(np.where(kk*exprho>1)[0]),"cells"
            print "sum min(k exp(rho) , 1) =", np.minimum(kk*exprho,1).sum()

    #.........................................................    loop on cells
    nQSO = 0
    if (not random_cond):
        nnQSO=0
    qsofits = FITS(out_file, 'rw', clobber=True)  # output file
    if not random_cond:
        names = ["Z_QSO_NO_RSD", "Z_QSO_RSD", "RA", "DEC", "HDU", "THING_ID", "PLATE", "MJD", "FIBERID", "PMF", "XX", "YY", "ZZ"] # , "XGRID", "YGRID", "ZGRID"]
    else:
        names = ["Z", "RA", "DEC", "HDU", "THING_ID", "PLATE", "MJD", "FIBERID", "PMF", "XX", "YY", "ZZ"] # , "XGRID", "YGRID", "ZGRID"]

    for mz in range(NZ):
        XX = x_axis
        YY = y_axis
        XY2 = (XX*XX).reshape(-1,1) + YY*YY	# broadcasting -> (NX,NY)
        ZZ = z_axis[mz]
        RR = np.sqrt(ZZ*ZZ + XY2)
        redshift = z_of_R(RR/h)
        delta_z = dz_interp[1] - dz_interp[0]
        iz = ((redshift - dz_interp[0]) / delta_z).round().astype(int)
        density = dn_cell[iz]

        # ==> should correct for the fact that   rnd1 < exp(rho)   not always true
        #  use reproducible random <==
        rnd1 = sp.random.ranf(size=(NX,NY))		#  float64
        if (not random_cond):
            cond1 = rnd1 < norm * exprho[:, :, mz]  # (NX,NY)
            # should be a Poisson of norm * np.exp(rho)  <==
            # sometime get 2 QSO in a cell
            nnQSO += np.size(np.where(cond1)[0])
        else:
            cond1 = rnd1 > (1.-0.006)

        cond2 = density_max * sp.random.ranf(size=(NX,NY)) < density	# (NX,NY)	

        # Draw random xyz in the cell
        XX = XX + np.random.uniform(-DX/2, DX/2, size=len(XX))
        YY = YY + np.random.uniform(-DY/2, DY/2, size=len(YY))

        XXX = XX.reshape(-1, 1) * np.ones(NY)  # (NX, NY)
        YYY = np.ones(NX).reshape(-1, 1) * YY
        ZZZ = ZZ + np.random.uniform(-DZ/2, DZ/2, size=[len(XX), len(YY)])
        ra, dec, RR = box.ComputeRaDecR2(XXX, YYY, ZZZ, np.radians(ra0), np.radians(dec0))
        ra = np.degrees(ra)
        dec = np.degrees(dec)
        redshift = z_of_R(RR/h)
        if redshift.min() > z_max +1.: continue
        if not random_cond and rsd:
            vpar = (XX.reshape(-1,1)*vx[:, :, mz]
                  + YY*vy[:, :, mz]
                  + ZZ*vz[:, :, mz]) / RR
            msk = np.where(redshift < z_max + 1.)  # don't go above z=5
            RR_RSD = RR.copy()
            RR_RSD[msk] += vpar[msk] * (1+redshift[msk]) * Dgrowth.interp(redshift[msk]) / (dgrowth0 * H0)
            redshift_RSD = z_of_R(RR_RSD/h)

        if rsd and not random_cond:
            cond3 = (util.diffmod(ra,ra0,360)<dra) * (util.diffmod(dec,dec0,180)<ddec) * (redshift_RSD>z_min) * (redshift_RSD<z_max)    # check norm still ok <==
        else:
            cond3 = (util.diffmod(ra,ra0,360)<dra) * (util.diffmod(dec,dec0,180)<ddec) * (redshift>z_min) * (redshift<z_max)    # check norm still ok <==

        if (fullBox):
            iqso = sp.where(cond1)
        else:
            iqso = sp.where(cond1 * cond2 * cond3)  # tuple (2,N_QSO)  ([iqso],[jqso])
            #            iqso = sp.where( cond1 * cond3 ) # prov
        if drawPlot:
            xx = x_axis[iqso[0]]
            yy = y_axis[iqso[1]]
            plt.plot(xx, yy, ls='none', marker='s')
            # plt.show()

        XX = XX[iqso[0]]
        YY = YY[iqso[1]]
        ZZ = ZZZ[iqso]
        RR = RR[iqso]
        zzz = z_of_R(RR/h)
        if RR.size == 0: continue
        if rsd and not random_cond:
            RR_RSD = RR_RSD[iqso]
            zzz_RSD = z_of_R(RR_RSD/h)
        else:
            zzz_RSD = zzz

        ra = ra[iqso]
        dec = dec[iqso]

        # # Draw random xyz in the cell
        # X = x_axis[iqso[0]] + np.random.uniform(-DX/2, DX/2, len(iqso[0]))
        # Y = y_axis[iqso[1]] + np.random.uniform(-DY/2, DY/2, len(iqso[0]))
        # Z = ZZ + np.random.uniform(-DZ/2, DZ/2, len(iqso[0]))
        # zzz = z_of_R(np.sqrt(X*X+Y*Y+Z*Z)/h)

        # ra, dec, _ = box.ComputeRaDecR2(XX, YY, ZZ, np.radians(ra0), np.radians(dec0))  # (NQS0_)
        # ra = np.degrees(ra)
        # dec = np.degrees(dec)

        # Select desi footprint
        if desi:
            msk = util.desi_footprint(ra, dec)
            ra = ra[msk]
            dec = dec[msk]
            zzz = zzz[msk]
            XX = XX[msk]
            YY = YY[msk]
            ZZ = ZZ[msk]
            zzz_RSD = zzz_RSD[msk]

        un = np.ones(len(zzz)).astype(int)
        # XGRID = np.array(iqso[0]) + i_slice*NX
        # YGRID = np.array(iqso[1])
        # ZGRID = un * mz
        #if (len(XGRID)>0):
        #    print ( np.mean(np.log(exprho[XGRID,YGRID,ZGRID])) )
        THING_ID = (chunk*1e9 + i_slice*1e6 + np.arange(nQSO, nQSO+len(zzz)) + 1).astype(int)  # start at 1
        HDU = un * i_slice  # QSOhdu
        plate = THING_ID
        mjd = np.random.randint(51608, high=57521, size=len(zzz))
        fiberid = np.random.randint(1,high=1001, size=len(zzz))
        # pmf = np.array([plate[i].astype(str)+'-'+mjd[i].astype(str)+'-'+fiberid[i].astype(str) for i in range(len(zzz))])
        pmf = np.empty((len(zzz)), dtype='|S21')
        for i in range(len(zzz)):
            l_tmp = plate[i].astype(str)+'-'+mjd[i].astype(str)+'-'+fiberid[i].astype(str)
            pmf[i] = l_tmp
        nQSO += len(zzz)
        # print(mz)
        # print(len(zzz), len(tanx), len(X))
        #print (len(zzz),len(ra),len(dec))
        if not random_cond:
            array_list = [np.float32(zzz), np.float32(zzz_RSD), np.float32(ra), np.float32(dec), np.int32(HDU), THING_ID, plate, np.int32(mjd), np.int32(fiberid), pmf, np.float32(XX), np.float32(YY), np.float32(ZZ)]  # , XGRID, YGRID, ZGRID]
        else:
            array_list = [np.float32(zzz), np.float32(ra), np.float32(dec), np.int32(HDU), THING_ID, plate, np.int32(mjd), np.int32(fiberid), pmf, np.float32(XX), np.float32(YY), np.float32(ZZ)]  # , XGRID, YGRID, ZGRID]

        # print("MZ : {}".format(mz))
        # print(array_list)

        if (len(qsofits) == 1):
            qsofits.write(array_list, names=names)
            #exit(0)
        else:
            qsofits[-1].append(array_list, names=names)
        #   this is not making a new HDU, this is appending to last hdu, i.e. hdu=1

    if len(qsofits) == 1:
        print("No QSO drawn, creating a null table")
        if not random_cond:
            null_table = [np.array([], dtype=np.float32), np.array([], dtype=np.float32), np.array([], dtype=np.float32), np.array([], dtype=np.float32), np.array([], dtype=np.int32), np.array([], dtype=np.int64), np.array([], dtype=np.int64), np.array([], dtype=np.int32), np.array([], dtype=np.int32), np.array([], dtype='|S21'), np.array([], dtype=np.float32), np.array([], dtype=np.float32), np.array([], dtype=np.float32)]
        else:
            null_table = [np.array([], dtype=np.float32), np.array([], dtype=np.float32), np.array([], dtype=np.float32), np.array([], dtype=np.int32), np.array([], dtype=np.int64), np.array([], dtype=np.int64), np.array([], dtype=np.int32), np.array([], dtype=np.int32), np.array([], dtype='|S21'), np.array([], dtype=np.float32), np.array([], dtype=np.float32), np.array([], dtype=np.float32)]

        qsofits.write(null_table, names=names)

    qsofits[1].write_key("seed", np.int32(seed), comment="seed used to generate randoms")
    qsofits[1].write_key("ra0", ra0, comment="right ascension of the box center")
    qsofits[1].write_key("dec0", dec0, comment="declination of the box center")
    qsofits[1].write_key("RA", None, comment="right ascension in degrees")
    qsofits[1].write_key("DEC", None, comment="declination in degrees")
    if not random_cond:
        qsofits[1].write_key("Z_QSO_NO_RSD", None, comment="redshift without rsd")
        qsofits[1].write_key("Z_QSO_RSD", None, comment="redshift with rsd")
    else:
        qsofits[1].write_key("Z", None, comment="redshift")

    qsofits[1].write_key("HDU", None, comment="Slice index of the QSO")
    qsofits[1].write_key("XX", None, comment="position on X axis in Mpc/h")
    qsofits[1].write_key("YY", None, comment="position on Y axis in Mpc/h")
    qsofits[1].write_key("ZZ", None, comment="position on Z axis in Mpc/h")
    # qsofits[1].write_key("XGRID", None, comment="index on X axis")
    # qsofits[1].write_key("YGRID", None, comment="index on Y axis")
    # qsofits[1].write_key("ZGRID", None, comment="index on Z axis")

    qsofits.close()
    print nQSO, "QSOs drawn"
    if (not random_cond):
        print nnQSO, "QSOs in the full box" # prov
    print out_file, "file written"
    print("Took {}s".format(time()-t_init))

    if drawPlot : plt.show()

if __name__ == "__main__":
    main()
