#!/usr/bin/env python
from __future__ import print_function, division
import numpy as np
from scipy.stats import norm
from scipy.interpolate import interp1d, interp2d
import astropy.table
import os, sys
import gc
import fitsio
import time
import glob
import argparse
from SaclayMocks import constant, util
# import cosmolopy.distance as dist
# from memory_profiler import profile


try:
    from pyigm.fN.fnmodel import FNModel
    fN_default = FNModel.default_model()
    fN_default.zmnx = (0.,4)
    fN_cosmo = fN_default.cosmo
    use_pyigm = True
    print("Using pyigm")
except:
    use_pyigm = False


def dz_of_z_func(cell_size=2.19, zmin=1.3, zmax=4., nbin=500):
    h = constant.h
    Om = constant.omega_M_0
    Ok = constant.omega_k_0
    # cosmo_fid = {'omega_M_0':Om, 'omega_lambda_0':OL, 'omega_k_0':Ok, 'h':h}
    # R_of_z, z_of_R = dist.quick_distance_function(dist.comoving_distance, return_inverse=True, **cosmo_fid)
    cosmo_fid = util.cosmo(Om, Ok=Ok, H0=100*h)
    R_of_z = cosmo_fid.r_comoving
    z_of_R = cosmo_fid.r_2_z
    rmin = R_of_z(zmin)*h
    rmax = R_of_z(zmax)*h
    r_vec = np.linspace(rmin, rmax, nbin)
    z_vec = z_of_R(r_vec/h)
    dz = z_of_R((r_vec + cell_size)/h) - z_of_R(r_vec/h)
    return interp1d(z_vec, dz)


def doloop(dlas_in_cell,velocity,zedges, dla_z, dla_skw_id, dla_rsd_dz, dla_count):
    """ Auxiliary function to perform the loop to populate the DLA cells"""
    for skw_id,dla in enumerate(dlas_in_cell):
        #Find cells that will be allocated at least one DLA
        dla_cells = np.where(dla>0)[0]
        #For each dla, assign it a redshift, a velocity and a column density.
        for cell in dla_cells:
            dla_z[dla_count:dla_count+dla[cell]] = zedges[cell]+(zedges[cell+1]-zedges[cell])*np.random.random(size=dla[cell])
            dla_skw_id[dla_count:dla_count+dla[cell]] = skw_id
            dla_rsd_dz[dla_count:dla_count+dla[cell]] = velocity[skw_id,cell]/constant.c
            dla_count = dla_count+dla[cell]
    return dla_z, dla_skw_id, dla_rsd_dz, dla_count

def nu_of_bD(b):
    """ Compute the Gaussian field threshold for a given bias"""
    nu = np.linspace(-10,100,500) # Generous range to interpolate
    # use something non linear for nu ? To get more points close to 0
    p_nu = norm.pdf(nu)
    galaxy_mean = norm.cdf(-nu)  # probability to be above threshold
                                # this gives zero for nu > 37.5
    b_nu = np.zeros(nu.shape)
    b_nu[galaxy_mean!=0] = p_nu[galaxy_mean!=0]/galaxy_mean[galaxy_mean!=0]
    # it means that to get a bias of 2, you need a value of the field
    # to be 2 times more probable than the probability to be above this value
    # b_nu[galaxy_mean==0] = nu[galaxy_mean==0]
        # approximation for nu > 37.5, better than 0.027, i.e 0.07%
    y = interp1d(b_nu,nu)
    return y(b)

# # Not used
# def get_bias_z(fname,dla_bias):
#     """ Given a path, read the z array there and return a bias inversely
#     proportional to the growth"""
#     colore_cosmo = fits.open(fname)[4].data
#     z = colore_cosmo['Z']
#     D = colore_cosmo['D']
#     y = interp1d(z,D)
#     bias = dla_bias/D*y(2.25)
#     return z, bias, D

# # Not used
# def get_sigma_g(fname, mode='SG'):
#     if mode=='SG':
#         # Biased as well
#         return fits.open(fname)[4].header['SIGMA_G']
#     if mode=='SKW':
#         # Approximation 2: Take the skewers (biased when QSOs are present)
#         skewers = fits.open(fname)[2].data
#         return np.std(skewers,axis=0)

def flag_DLA(zq,z_cells,deltas,nu_arr,sigma_g,zlow,dz_of_z, rand=False):
    """ Flag the pixels in a skewer where DLAs are possible"""
    # find cells with density above threshold
    if not rand:
        flag = deltas > nu_arr*sigma_g  # (nspec, npix)
    else:
        flag = np.bool_(np.ones_like(deltas))  # (nspec, npix)
    # mask cells with z > z_qso, where DLAs would not be observed
    Nq=len(zq)
    for i in range(Nq):
        dz = dz_of_z(zq[i])  # avoid drawing DLA in same cell as QSO
        low_z = (z_cells < zq[i]-dz) & (z_cells > zlow)
        flag[i,:] *= low_z
    return flag

#number per unit redshift from minimum lg(N) in file (17.2) to argument
# Reading file from https://arxiv.org/pdf/astro-ph/0407378.pdf

def dnHD_dz_cumlgN(z,logN):
    tab = astropy.table.Table.read(os.path.abspath('LyaCoLoRe/example_data/zheng_cumulative.overz'),format='ascii')
    y = interp2d(tab['col1'],tab['col2'],tab['col3'],fill_value=None)
    return y(z,logN)


def dNdz(z, Nmin=20.0, Nmax=22.5):
    """ Get the column density distribution as a function of z,
    for a given range in N"""
    if use_pyigm:
        print("Use pyigm to compute dNdz.")
        # get incidence rate per path length dX (in comoving coordinates)
        dNdX = fN_default.calculate_lox(z,Nmin,Nmax)
        # convert dX to dz
        dXdz = fN_cosmo.abs_distance_integrand(z)
        dndz = dNdX * dXdz
        return dndz
    else:
        return dnHD_dz_cumlgN(z,Nmax)-dnHD_dz_cumlgN(z,Nmin)


# def dNdz(z, Nmin=20.0, Nmax=22.5, nsamp=100):
#     """ Get the column density distribution as a function of z,
#     for a given range in N"""
#     # get incidence rate per path length dX (in comoving coordinates)
#     nn = np.linspace(Nmin,Nmax,nsamp)
#     aux = fN_default.evaluate(nn, z)
#     dNdz = np.sum(np.exp(aux)*(nn[1]-nn[0]))
#     return dNdz


def get_N(z, Nmin=20.0, Nmax=22.5, nsamp=100):
    """ Get random column densities for a given z
    """

    # number of DLAs we want to generate
    Nz = len(z)
    nn = np.linspace(Nmin,Nmax,nsamp)
    probs = np.zeros([Nz,nsamp])
    if use_pyigm:
        auxfN = fN_default.evaluate(nn,z)
        #auxfN = (np.cumsum(10**auxfN, axis=0)/np.sum(10**auxfN, axis=0)).T
        probs = (np.exp(auxfN)/np.sum(np.exp(auxfN), axis=0)).T
        #plt.plot(nn,auxfN.T)
    else:
        probs_low = dnHD_dz_cumlgN(z,nn[:-1]).T
        probs_high = dnHD_dz_cumlgN(z,nn[1:]).T
        probs[:,1:] = probs_high-probs_low
    NHI = np.zeros(Nz)
    for i in range(Nz):
        #if use_pyigm:
        #    nfunc = interp1d(auxfN[i],nn,fill_value='extrapolate')
        #    NHI[i] = nfunc(np.random.uniform())
        #else:
        #    NHI[i] = np.random.choice(nn,size=1,p=probs[i]/np.sum(probs[i]))+(nn[1]-nn[0])*np.random.random(size=1)
        NHI[i] = np.random.choice(nn,size=1,p=probs[i]/np.sum(probs[i]))+(nn[1]-nn[0])*np.random.random(size=1)
    return NHI


def get_NHI(z, NHI_min=17.2, NHI_max=22.5, NHI_nsamp=100):
    """ Get random column densities for a given z
    This function is taken from LyaColoRe
    https://github.com/igmhub/LyaCoLoRe/blob/master/py/DLA.py
    """
    times = []
    t = time.time()

    # number of DLAs we want to generate
    Nz = len(z)

    # Set up the grid in NHI, and define its edges/widths.
    # First in log space.
    log_NHI_edges = np.linspace(NHI_min,NHI_max,NHI_nsamp+1)
    log_NHI = (log_NHI_edges[1:] + log_NHI_edges[:-1])/2.
    log_NHI_widths = log_NHI_edges[1:] - log_NHI_edges[:-1]
    # Then in linear space.
    NHI_edges = 10**log_NHI_edges
    NHI_widths = NHI_edges[1:] - NHI_edges[:-1]

    times += [time.time()-t]
    t = time.time()

    probs = np.zeros([Nz,NHI_nsamp])

    if use_pyigm:

        #Evaluate f at the points of the NHI grid and each redshift.
        f = 10**fN_default.evaluate(log_NHI,z)

        times += [time.time()-t]
        t = time.time()
        
        #Calculate the probaility of each NHI bin.
        aux = f*np.outer(NHI_widths,np.ones(z.shape))
        
        times += [time.time()-t]
        t = time.time()
    
        probs = (aux/np.sum(aux,axis=0)).T

        times += [time.time()-t]
        t = time.time()

    else:

        # TODO: test this
        probs_low = dnHD_dz_cumlgN(z,nn[:-1]).T
        probs_high = dnHD_dz_cumlgN(z,nn[1:]).T
        probs[:,1:] = probs_high-probs_low

    times += [time.time()-t]
    t = time.time()

    #Calculate the cumulative distribution
    """
    cumulative = np.zeros(probs.shape)
    for i in range(NHI_nsamp):
        cumulative[:,i] = np.sum(probs[:,:i],axis=1)
    """
    cumulative = np.cumsum(probs,axis=1)

    #Add the top and bottom edges on to improve interpolation.
    """
    log_NHI_interp = np.concatenate([[log_NHI_edges[0]],log_NHI,[log_NHI_edges[1]]])
    end_0 = np.zeros((z.shape[0],1))
    end_1 = np.ones((z.shape[0],1))
    cumulative_interp = np.concatenate([end_0,cumulative,end_1],axis=1)
    """

    log_NHI_interp = log_NHI_edges
    end_0 = np.zeros((z.shape[0],1))
    cumulative_interp = np.concatenate([end_0,cumulative],axis=1)

    times += [time.time()-t]
    t = time.time()

    #Assign NHI values by choosing a random number in [0,1] and interpolating
    #the cumulative distribution to get a value of NHI.
    log_NHI_values = np.zeros(Nz)
    for i in range(Nz):
        p = np.random.uniform()
        log_NHI_values[i] = np.interp(p,cumulative_interp[i,:],log_NHI_interp)

    times += [time.time()-t]
    t = time.time()

    times = np.array(times)
    times /= np.sum(times)
    #print(times.round(4))

    return log_NHI_values


# @profile
def add_DLA_table_to_object_Saclay(hdulist,dNdz_arr,dz_of_z,dla_bias=2.0,extrapolate_z_down=None,Nmin=20.0,Nmax=22.5,zlow=1.8, rand=False, nhi_low_cut=None, nhi_high_cut=None):
    qso = hdulist['METADATA'].read() # Read the QSO table
    lam = hdulist['LAMBDA'].read() # Read the vector with the wavelenghts corresponding to each cell
    deltas = hdulist['DELTA'].read()  # (nspec, npix)
    velocity = hdulist['VELO_PAR'].read()  # (nspec, npix)
    #Linear growth rate of each cell in the skewer
    D_cell = hdulist['GROWTHF'].read()  # (npix)
    hdulist.close()
    #Quasar redshift for each skewer
    zq = qso['Z']  # (nspec)
    #Redshift of each cell in the skewer
    z_cell = lam / constant.lya - 1  # (npix)
    # # Not use
    # # Read cosmo
    # cosmo_hdu = fitsio.read_header(fname_cosmo, ext=1) # Reading the cosmological parameters used for the simulation
    # Oc = cosmo_hdu['OM']-cosmo_hdu['OB'] # Omega_c
    # Ob = cosmo_hdu['OB'] # Omega_b
    # h = cosmo_hdu['H'] # h
    # Ok = cosmo_hdu['OK'] # Omega_k
    #Setup bias as a function of redshift
    # y = interp1d(z_cell,D_cell)
    # bias = dla_bias/(D_cell)*y(2.25)  # (npix)
    # sigma_g = fitsio.FITS(fname_sigma)[0].read_header()['SIGMA']
    sigma_g = constant.sigma_l
    # Gaussian field threshold:
    nu_arr = nu_of_bD(dla_bias*sigma_g*D_cell)  # (npix)
    #Figure out cells that could host a DLA, based on Gaussian fluctuation
    flagged_cells = flag_DLA(zq,z_cell,deltas,nu_arr,sigma_g,zlow, dz_of_z, rand)
    flagged_cells[deltas==-1e6]=False  # don't draw DLA outside forest
    #Edges of the z bins
    if extrapolate_z_down and extrapolate_z_down<z_cell[0]:
        zedges = np.concatenate([[extrapolate_z_down],(z_cell[1:]+z_cell[:-1])*0.5,[z_cell[-1]+(-z_cell[-2]+z_cell[-1])*0.5]]).ravel()
    else:
        zedges = np.concatenate([[z_cell[0]],(z_cell[1:]+z_cell[:-1])*0.5,[z_cell[-1]+(-z_cell[-2]+z_cell[-1])*0.5]]).ravel()
    z_width = zedges[1:]-zedges[:-1]

    #Get the average number of DLAs per cell, from the column density dist.
    mean_N_per_cell = z_width*dNdz_arr  # (npix)
    #For a given z, probability of having the density higher than the threshold
    p_nu_z = 1.0-norm.cdf(nu_arr)

    #Define mean of the Poisson distribution (per cell)
    mu = mean_N_per_cell * np.ones_like(flagged_cells)
    # correct for the threshold dependency
    if not rand:
        mu /= p_nu_z
    # mu *= 20000  # incrase number of DLA
    # mu *= 6.4
    # mu = mean_N_per_cell*(1+bias*deltas)

    mu[~flagged_cells]=0
    #Select cells that will hold a DLA, drawing from the Poisson distribution
    pois = np.random.poisson(mu)#,size=(len(zq),len(mu)))
    #Number of DLAs in each cell (mostly 0, several 1, not many with >1)
    dlas_in_cell = pois*flagged_cells
    ndlas = np.sum(dlas_in_cell)
    #Store information for each of the DLAs that will be added
    dla_z = np.zeros(ndlas)
    dla_skw_id = np.zeros(ndlas, dtype='int64')
    dla_rsd_dz = np.zeros(ndlas)
    dla_count = 0
    dla_z, dla_skw_id, dla_rsd_dz, dla_count = doloop(dlas_in_cell, velocity, zedges, dla_z, dla_skw_id, dla_rsd_dz, dla_count)
    # dla_NHI = get_N(dla_z,Nmin=Nmin,Nmax=Nmax)
    dla_NHI = get_NHI(dla_z,NHI_min=Nmin,NHI_max=Nmax)

    #global id for the skewers
    MOCKIDs = qso['THING_ID'][dla_skw_id]
    ZQSO = zq[dla_skw_id]
    ra = qso['RA'][dla_skw_id]
    dec = qso['DEC'][dla_skw_id]
    z_norsd = qso['Z_noRSD'][dla_skw_id]
    #Make the data into a table HDU
    # dla_table = astropy.table.Table([MOCKIDs,dla_z,dla_rsd_dz,dla_NHI,ZQSO, z_norsd, ra, dec],names=('MOCKID','Z_DLA','DZ_DLA','N_HI_DLA','Z_QSO', 'Z_QSO_NO_RSD', 'RA', 'DEC'))
    # return dla_table, ndlas

    if nhi_low_cut is not None and nhi_high_cut is not None:
        msk = (dla_NHI > nhi_low_cut) & (dla_NHI < nhi_high_cut)
        MOCKIDs = MOCKIDs[msk]
        dla_z = dla_z[msk]
        dla_rsd_dz = dla_rsd_dz[msk]
        dla_NHI = np.ones(len(dla_z))*np.mean([nhi_low_cut, nhi_high_cut])
        ZQSO = ZQSO[msk]
        z_norsd = z_norsd[msk]
        ra = ra[msk]
        dec = dec[msk]

    return [MOCKIDs, dla_z, dla_rsd_dz, dla_NHI, ZQSO, z_norsd, ra, dec], ndlas

######

# Options and main

######

# @profile
def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input_path', type = str, default = None, required = True,
                        help='Path to input directory tree to explore, e.g., /global/cscratch1/sd/*/spectra/*')
    parser.add_argument('--output_path', type = str, default = None, required = True,
                        help='Output path')
    parser.add_argument('--nmin', type = float, default=17.2,
                        help='Minimum value of log(NHI) to consider')
    parser.add_argument('--nmax', type = float, default=22.5,
                        help='Maximum value of log(NHI) to consider')
    parser.add_argument('--nhi-low-cut', type = float, default=None,
                        help='cut HCDs with log(n_HI) < nhi-low_cut')
    parser.add_argument('--nhi-high-cut', type = float, default=None,
                        help='cut HCDs with log(n_HI) > nhi-high_cut')
    parser.add_argument('--dla_bias', type = float, default=2.,
                        help='DLA bias at z=2.25')
    parser.add_argument('--cell_size', type = float, default=2.19,
                        help='size of voxcell')
    parser.add_argument('-seed', type = int, default=None,
                        help='set seed')
    parser.add_argument("-random",type = str, default='False',
                        help="If True, generate randoms")
    parser.add_argument("--random_factor",type = float, default=3.,
                        help="Factor x thus that n_rand = x * n_data")
    args = parser.parse_args()
    
    t0 = time.time()
    random_cond = util.str2bool(args.random)
    if random_cond:
        print("Generating random DLAs...")
    seed = args.seed
    if seed is None:
        seed = np.random.randint(2**31 -1, size=1)[0]
        np.random.seed(seed)
        print("Seed has not been specified. Seed is set to {}".format(seed))
    else:
        np.random.seed(seed)
        print("Specified seed is {}".format(seed))
    
    print("Files will be read from {}".format(args.input_path))
    if os.path.isdir(args.output_path):
        if not random_cond:
            filename = args.output_path + "/dla"
        else:
            filename = args.output_path + "/dla_randoms"
        if args.nhi_low_cut is not None and args.nhi_high_cut is not None:
            filename += "_nhi_{}_{}".format(args.nhi_low_cut, args.nhi_high_cut)
        filename += ".fits"
    else:
        filename = args.output_path
    try:
        outfits = fitsio.FITS(filename, 'rw', clobber=True)
    except IOError:
        print("Can't create or open file {}".format(filename))
        print("Exiting!")
        sys.exit()
    print("Output will be written in {}".format(filename))
    flist = glob.glob(args.input_path+"/*")
    print('Will read', len(flist),' files')
    hdulist = fitsio.FITS(flist[0])
    lam = hdulist[2].read()
    # cosmo_hdu = fitsio.FITS(args.fname_cosmo)[1].read_header()
    z_cell = lam / constant.lya - 1.
    dNdz_arr = dNdz(z_cell, Nmin=args.nmin, Nmax=args.nmax)
    # dNdz_arr *= 3.06  # prov: increase number of DLAs to match observed b_hcd
    # dNdz_arr *= 3.
    # dNdz_arr *= 20000.
    # dNdz_arr *= 6.4
    # dNdz_arr *= 5.9
    if random_cond:
        dNdz_arr *= args.random_factor
    # dNdz_arr /= (-0.01534254*z_cell + 0.0597803)*6.4 / 0.186  # correct the z dependency
    # dz_of_z = dz_of_z_func(args.cell_size)
    dz_of_z = dz_of_z_func(0)  # don't remove DLA close to QSO
    ndlas = 0
    mockid = []
    z_dla = []
    dz_dla = []
    n_hi_dla = []
    z_qso = []
    z_qso_no_rsd = []
    ra = []
    dec = []
    names=['MOCKID','Z_DLA','DZ_DLA','N_HI_DLA','Z_QSO', 'Z_QSO_NO_RSD', 'RA', 'DEC']

    t_loop = time.time()
    for i, fname in enumerate(flist):
        try:
            hdulist = fitsio.FITS(fname)
            table, n = add_DLA_table_to_object_Saclay(hdulist, dNdz_arr,dz_of_z, args.dla_bias, Nmin=args.nmin, Nmax=args.nmax, rand=random_cond, nhi_low_cut=args.nhi_low_cut, nhi_high_cut=args.nhi_high_cut)
            hdulist.close()
        except IOError:
            print("WARNING: can't read {}".format(fname))
        ndlas += n
        # if i==0:
        #     out_table = table
        # else:
        #     out_table = astropy.table.vstack([out_table, table])
        mockid.append(table[0])
        z_dla.append(table[1])
        dz_dla.append(table[2])
        n_hi_dla.append(table[3])
        z_qso.append(table[4])
        z_qso_no_rsd.append(table[5])
        ra.append(table[6])
        dec.append(table[7])
        del table
        gc.collect()
        if i%500==0:
            print('Read %d of %d' %(i,len(flist)))

    print("Loop on files done. Took {} s".format(time.time()-t_loop))
    print("Writting output file ...")
    t_writting = time.time()
    mockid = np.concatenate(mockid)
    z_dla = np.concatenate(z_dla)
    dz_dla = np.concatenate(dz_dla)
    n_hi_dla = np.concatenate(n_hi_dla)
    z_qso = np.concatenate(z_qso)
    z_qso_no_rsd = np.concatenate(z_qso_no_rsd)
    ra = np.concatenate(ra)
    dec = np.concatenate(dec)
    table = [mockid, z_dla, dz_dla, n_hi_dla, z_qso, z_qso_no_rsd, ra, dec]
    outfits.write(table, names=names)
    outfits.close()
    # out_table.write(filename, overwrite=True)
    print("Fits table written. Took {} s".format(time.time() - t_writting))
    print("Draw {} DLAs".format(ndlas))
    print("Took {} s".format(time.time() - t0))


if __name__ == "__main__":
    main()
