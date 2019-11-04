#!/usr/bin/env python
import argparse
import subprocess
import numpy as np
from SaclayMocks import fit_az
from time import time

'''
This code is used to first fit the parameter a(z) (for a given redshift) on 
the P1D data ; then it tunes the shape of the missing 1D power spectrum
to get the proper 1D power spectrum.

For the moment, you need to compute spectra for each different redshift bin
'''

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--in-dir", type=str, default=None, required=True,
        help="directory where all files are saved")

parser.add_argument("--z", type=float, required=True,
        help="redshift at which a is fitted")

parser.add_argument("--c", type=float, required=True,
        help="value of c parameter, corresponding to beta")

parser.add_argument("--a", type=float, default=0.1, required=False,
        help="reference value of a")

parser.add_argument("--b", type=float, default=1.58, required=False,
        help="value of b")

parser.add_argument("--n-iter", type=int, default=10, required=False,
        help="number of iteration to tune the 1D power spectrum shape")

parser.add_argument("--convergence-factor", type=float, default=1, required=False,
        help="convergence factor to tune p1d shape")

parser.add_argument("--convergence-criterium", type=float, default=None, required=False,
        help="convergence criterium to stop iteration procedure when reached")

parser.add_argument("--seed", type=int, default=42, required=False,
        help="value of b")

parser.add_argument("--compute-spectra", action='store_true', required=False,
        help="Compute spectra files using submit_mocks.py")

parser.add_argument("--fit-az", action='store_true', required=False,
        help="Fit the parameter a(z) using the P1D data")

parser.add_argument("--fit-p1d", action='store_true', required=False,
        help="Tune the shape of the 1D power spectrum")

parser.add_argument("--check-plots", action='store_true', required=False,
        help="Do checking plots")

args = parser.parse_args()

indir = args.in_dir
z = args.z
a = args.a
b = args.b
c = args.c
seed = args.seed
convergence_factor = args.convergence_factor
convergence_criterium = args.convergence_criterium
do_plots = args.check_plots

if args.compute_spectra:
    # Create directories and scripts
    command = 'submit_mocks.py'
    command += ' --mock-dir ' + indir
    command += ' --fit-p1d {} {} {} {} {}'.format(z, a, b, c, seed)
    print("Running {} ...\n".format(command))
    subprocess.check_call(command, shell=True)

    # Compute P1D missing
    command = 'p1d_missing.py '
    command += '--out-file ' + indir + '/mock_0/chunk_1/p1dmiss.fits '
    command += '--beta {} '.format(c)
    command += '--zref {} '.format(z)
    if do_plots:
        command += '--plot-p1d'
    print("\n\nRunning {} ...\n".format(command))
    subprocess.check_call(command, shell=True)

    # Run submit.sh
    command = indir + '/mock_0/output/runs/submit.sh'
    print("\n\nRunning {} ...\n".format(command))
    subprocess.check_call(command, shell=True)

indir += '/mock_0/chunk_1/'
# Fitting a(z)
if args.fit_az:
    print("\n\nFitting a...\n")
    t0 = time()
    fitter = fit_az.Fitter(indir, z, a, c, bb=b, Nreg=1, pixel=0.2,
        convergence_factor=convergence_factor, convergence_criterium=convergence_criterium)
    fitter.read_data()
    fitter.read_mock()
    fitter.minimize()
    fitter.export(indir)
    if do_plots:
        fitter.check_p1d()
    print("\nDone. For z = {}, a = {}. Took {} s\n".format(z, fitter.fit['a'], time() - t0))

# Tuning the P1D shape
if args.fit_p1d:
    print("\nTunning the shape of 1D power spectrum...")
    if not args.fit_az:
        print("Tunning of P1D shape is done using a={}".format(a))
        fitter = fit_az.Fitter(indir, z, a, c, bb=b, Nreg=1, pixel=0.2,
            convergence_factor=convergence_factor, convergence_criterium=convergence_criterium)
        fitter.read_mock()
    else:
        a = fitter.fit['a']
    t0 = time()
    k = np.concatenate((np.arange(0, 3, 0.1), np.arange(3, 20, 0.5)))
    fitter.read_model()
    fitter.read_p1dmiss()
    fitter.compute_p1d(a, bins=k)
    fitter.smooth_p1d()
    print("nspec = {}\n".format(fitter.mock['spectra'].shape))
    if convergence_criterium:
        while not fitter.converged:
            fitter.iterate(a=a, bins=k, plot=do_plots)
    else:
        for n in range(args.n_iter):
            fitter.iterate(a=a, bins=k, plot=do_plots)
    print("\nTunning done. Took {} s".format(time() - t0))
