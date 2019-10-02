#!/usr/bin/env python
import os, sys
import argparse
import subprocess
from SaclayMocks import fit_az


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

parser.add_argument("--seed", type=int, default=42, required=False,
        help="value of b")

parser.add_argument("--compute-spectra", action='store_true', required=False,
        help="Compute spectra files using submit_mocks.py")

parser.add_argument("--check-plots", action='store_true', required=False,
        help="Do checking plots")

args = parser.parse_args()

indir = args.in_dir
z = args.z
a = args.a
b = args.b
c = args.c
seed = args.seed
do_plots = args.check_plots

if args.compute_spectra:
    # Create directories and scripts
    command = 'submit_mocks.py'
    command += ' --mock-dir ' + indir
    command += ' --fit-p1d {} {} {} {} {}'.format(z, a, b, c, seed)
    print("Running {} ...".format(command))
    subprocess.check_call(command, shell=True)

    # Compute P1D missing
    command = 'p1d_missing.py '
    command += '--out-file ' + indir + '/mock_0/chunk_1/p1dmiss.fits '
    command += '--beta {} '.format(c)
    if do_plots:
        command += '--plot-p1d'
    print("Running {} ...".format(command))
    subprocess.check_call(command, shell=True)

    # Run submit.sh
    command = indir + '/mock_0/output/runs/submit.sh'
    print("Running {} ...".format(command))
    subprocess.check_call(command, shell=True)

print("Fitting a...")
indir += '/mock_0/chunk_1/'
fit_az = fit_az.Fitter(indir, z, a, c, bb=b, Nreg=1, xx=100., pixel=0.2)
fit_az.read_data()
fit_az.read_mock()
fit_az.minimize()
fit_az.export(indir)
if do_plots:
    fit_az.check_p1d()
