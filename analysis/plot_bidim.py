import matplotlib.pyplot as plt
import fitsio
import h5py
import scipy as sp
from SaclayMocks import util
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, help="ex: Out/v2.7.1/from_transmission")
parser.add_argument("--to-do", type=str, nargs="*", help="ex: cf xcf")
parser.add_argument("--title", default="")
args = parser.parse_args()

indir = args.i
if 'cf' in args.to_do:
    cf_file = 'e_cf.fits'
    fit_file = 'result_cf.h5'
    nrp = 50
    nrt = 50
    fit_name = 'LYA(LYA)xLYA(LYA)'
if 'xcf' in args.to_do:
    cf_file = 'e_xcf.fits'
    fit_file = 'result_xcf.h5'
    nrp = 100
    nrt = 50
    fit_name = 'LYA(LYA)xQSO'

cf = fitsio.FITS(indir+"/Correlations/"+cf_file)
fit = h5py.File(indir+"/Fit/"+fit_file)

rp = cf[1].read()['RP']
rt = cf[1].read()['RT']
da = cf[1].read()['DA']
co = cf[1].read()['CO']
da_2d = util.convert1DTo2D(da, nrp, nrt)
co_2d = util.convert1DTo2D(sp.diag(co), nrp, nrt)
da_fit = fit[fit_name]['fit'].value
da_fit_2d = util.convert1DTo2D(da_fit, nrp, nrt)
extent = [rt.min(), rt.max(), rp.min(), rp.max()]

plt.imshow((da_2d - da_fit_2d) / sp.sqrt(co_2d), origin='lower',cmap='seismic', vmin=-2., vmax=2., extent=extent)
plt.colorbar()
plt.xlabel(r'$rt [Mpc.h^{-1}]$')
plt.ylabel(r'$rp [Mpc.h^{-1}]$')
plt.title(args.title)
plt.show()
