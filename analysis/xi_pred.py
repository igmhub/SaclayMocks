import fitsio
from SaclayMocks import powerspectrum, util, constant
import numpy as np
import scipy as sp
import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, help="ex: Out/v2.7.1/from_transmission")
parser.add_argument("-o", type=str, help="ex: Out/v2.7.1/from_transmission")
parser.add_argument("--kind", type=str, help="specify if the prediction is Kaiser or FGPA. Default Kaiser", default='Kaiser')
# parser.add_argument("--to-do", type=str, nargs="*", help="ex: cf xcf")
# parser.add_argument("-aa", type=float)
# parser.add_argument("-bb", type=float)
# parser.add_argument("-cc", type=float)
# parser.add_argument("-growthf", type=float)
# parser.add_argument("-sigma", type=float)

args = parser.parse_args()
# a=args.aa
# b=args.bb
# c=args.cc
# G=args.growthf
# sigma_g=args.sigma
# infile = args.i + "/Correlations/e_cf.fits"
# outfile = args.i + "/Correlations/e_cf_pred.fits"

infile = args.i
outfile = args.o
print("Reading {} ...".format(infile))
ecf = fitsio.read(infile, ext=1)
head = fitsio.read_header(infile, ext=1)

zeff = util.zeff(infile)
z0 = 2.2466318099484273
print("zeff = {}".format(zeff))
filename = os.path.expandvars("$SACLAYMOCKS_BASE/etc/params.fits")
a_of_z = util.InterpFitsTable(filename, 'z', 'a')
c_of_z = util.InterpFitsTable(filename, 'z', 'c')
Om = constant.omega_M_0
growthf_24 = util.fgrowth(2.4, Om)
# a = a_of_z.interp(zeff)
a = a_of_z.interp(z0)
# a = 1
b = 1.58
# c = c_of_z.interp(zeff)
c = c_of_z.interp(z0)
G = growthf_24*(1+2.4)/(1+zeff)
sigma_g = 8.2  # sigma = 7.6 mesured at z=2.25
# sigma_g = 2.55
print("a = {}".format(a))
print("b = {}".format(b))
print("c = {}".format(c))
print("G = {}".format(G))
print("sigma_g = {}".format(sigma_g))

### Compute xi pred
# FGPA
if args.kind == 'FGPA':
    xipred = powerspectrum.xi_prediction(a=a,b=b,G=G,sigma_g=sigma_g,c=c)
    rr = np.sqrt(ecf['RP']**2 + ecf['RT']**2)
    mu = ecf['RP'] / rr
    xi = xipred.xi_F(rr, mu)

# Kaiser
if args.kind == 'Kaiser':
    dx = 2.19
    k_ny = np.pi / dx
    rmax = 300
    # fgrowth = util.fgrowth(zeff, constant.omega_M_0)
    fgrowth = 1.
    bias = b*G
    beta = c
    f = beta

    k = np.linspace(0, 10, 100000)
    filename = os.path.expandvars("$SACLAYMOCKS_BASE/etc/PlanckDR12.fits")
    pk = np.float32(powerspectrum.P_0(filename).P(k))
    pk *= np.exp(-k*k*dx*dx/2)**2  # gaussian smoothing
    r, xi = powerspectrum.xi_from_pk(k, pk)
    xi_ham = powerspectrum.xi_Hamilton(r, xi, rmax)

    rp = ecf['RP']
    rt = ecf['RT']
    rr = np.sqrt(rp**2 + rt**2)
    mu = rp / rr
    xi = bias**2 * xi_ham.xi(f, rr, mu)

table = [ecf['RP'], ecf['RT'], ecf['Z'], xi, ecf['CO'], ecf['DM'], ecf['NB']]
names = ['RP', 'RT', 'Z', 'DA', 'CO', 'DM', 'NB']
outfits = fitsio.FITS(outfile, 'rw', clobber=True)
outfits.write(table, names=names, header=head)
outfits[1].write_key('a', a)
outfits[1].write_key('b', b)
outfits[1].write_key('c', c)
outfits[1].write_key('G', G)
outfits[1].write_key('sigma_g', sigma_g)
outfits.close()
print("File {} written.".format(outfile))
