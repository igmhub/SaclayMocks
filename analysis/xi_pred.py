import fitsio
from SaclayMocks import powerspectrum, util, constant
import numpy as np
from functools import partial
import scipy as sp
import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str, help="input file")
parser.add_argument("-o", type=str, help="output file")
parser.add_argument("--kind", type=str, help="specify if the prediction is Kaiser or FGPA. Default FGPA", default='FGPA')
parser.add_argument("--zeff", type=float, help="specify a given redshift. If not specify, takes zeff from input file.", default=None)
# parser.add_argument("--to-do", type=str, nargs="*", help="ex: cf xcf")
parser.add_argument("--a", type=float, default=None)
parser.add_argument("--b", type=float, default=None)
parser.add_argument("--c", type=float, default=None)
parser.add_argument("--growthf", type=float, default=None)
parser.add_argument("--sigma", type=float, default=None)
parser.add_argument("--p1d-file", type=str, default=None)
parser.add_argument("--no-proj", action='store_true', help='do not do the projection with the distorsion matrix')
# parser.add_argument("--cf-1d", action='store_true', help='return the 1D correlation function')

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

filename = os.path.expandvars("$SACLAYMOCKS_BASE/etc/params.fits")
a_of_z = util.InterpFitsTable(filename, 'z', 'a')
b_of_z = util.InterpFitsTable(filename, 'z', 'b')
c_of_z = util.InterpFitsTable(filename, 'z', 'c')

# Read parameters
zeff = args.zeff
if zeff is None:
    zeff = util.zeff(infile)

a = args.a
if a is None:
    a = a_of_z.interp(zeff)

b = args.b
if b is None:
    b = b_of_z.interp(zeff)

c = args.c
if c is None:
    c = c_of_z.interp(zeff)

G = args.growthf
if G is None:
    Om = constant.omega_M_0
    growthf_24 = util.fgrowth(2.4, Om)
    G = growthf_24*(1+2.4)/(1+zeff)

sigma_g = args.sigma
if sigma_g is None:
    if args.p1d_file is not None:
        print("sigma_g is computed using {}".format(args.p1d_file))
        try:
            data_p1d = fitsio.read(args.p1d_file, ext=1)
            p1dmiss = sp.interpolate.interp1d(data_p1d['k'], data_p1d['P1DmissRSD'])
        except:
            p1dmiss = partial(util.InterpP1Dmissing(args.p1d_file), redshift=zeff)
        sigma_g = util.sigma_g(zeff, p1dmiss=p1dmiss, c=c)
    else:
        print("sigma_g is computed using standard p1dmissing")
        sigma_g = util.sigma_g(zeff, c=c)


### Compute xi pred
# FGPA
if args.kind == 'FGPA':
    print("Computing FGPA model...")
    print("zeff = {}".format(zeff))
    print("a = {}".format(a))
    print("b = {}".format(b))
    print("c = {}".format(c))
    print("G = {}".format(G))
    print("sigma_g = {}".format(sigma_g))
    xipred = powerspectrum.xi_prediction(a=a,b=b,G=G,sigma_g=sigma_g,c=c)
    rr = np.sqrt(ecf['RP']**2 + ecf['RT']**2)
    mu = ecf['RP'] / rr
    xi = xipred.xi_F(rr, mu)
    if args.no_proj is not None:
        np.dot(xi, ecf['DM'], out=xi)

# Kaiser
if args.kind == 'Kaiser':
    print("Computing Kaiser model...")
    binsize = 4  # binsize of correlation function
    dx = 2.19
    k_ny = np.pi / dx
    rmax = 300
    # fgrowth = util.fgrowth(zeff, constant.omega_M_0)
    # fgrowth = 1.
    bias = b
    beta = c
    # f = beta
    print("bias = {}".format(bias))
    print("G = {}".format(G))
    print("beta = {}".format(beta))

    k = np.linspace(0, 10, 100000)
    filename = os.path.expandvars("$SACLAYMOCKS_BASE/etc/PlanckDR12.fits")
    pk = np.float32(powerspectrum.P_0(filename).P(k))
    # pk *= np.exp(-k*k*dx*dx/2)**2  # gaussian smoothing effect
    # pk *= np.sinc(-k*dx/(2*np.pi))**2  # voxel size effect
    pk *= np.sinc(-k*binsize/(2*np.pi))**2  # binsize of correlation function effect
    r, xi = powerspectrum.xi_from_pk(k, pk)
    xi_ham = powerspectrum.xi_Hamilton(r, xi, rmax)

    rp = ecf['RP']
    rt = ecf['RT']
    rr = np.sqrt(rp**2 + rt**2)
    mu = rp / rr
    xi = bias**2 * G**2 * xi_ham.xi(beta, rr, mu)

z = np.ones_like(ecf['Z'])*zeff
table = [ecf['RP'], ecf['RT'], z, xi, ecf['CO'], ecf['DM'], ecf['NB']]
names = ['RP', 'RT', 'Z', 'DA', 'CO', 'DM', 'NB']
#table = [ecf['RP'], ecf['RT'], z, xi, ecf['NB']]
#names = ['RP', 'RT', 'Z', 'DA', 'NB']
outfits = fitsio.FITS(outfile, 'rw', clobber=True)
outfits.write(table, names=names, header=head)
outfits[1].write_key('a', str(a))
outfits[1].write_key('b', str(b))
outfits[1].write_key('c', str(c))
outfits[1].write_key('G', str(G))
outfits[1].write_key('sigma_g', str(sigma_g))
outfits.close()
print("File {} written.".format(outfile))
