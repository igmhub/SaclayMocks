import fitsio
import argparse
import healpy
import scipy as sp


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--to-do", type=str, help='QSO or DLA')
parser.add_argument("--drq", type=str, default=None, required=True, help='path to drq file')
parser.add_argument("--randoms", type=str, default=None, required=True, help='path to randoms file')
parser.add_argument("--out-dir", type=str, default=None, required=True, help='output directory')
parser.add_argument("--nb-data", type=int, default=400000, required=False, help='number of objects')
parser.add_argument("--z-min", type=float, default=None, required=False, help='minimum redshift')
parser.add_argument("--z-max", type=float, default=None, required=False, help='maximum redshift')
parser.add_argument("--nside", type=int, default=16, required=False, help='nside parameter for healpix pixels')
args = parser.parse_args()

nside = args.nside
obj = args.to_do
print("drq is {}".format(args.drq))
print("randoms is {}".format(args.randoms))

### Data
h = fitsio.FITS(args.drq)
ra = h[1]['RA'][:]
dec = h[1]['DEC'][:]
z = h[1]['Z'][:]
phi = ra*sp.pi/180.
th = sp.pi/2.-dec*sp.pi/180.
pix = healpy.ang2pix(nside,th,phi)
print(' There are {} objects in the data catalog'.format(ra.size) )
print(' Across {} different healpix pixels of nside={}'.format(sp.unique(pix).size,nside) )
print(' with redshift in {} < z < {}'.format(z.min(), z.max()))
print('\n')
h.close()

### Randoms
h = fitsio.FITS(args.randoms)
ra = h[1]['RA'][:]
dec = h[1]['DEC'][:]
z = h[1]['Z'][:]
phi = ra*sp.pi/180.
th = sp.pi/2.-dec*sp.pi/180.
pix = healpy.ang2pix(nside,th,phi)
print(' There are {} objects in the random catalog'.format(ra.size) )
print(' Across {} different healpix pixels of nside={}'.format(sp.unique(pix).size,nside) )
print(' with redshift in {} < z < {}'.format(z.min(), z.max()))
print('\n')
h.close()

nbData = args.nb_data
coef = 1
nbRandoms = coef*nbData

### Data
h = fitsio.FITS(args.drq)
data = {}
for k in ['RA','DEC','Z','THING_ID','PLATE','MJD','FIBERID']:
    data[k] = h[1][k][:]
h.close()
w = data['Z']>args.z_min
for k in data.keys():
    data[k] = data[k][w]
w = data['Z']<args.z_max
for k in data.keys():
    data[k] = data[k][w]
phi = data['RA']*sp.pi/180.
th = sp.pi/2.-data['DEC']*sp.pi/180.
pix = healpy.ang2pix(nside,th,phi)
data['PIX'] = pix

### Randoms
h = fitsio.FITS(args.randoms)
rand = {}
for k in ['RA','DEC','Z']:
    rand[k] = h[1][k][:]
for k in ['THING_ID','PLATE','MJD','FIBERID']:
    rand[k] = h[1]['MOCKID'][:]
h.close()
w = rand['Z']>args.z_min
for k in rand.keys():
    rand[k] = rand[k][w]
w = rand['Z']<args.z_max
for k in rand.keys():
    rand[k] = rand[k][w]
phi = rand['RA']*sp.pi/180.
th = sp.pi/2.-rand['DEC']*sp.pi/180.
pix = healpy.ang2pix(nside,th,phi)
rand['PIX'] = pix

### Same HEALpix
w = sp.in1d(data['PIX'],rand['PIX'])
for k,v in data.items():
    v = v[w]
w = sp.in1d(rand['PIX'],data['PIX'])
for k,v in rand.items():
    v = v[w]

### Save data
if nbData > data['RA'].size:
    nbData_tmp = data['RA'].size
    w = sp.random.choice(sp.arange(data['RA'].size), size=nbData_tmp, replace=False)
else:
    w = sp.random.choice(sp.arange(data['RA'].size), size=nbData, replace=False)    
# assert nbData<=data['RA'].size
out = fitsio.FITS(args.out_dir+'/'+obj+'_D_'+str(nbData)+'.fits','rw',clobber=True)
cols = [ v[w] for k,v in data.items() if k not in ['PIX'] ]
names = [ k for k in data.keys() if k not in ['PIX'] ]
out.write(cols,names=names)
out.close()
print(args.out_dir+'/'+obj+'_D_'+str(nbData)+'.fits written')

### Save randoms
if nbRandoms > rand['RA'].size:
    nbRandoms_tmp = coef * nbData_tmp
    w = sp.random.choice(sp.arange(rand['RA'].size), size=nbRandoms_tmp, replace=False)
else:
    w = sp.random.choice(sp.arange(rand['RA'].size), size=nbRandoms, replace=False)
# assert nbRandoms<=rand['RA'].size
out = fitsio.FITS(args.out_dir+'/'+obj+'_R_'+str(nbRandoms)+'.fits','rw',clobber=True)
cols = [ v[w] for k,v in rand.items() if k not in ['PIX'] ]
names = [ k for k in rand.keys() if k not in ['PIX'] ]
out.write(cols,names=names)
out.close()
print(args.out_dir+'/'+obj+'_R_'+str(nbRandoms)+'.fits written')
