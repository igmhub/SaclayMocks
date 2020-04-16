import numpy as np
import fitsio
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str)
parser.add_argument("-zeff", type=float)
args = parser.parse_args()

# Read data
print("Reading {}".format(args.i))
data1 = fitsio.read(args.i, ext=1)
head1 = fitsio.read_header(args.i, ext=1)

redshift = np.ones_like(data1['Z']) * args.zeff
table1 = [data1['RP'], data1['RT'], redshift, data1['DA'], data1['CO'], data1['DM'], data1['NB']]
names1 = ['RP', 'RT', 'Z', 'DA', 'CO', 'DM', 'NB']

try:
    data2 = fitsio.read(args.i, ext=2)
    table2 = [data2['DMRP'], data2['DMRT'], redshift]
    names2 = ['DMRP', 'DMRT', 'DMZ']
except OSError:
    table2 = [data1['RP'], data1['RT'], redshift]
    names2 = ['DMRP', 'DMRT', 'DMZ']

# Make a save of the input CF
command = "cp --backup=numbered {i} {i}.save".format(i=args.i)
subprocess.check_call(command, shell=True)
print("Input file copied to {}.save".format(args.i))

# Write the new CF
fits = fitsio.FITS(args.i, 'rw', clobber=True)
fits.write(table1, names=names1, header=head1)
fits.write(table2, names=names2)
print("New file {} written with zeff = {}".format(args.i, args.zeff))
