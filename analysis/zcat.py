import picca.utils
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-i", type=str)
parser.add_argument("-o", type=str)
parser.add_argument("--downsampling-z-cut", default=None)
parser.add_argument("--downsampling-nb", default=None)
args = parser.parse_args()

path_in = args.i
path_out = args.o
print("Reading from {}".format(path_in))
print("Writting in {}".format(path_out))
print("downsampling_z_cut : {}".format(args.downsampling_z_cut))
print("downsampling_nb : {}".format(args.downsampling_nb))
picca.utils.desi_from_ztarget_to_drq(ztarget=path_in,drq=path_out,spectype="QSO", downsampling_z_cut=args.downsampling_z_cut, downsampling_nb=args.downsampling_nb)
