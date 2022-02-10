import picca.converters
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--in-dir")
parser.add_argument("--zcat")
parser.add_argument("--out-dir")
parser.add_argument("--nspec", type=int, default=None)
args = parser.parse_args()

print("Args are :\nzcat: {}\nindir: {}\noutdir: {}\nnpec: {}".format(args.zcat, args.in_dir, args.out_dir, args.nspec))
picca.converters.desi_convert_transmission_to_delta_files(args.zcat, args.out_dir, in_dir=args.in_dir, max_num_spec=args.nspec)
