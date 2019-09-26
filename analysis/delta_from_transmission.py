import picca.utils
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i")
parser.add_argument("-zcat")
parser.add_argument("-o")
parser.add_argument("--nspec", type=int, default=None)
args = parser.parse_args()

zcat = args.zcat
indir = args.i
outdir = args.o
nspec = args.nspec
print("Args are :\nzcat: {}\nindir: {}\noutdir: {}\nnpec: {}".format(zcat, indir, outdir, nspec))
picca.utils.desi_convert_transmission_to_delta_files(zcat=zcat,indir=indir,outdir=outdir, nspec=nspec)
