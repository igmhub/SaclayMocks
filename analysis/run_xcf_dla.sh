#!/bin/bash -l

#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -p debug
#SBATCH -J do_xcf
#SBATCH -t 00:30:00
#SBATCH -A eboss
#SBATCH --output=logs/xcf_dla.log

outdir=Out/debug_2560_4/
echo "RUNNING XCF"

do_xcf.py --out $outdir/DLA_from_transmission/Correlations/xcf_DLA.fits.gz --in-dir $outdir/DLA_from_transmission/deltas/ --drq $outdir/DLA_from_transmission/Catalogs/DLA.fits.gz &

echo "RUNNING XDMAT"

do_xdmat.py --out $outdir/DLA_from_transmission/Correlations/xdmat_DLA.fits.gz --in-dir $outdir/DLA_from_transmission/deltas/ --drq $outdir/DLA_from_transmission/Catalogs/DLA.fits.gz --rej 0.98 &

wait


echo "RUNNIN EXPORT"

export.py --out $outdir/DLA_from_transmission/Correlations/xcf-exp_DLA.fits.gz --data $outdir/DLA_from_transmission/Correlations/xcf_DLA.fits.gz --dmat $outdir/DLA_from_transmission/Correlations/xdmat_DLA.fits.gz

