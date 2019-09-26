#!/bin/bash -l

#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -p regular
#SBATCH -J metal_dmat
#SBATCH -t 02:00:00
#SBATCH -A eboss
#SBATCH --output=logs/metal_xdmat_dr14.log

inpath="Analysis/dr14/DR14Q/"
nproc="--nproc 64"  # comment if you don't want to limit nproc
zmin=0
zmax=10
# metals="SiIII(1207) SiII(1193) SiII(1190) SiII(1260) CIV(1548) CIV(1551)"
metals="SiII(1190) SiII(1193) SiIII(1207) SiII(1260)"
Om=0.3147
rej=0.999
start=$SECONDS
echo "Begining metal_dmat..."
command="srun -n 1 -c 64 picca_metal_xdmat.py --in-dir ${inpath}/deltas/ --drq /global/projecta/projectdirs/sdss/data/sdss/dr14/eboss/qso/DR14Q/DR14Q_v4_4.fits --out ${inpath}/metal_xdmat_z_${zmin}_${zmax}.fits  --rej $rej --z-cut-min $zmin --z-cut-max $zmax --fid-Om $Om $nproc --abs-igm $metals --z-evol-obj 1.44"
echo "command: $command"
$command

echo "Done. $(( SECONDS - start )) s\n\n"
