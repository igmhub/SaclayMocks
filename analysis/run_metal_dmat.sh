#!/bin/bash -l

#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -p regular
#SBATCH -J metal_dmat
#SBATCH -t 02:00:00
#SBATCH -A eboss
#SBATCH --output=logs/metal_dmat_dr14.log3

inpath="Analysis/dr14/DR14Q/"
nproc="--nproc 64"  # comment if you don't want to limit nproc
zmin=0
zmax=10
# metals="SiIII(1207) SiII(1193) SiII(1190) SiII(1260) CIV(1548) CIV(1551)"
metals="SiII(1190) SiII(1193) SiIII(1207) SiII(1260) CIV(eff)"
Om=0.3147
rej=0.999
start=$SECONDS
echo "Begining metal_dmat..."
echo "command: srun -n 1 -c 64 picca_metal_dmat.py --in-dir ${inpath}/deltas/ --out ${inpath}/metal_dmat_z_${zmin}_${zmax}_remove_pairs.fits  --rej $rej --z-cut-min $zmin --z-cut-max $zmax --fid-Om $Om $nproc --abs-igm $metals --remove-same-half-plate-close-pairs"
srun -n 1 -c 64 /usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " picca_metal_dmat.py --in-dir ${inpath}/deltas/ --out ${inpath}/metal_dmat_z_${zmin}_${zmax}_remove_pairs.fits  --rej $rej --z-cut-min $zmin --z-cut-max $zmax --fid-Om $Om $nproc --abs-igm $metals --remove-same-half-plate-close-pairs

# picca_metal_dmat.py --in-dir ${inpath}/deltas/ --out ${inpath}/metal_dmat_z_$zmin_$zmax.fits  --rej $rej --z-cut-min $zmin --z-cut-max $zmax --fid-Om $Om $nproc --abs-igm $metals

# zmin=0
# zmax=2.35
# picca_metal_dmat.py --in-dir ${inpath}/deltas/ --out ${inpath}/metal_dmat_z_${zmin}_${zmax}.fits  --rej $rej --z-cut-min $zmin --z-cut-max $zmax --fid-Om $Om $nproc --abs-igm $metals

# zmin=2.35
# zmax=2.65
# picca_metal_dmat.py --in-dir ${inpath}/deltas/ --out ${inpath}/metal_dmat_z_${zmin}_${zmax}.fits  --rej $rej --z-cut-min $zmin --z-cut-max $zmax --fid-Om $Om $nproc --abs-igm $metals

# zmin=2.65
# zmax=3.05
# picca_metal_dmat.py --in-dir ${inpath}/deltas/ --out ${inpath}/metal_dmat_z_${zmin}_${zmax}.fits  --rej $rej --z-cut-min $zmin --z-cut-max $zmax --fid-Om $Om $nproc --abs-igm $metals

# zmin=3.05
# zmax=10
# picca_metal_dmat.py --in-dir ${inpath}/deltas/ --out ${inpath}/metal_dmat_z_${zmin}_${zmax}.fits  --rej $rej --z-cut-min $zmin --z-cut-max $zmax --fid-Om $Om $nproc --abs-igm $metals

# wait
echo "Done. $(( SECONDS - start )) s\n\n"
