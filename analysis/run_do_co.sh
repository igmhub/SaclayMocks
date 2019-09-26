#!/bin/bash -l

#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -p regular
#SBATCH -J do_co
#SBATCH -t 02:00:00
#SBATCH -A eboss
#SBATCH --output=logs/do_co_QSO_v4.5.0.log7

####
# This code only needs the master.fits and master_randoms.fits.gz
####

# Parameters :
# indir=/global/projecta/projectdirs/desi/mocks/lya_forest/develop/saclay/debug/1024_1024_1536/mock_0/output/
# indir=/global/projecta/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.4/v4.4.0/
# indir=/global/projecta/projectdirs/desi/mocks/lya_forest/saclay/v4.1/
indir=/global/cscratch1/sd/tetourne/DesiMocks/v4.5.0/mock_0/output/
# indir=/global/cscratch1/sd/tetourne/DesiMocks/debug/z_dependence_qso4/mock_0/output/

version=v4.5.0_7
nproc=64  # 64 / 4
zref="--z-ref 2.84"  # default in picca is 2.25
rpmin=0
rpmax=200
rtmax=200
np=50
nt=50

todo=QSO  # QSO or DLA
rsd=True

do_drq=1  # drq from master
downsampling_z_cut=10000000
downsampling_nb=0

do_catalogs=1  # catalogs from drq and master_randoms.fits.gz
nb_data=700000
zmin=2.5
zmax=3.6

start=$SECONDS
# Create directories
echo "Creating directories"
if [ ! -e Out/${version}/ ]; then mkdir Out/${version}/; fi
if [ ! -e Out/$version/CO ]; then mkdir Out/$version/CO; fi
inpath=Out/$version/CO
if [ ! -e ${inpath}/Catalog ]; then mkdir ${inpath}/Catalog; fi
if [ ! -e ${inpath}/Correlations ]; then mkdir ${inpath}/Correlations; fi
# if [ ! -e ${inpath}/Correlations/Rands ]; then mkdir ${inpath}/Correlations/Rands; fi
if [ ! -e ${inpath}/logs ]; then mkdir ${inpath}/logs; fi
echo -e "Done. $(( SECONDS - start )) s"
echo -e "Working directory is ${inpath}\n\n"

### Do drq
if [ $do_drq -gt 0 ]; then
    echo "Running zcat_from_master.py"
    if [ $todo == "QSO" ]; then
	command="python zcat_from_master.py -i ${indir}/master.fits -o ${inpath}/Catalog/zcat_desi_drq.fits --to-do QSO --downsampling-z-cut ${downsampling_z_cut} --downsampling-nb ${downsampling_nb} --rsd $rsd"
    elif [ $todo == "DLA" ]; then
	command="python zcat_from_master.py -i ${indir}/master_DLA.fits -o ${inpath}/Catalog/zcat_desi_dla.fits --to-do DLA --downsampling-z-cut ${downsampling_z_cut} --downsampling-nb ${downsampling_nb} --rsd $rsd"
    fi
    echo $command
    /usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/zcat_from_master_$todo.log
    echo -e "Done. $(( SECONDS - start )) s\n\n"
fi

### Do catalogs
if [ $do_catalogs -gt 0 ]; then
    echo "Running randoms_to_drq.py"
    if [ $todo == "QSO" ]; then
    command="python randoms_to_drq.py --to-do QSO --drq $inpath/Catalog/zcat_desi_drq.fits --randoms $indir/master_randoms.fits --out-dir ${inpath}/Catalog/ --nb-data $nb_data --z-min $zmin --z-max $zmax"
    elif [ $todo == "DLA" ]; then
    command="python randoms_to_drq.py --to-do DLA --drq $inpath/Catalog/zcat_desi_dla.fits --randoms $indir/master_DLA_randoms.fits --out-dir ${inpath}/Catalog/ --nb-data $nb_data --z-min $zmin --z-max $zmax"
    fi
    echo $command
    /usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/randoms_to_drq_$todo.log
    echo -e "Done. $(( SECONDS - start )) s\n\n"
fi

### Run do_co.py
echo "Running do_co.py"
echo "DD:"
command="picca_co.py --drq ${inpath}/Catalog/${todo}_D_${nb_data}.fits --out ${inpath}/Correlations/${todo}_DD.fits.gz --type-corr DD --rp-max $rpmax --rt-max $rtmax --rp-min $rpmin --nt $nt --np $np --z-cut-min $zmin --z-cut-max $zmax --nproc $nproc $zref"
echo $command
/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/do_co_${todo}_DD.log &
# echo -e "Done. $(( SECONDS - start )) s\n\n"

echo "RR:"
command="picca_co.py --drq ${inpath}/Catalog/${todo}_R_${nb_data}.fits --out ${inpath}/Correlations/${todo}_RR.fits.gz --type-corr RR --rp-max $rpmax --rt-max $rtmax --rp-min $rpmin --nt $nt --np $np --z-cut-min $zmin --z-cut-max $zmax --nproc $nproc $zref"
echo $command
/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/do_co_${todo}_RR.log &
# echo -e "Done. $(( SECONDS - start )) s\n\n"

echo "DR:"
command="picca_co.py --drq ${inpath}/Catalog/${todo}_D_${nb_data}.fits --drq2 ${inpath}/Catalog/${todo}_R_${nb_data}.fits --out ${inpath}/Correlations/${todo}_DR.fits.gz --type-corr DR --rp-max $rpmax --rt-max $rtmax --rp-min $rpmin --nt $nt --np $np --z-cut-min $zmin --z-cut-max $zmax --nproc $nproc $zref"
echo $command
/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/do_co_${todo}_DR.log &
# echo -e "Done. $(( SECONDS - start )) s\n\n" 

echo "RD:"
command="picca_co.py --drq ${inpath}/Catalog/${todo}_R_${nb_data}.fits --drq2 ${inpath}/Catalog/${todo}_D_${nb_data}.fits --out ${inpath}/Correlations/${todo}_RD.fits.gz --type-corr DR --rp-max $rpmax --rt-max $rtmax --rp-min $rpmin --nt $nt --np $np --z-cut-min $zmin --z-cut-max $zmax --nproc $nproc $zref"
echo $command
/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/do_co_${todo}_RD.log &

wait
echo -e "Done. $(( SECONDS - start )) s\n\n"

### Run export.py
echo "Running export.py"
command="picca_export_co.py --DD-file ${inpath}/Correlations/${todo}_DD.fits.gz --RR-file ${inpath}/Correlations/${todo}_RR.fits.gz --RD-file ${inpath}/Correlations/${todo}_RD.fits.gz --DR-file ${inpath}/Correlations/${todo}_DR.fits.gz --out ${inpath}/Correlations/e_co_${todo}.fits"
echo $command
/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/export_co_${todo}.log
echo -e "Done. $(( SECONDS - start )) s\n\n"
