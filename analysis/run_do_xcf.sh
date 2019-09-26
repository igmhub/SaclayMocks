#!/bin/bash -l
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -p regular
#SBATCH -J do_xcf
#SBATCH -t 08:00:00
#SBATCH -A eboss
#SBATCH --output=logs/do_xcf_v4.6_38_1.log

# OMP_NUM_THREADS=1
sbatch=1
# Parameters :
# indir=DesiMocks/debug/1024-1024-1536/mock_0/output/
# indir=DesiMocks/debug/256-256-1536_iso/mock_1/output/
# indir=DesiMocks/prod/2560-2560-1536_v3.2.0/mock_0/output/
# indir=DesiMocks/prod/512-512-1536_v2.7/mock_1/output/
# indir=/global/cscratch1/sd/tetourne/DesiMocks/v4.6/mock_1/output/
# indir=/global/cscratch1/sd/legoff/DesiMocks/prod/256/mock_0/output/
indir=/global/cscratch1/sd/tetourne/DesiMocks/debug2/v4.6_38/mock_1/output/
# indir=/global/cscratch1/sd/tetourne/DesiMocks/v4.7/mock_1/output/

quick_folder=quick-0.1/
version=debug_v4.6_38
# version=v4.7_1

zmin=0
zmax=2.25
rpmax=200
rtmax=200
Om=0.3147

shuffle_seed=42

downsampling_z_cut=2.1
downsampling_nb=1000000
# nspec="--nspec 700000"  # comment if you don't want to set npec
nproc="--nproc 32"
compute_deltas=0
delta_from_do_deltas=0  # 0 is delta from transmission; use no-project in this case
no_project="--no-project"  # comment if you don't want to use --no-project option
use_xdmat=0
do_xdmat=0

coadd=1

if [ $coadd -gt 0 ]; then
    exp="e_xcf_z_${zmin}_${zmax}"
    xcf="xcf_z_${zmin}_${zmax}"
    xdmat="xdmat_z_${zmin}_${zmax}"
    zmin_delta=0
    zmax_delta=10
else
    exp="e_xcf"
    cf="xcf"
    dmat="xdmat"
    zmin_delta=$zmin
    zmin_delta=$zmax
fi

start=$SECONDS
### Create directories
echo "Creating directories"
if [ ! -e Out/${version}/ ]; then mkdir Out/${version}/; fi
if [ $delta_from_do_deltas -gt 0 ]
then
    if [ ! -e Out/${version}/from_quickquasars ]; then mkdir Out/${version}/from_quickquasars; fi
    inpath=Out/${version}/from_quickquasars
else
    if [ ! -e Out/${version}/from_transmissions ]; then mkdir Out/${version}/from_transmissions; fi
    inpath=Out/${version}/from_transmissions
fi
if [ ! -e ${inpath}/Catalog ]; then mkdir ${inpath}/Catalog; fi
if [ ! -e ${inpath}/Correlations ]; then mkdir ${inpath}/Correlations; fi
if [ ! -e ${inpath}/Delta_LYA ]; then mkdir ${inpath}/Delta_LYA; fi
if [ ! -e ${inpath}/Delta_LYA/Delta ]; then mkdir ${inpath}/Delta_LYA/Delta; fi
if [ ! -e ${inpath}/logs ]; then mkdir ${inpath}/logs; fi
echo -e "Done. $(( SECONDS - start )) s"
echo -e "Working directory is ${inpath}\n\n"

if [ $compute_deltas -gt 0 ]
then
    if [ $delta_from_do_deltas -gt 0 ]
    then
	### Run zcat.py
	echo "Running zcat.py"
	command="python zcat.py -i ${indir}/${quick_folder}/zcat.fits -o ${inpath}/Catalog/zcat_desi_drq.fits  --downsampling-z-cut ${downsampling_z_cut} --downsampling-nb ${downsampling_nb}"
	echo $command
	/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/zcat.log
	echo -e "Done. $(( SECONDS - start )) s\n\n"
	### Run do_deltas.py
	echo "Running do_deltas.py"
	if [ $sbatch -gt 0 ]; then
	    command="srun -n 1 -c 64 picca_deltas.py --drq ${inpath}/Catalog/zcat_desi_drq.fits --in-dir ${indir}/${quick_folder}/spectra-16/ --out-dir ${inpath}/Delta_LYA/Delta/ --mode desi --iter-out-prefix Delta_LYA/Log/delta_attributes --zqso-min $zmin_delta --zqso-max $zmax_delta $nspec $nproc"
	else
	    command="picca_deltas.py --drq ${inpath}/Catalog/zcat_desi_drq.fits --in-dir ${indir}/${quick_folder}/spectra-16/ --out-dir ${inpath}/Delta_LYA/Delta/ --mode desi --iter-out-prefix Delta_LYA/Log/delta_attributes --zqso-min $zmin_delta --zqso-max $zmax_delta $nspec"
	fi
	echo $command
	/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/do_deltas.log
	echo -e "Done. $(( SECONDS - start )) s\n\n"

    else
	# Run zcat_from_master.py
	echo "Running zcat_from_master.py"
	command="python zcat_from_master.py -i ${indir}/master.fits -o ${inpath}/Catalog/zcat_desi_drq.fits --downsampling-z-cut ${downsampling_z_cut} --downsampling-nb ${downsampling_nb}"
	echo $command
	/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/zcat_from_master.log
	echo -e "Done. $(( SECONDS - start )) s\n\n"
	# # Run zcat_from_master.py
	# echo "Running zcat_from_master.py"
	# command="python zcat_from_master.py -i ${indir}/master_randoms.fits -o ${inpath}/Catalog/zcat_desi_drq_rand.fits --downsampling-z-cut ${downsampling_z_cut} --downsampling-nb ${downsampling_nb} --rand True"
	# echo $command
	# /usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/zcat_from_master.log
	# echo -e "Done. $(( SECONDS - start )) s\n\n"

	### Run delta_from_transmission.py
	echo "Running delta_from_transmission.py"
	command="python delta_from_transmission.py -i ${indir} -zcat ${inpath}/Catalog/zcat_desi_drq.fits -o ${inpath}/Delta_LYA/Delta/ $nspec"
	echo $command
	/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/delta_from_transmission.log
	echo -e "Done. $(( SECONDS - start )) s\n\n"
    fi
fi

### Run do_xcf.py
echo "Running do_xcf.py :"
command="picca_xcf.py --in-dir ${inpath}/Delta_LYA/Delta/ --drq ${inpath}/Catalog/zcat_desi_drq.fits --out ${inpath}/Correlations/${xcf}.fits --rp-max $rpmax --rt-max $rtmax --z-min-obj $zmin --z-max-obj $zmax $nspec $no_project $nproc --fid-Om $Om --z-evol-obj 1.44"
echo $command
/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/${xcf}.log &

echo "Running do_xcf.py :"
command="picca_xcf.py --in-dir ${inpath}/Delta_LYA/Delta/ --drq ${inpath}/Catalog/zcat_desi_drq.fits --out ${inpath}/Correlations/${xcf}_shuffle_$shuffle_seed.fits --rp-max $rpmax --rt-max $rtmax --z-min-obj $zmin --z-max-obj $zmax $nspec $no_project $nproc --fid-Om $Om --z-evol-obj 1.44 --shuffle-distrib-obj-seed $shuffle_seed"
echo $command
/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/${xcf}_shuffle.log &

wait
echo -e "Done. $(( SECONDS - start )) s\n\n"

### Run do_xdmat.py
if [ ${do_xdmat} -gt 0 ]
then
    echo -e "Running do_xdmat.py"
    if [ $sbatch -gt 0 ]; then
	command="srun -n 1 -c 64 picca_xdmat.py --in-dir ${inpath}/Delta_LYA/Delta/ --drq ${inpath}/Catalog/zcat_desi_drq.fits --out ${inpath}/Correlations/${xdmat}.fits  --rej 0.98 $nproc"
    else
	command="picca_xdmat.py --in-dir ${inpath}/Delta_LYA/Delta/ --drq ${inpath}/Catalog/zcat_desi_drq.fits --out ${inpath}/Correlations/${xdmat}.fits  --rej 0.98 $nproc"
    fi
    echo $command
    /usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/${xdmat}.log
    echo -e "Done. $(( SECONDS - start )) s\n\n"
fi

# ### Run export.py
# if [ ${do_xexport} -gt 0 ]
# then
#     echo -e "Running export.py"
#     command="export.py --data ${inpath}/Correlations/xcf.fits --dmat ${inpath}/Correlations/xdmat.fits --out ${inpath}/Correlations/e_xcf_dmat.fits"
#     /usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/export_xdmat.log
#     echo -e "Done. $(( SECONDS - start )) s\n\n"
# fi

### Run export.py
echo "Running export.py"
if [ $use_xdmat -gt 0 ]; then
    command="picca_export.py --data ${inpath}/Correlations/${xcf}.fits --dmat ${inpath}/Correlations/${xdmat}.fits --out ${inpath}/Correlations/${exp}.fits --remove-shuffled-correlation $inpath/Correlations/${xcf}_shuffle_$shuffle_seed.fits"
else
    command="picca_export.py --data ${inpath}/Correlations/${xcf}.fits --out ${inpath}/Correlations/${exp}.fits --remove-shuffled-correlation $inpath/Correlations/${xcf}_shuffle_$shuffle_seed.fits"
fi
/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/${exp}.log
echo -e "Done. $(( SECONDS - start )) s\n\n"
