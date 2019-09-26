#!/bin/bash -l

#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -p regular
#SBATCH -J do_cf
#SBATCH -t 04:00:00
#SBATCH -A eboss
#SBATCH --output=logs/do_cf_v4.7_1_3.log

# OMP_NUM_THREADS=1
sbatch=1
# Parameters :
# indir=DesiMocks/debug/1024-1024-1536/mock_70/output/
# indir=/global/projecta/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.1/
# indir=/global/projecta/projectdirs/desi/mocks/lya_forest/develop/london/v8.0/v8.0.0/
# indir=/global/cscratch1/sd/tetourne/DesiMocks/v4.6/mock_1/output/
# indir=/global/cscratch1/sd/tetourne/DesiMocks/debug/v4.0/mock_0/output/
# indir=/global/cscratch1/sd/tetourne/DesiMocks/debug/4705a65/mock_0/output/
# indir=/global/cscratch1/sd/tetourne/DesiMocks/debug2/dc339d72/mock_0/output/
# indir=/global/cscratch1/sd/tetourne/DesiMocks/debug2/v4.6_38/mock_1/output/
# indir=/global/cscratch1/sd/legoff/DesiMocks/prod/256/mock_0/output/
indir=/global/cscratch1/sd/tetourne/DesiMocks/v4.7/mock_1/output/

quick_folder=quick-0.2/
# version=debug_1024_71
# version=debug_256_17
# version=debug_1024_71
# version=debug_v4.6_38
version=v4.7_1

Om=0.3147
zmin=2.45
zmax=2.75
rpmax=200
rtmax=200
downsampling_z_cut_min=0
downsampling_z_cut_max=10
downsampling_nb=1000000
# nspec="--nspec 20000"  # comment if you don't want to set npec
nproc="--nproc 32"  # comment if you don't want to limit number of proc
compute_deltas=0
delta_from_do_deltas=0  # 0 is delta from transmission; use no-project in this case
no_project="--no-project"  # comment if you don't want to use --no-project option
do_cf=1
do_dmat=0
use_dmat=0
do_exp=1

coadd=1

if [ $coadd -gt 0 ]; then
    exp="e_cf_z_${zmin}_${zmax}"
    cf="cf_z_${zmin}_${zmax}"
    dmat="dmat_z_${zmin}_${zmax}"
    zmin_delta=0
    zmax_delta=10
else
    exp="e_cf"
    cf="cf"
    dmat="dmat"
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

### Run codes
if [ $compute_deltas -gt 0 ]
then
    if [ $delta_from_do_deltas -gt 0 ]
    then
	# Run zcat.py
	echo "Running zcat.py"
	command="python zcat.py -i ${indir}/${quick_folder}/zcat.fits -o ${inpath}/Catalog/zcat_desi_drq.fits --downsampling-z-cut ${downsampling_z_cut_min} --downsampling-z-cut-max ${downsampling_z_cut_max} --downsampling-nb ${downsampling_nb}"
	echo $command
	/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/zcat.log
	echo -e "Done. $(( SECONDS - start )) s\n\n"

	# Run do_deltas.py
	echo "Running do_deltas.py"
	if [ $sbatch -gt 0 ]; then
	    command="srun -n 1 -c 64 picca_deltas.py --drq ${inpath}/Catalog/zcat_desi_drq.fits --in-dir ${indir}/${quick_folder}/spectra-16/ --out-dir ${inpath}/Delta_LYA/Delta/ --mode desi --zqso-min $zmin_delta --zqso-max $zmax_delta $nspec $nproc"
	else
	    command="picca_deltas.py --drq ${inpath}/Catalog/zcat_desi_drq.fits --in-dir ${indir}/${quick_folder}/spectra-16/ --out-dir ${inpath}/Delta_LYA/Delta/ --mode desi --zqso-min $zmin_delta --zqso-max $zmax_delta $nspec $nproc"
	fi
	echo $command
	/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/do_deltas.log
	echo -e "Done. $(( SECONDS - start )) s\n\n"

    else
	# Run zcat_from_master.py
	echo "Running zcat_from_master.py"
	command="python zcat_from_master.py -i ${indir}/master.fits -o ${inpath}/Catalog/zcat_desi_drq.fits --downsampling-z-cut ${downsampling_z_cut_min} --downsampling-z-cut-max ${downsampling_z_cut_max} --downsampling-nb ${downsampling_nb} "
	echo $command
	/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/zcat_from_master.log
	echo -e "Done. $(( SECONDS - start )) s\n\n"
	
	# Run delta_from_transmission.py
	echo "Running delta_from_transmission.py"
	command="python delta_from_transmission.py -i ${indir} -zcat ${inpath}/Catalog/zcat_desi_drq.fits -o ${inpath}/Delta_LYA/Delta/ $nspec"
	echo $command
	/usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/delta_from_transmission.log
	echo -e "Done. $(( SECONDS - start )) s\n\n"
    fi
fi

# Run do_cf.py
if [ ${do_cf} -gt 0 ]
then
    echo "Running do_cf.py :"
    if [ $sbatch -gt 0 ]; then
        command="srun -n 1 -c 64 picca_cf.py --in-dir ${inpath}/Delta_LYA/Delta/ --out ${inpath}/Correlations/${cf}.fits --rp-max $rpmax --rt-max $rtmax $nspec ${no_project} --z-cut-min $zmin --z-cut-max $zmax --fid-Om $Om $nproc"
    else
        command="picca_cf.py --in-dir ${inpath}/Delta_LYA/Delta/ --out ${inpath}/Correlations/${cf}.fits --rp-max $rpmax --rt-max $rtmax $nspec ${no_project} --z-cut-min $zmin --z-cut-max $zmax --fid-Om $Om $nproc"
    fi
    echo $command
    /usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/${cf}.log
    echo -e "Done. $(( SECONDS - start )) s\n\n"
fi

### Run do_dmat.py
if [ ${do_dmat} -gt 0 ]
then
    echo -e "Running do_dmat.py"
    if [ $sbatch -gt 0 ]; then
	command="srun -n 1 -c 64 picca_dmat.py --in-dir ${inpath}/Delta_LYA/Delta/ --out ${inpath}/Correlations/${dmat}.fits  --rej 0.99 --z-cut-min $zmin --z-cut-max $zmax --fid-Om $Om $nproc"
    else
	command="picca_dmat.py --in-dir ${inpath}/Delta_LYA/Delta/ --out ${inpath}/Correlations/${dmat}.fits  --rej 0.99 --z-cut-min $zmin --z-cut-max $zmax --fid-Om $Om $nproc"
    fi
    echo $command
    /usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/${dmat}.log
    echo -e "Done. $(( SECONDS - start )) s\n\n"
fi

# Run export.py
if [ ${do_exp} -gt 0 ]
then
    echo "Running export.py"
    if [ $sbatch -gt 0 ]; then
        if [ $use_dmat -gt 0 ]; then
    	command="picca_export.py --data ${inpath}/Correlations/${cf}.fits --dmat ${inpath}/Correlations/${dmat}.fits --out ${inpath}/Correlations/${exp}.fits"
        else
    	command="srun -n 1 -c 64 picca_export.py --data ${inpath}/Correlations/${cf}.fits --out ${inpath}/Correlations/${exp}.fits"
        fi
    else
        if [ $use_dmat -gt 0 ]; then
    	command="picca_export.py --data ${inpath}/Correlations/${cf}.fits --dmat ${inpath}/Correlations/${dmat}.fits --out ${inpath}/Correlations/${exp}.fits"
        else
    	command="picca_export.py --data ${inpath}/Correlations/${cf}.fits --out ${inpath}/Correlations/${exp}.fits"
        fi
    fi
    echo $command
    /usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/${exp}.log
    echo -e "Done. $(( SECONDS - start )) s\n\n"
fi
