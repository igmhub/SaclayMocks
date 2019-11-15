#!/bin/bash -l
#SBATCH -N 1
#SBATCH -C haswell
#SBATCH -p debug
#SBATCH -J fit_cf
#SBATCH -t 00:30:00
#SBATCH -A eboss
#SBATCH --output=logs/fit.log

root=/global/cscratch1/sd/tetourne/Out/

# Parameters :
fit_pred=1
sbatch=0
do_deltas=0  # 1 is delta from do_delta, 0 is from transmission

# version=debug_v4.1_2
# version=debug_v4.4
# version=debug_z_dep_qso50
# version=debug_v4.6_38
version=fit_z2.2_check_1

# do_dmat=0  # run dmat only if continuum fitting
# do_export=0  # if run dmat, then run export
# do_xdmat=0
# do_xexport=0
use_dmat=0
use_xdmat=0

fit_cf=1
fit_xcf=0
fit_co=0
object=QSO  # QSO or DLA


zeff=2.2

zbins=0  # set to 1 if you want to fit a particular redshift bin
zmin=2.75
zmax=3.6

if [ $do_deltas -gt 0 ]
then
    inpath=$root/${version}/from_quickquasars/
else
    inpath=$root/${version}/from_transmissions/
fi
if [ $fit_co -gt 0 ]; then
    inpath=$root/${version}/CO/
fi
if [ ! -e ${inpath}/Fit ]; then mkdir ${inpath}/Fit; fi

start=$SECONDS

# ### Run do_dmat.py
# if [ ${do_dmat} -gt 0 ]
# then
#     echo -e "Running do_dmat.py"
#     if [ $sbatch -gt 0 ]; then
# 	command="srun -n 1 -c 64 do_dmat.py --in-dir ${inpath}/Delta_LYA/Delta/ --out ${inpath}/Correlations/dmat.fits  --rej 0.99 --nproc 32"
#     else
# 	command="do_dmat.py --in-dir ${inpath}/Delta_LYA/Delta/ --out ${inpath}/Correlations/dmat.fits  --rej 0.99"
#     fi
#     /usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/do_dmat.log
#     echo -e "Done. $(( SECONDS - start )) s\n\n"
# fi

# ### Run export.py
# if [ ${do_export} -gt 0 ]
# then
#     echo -e "Running export.py"
#     command="export.py --data ${inpath}/Correlations/cf.fits --dmat ${inpath}/Correlations/dmat.fits --out ${inpath}/Correlations/e_cf_dmat.fits"
#     /usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/export_dmat.log
#     echo -e "Done. $(( SECONDS - start )) s\n\n"
# fi

# ### Run do_xdmat.py
# if [ ${do_xdmat} -gt 0 ]
# then
#     echo -e "Running do_xdmat.py"
#     if [ $sbatch -gt 0 ]; then
# 	command="srun -n 1 -c 64 do_xdmat.py --in-dir ${inpath}/Delta_LYA/Delta/ --drq ${inpath}/Catalog/zcat_desi_drq.fits --out ${inpath}/Correlations/xdmat.fits  --rej 0.99 --nproc 64"
#     else
# 	command="do_xdmat.py --in-dir ${inpath}/Delta_LYA/Delta/ --drq ${inpath}/Catalog/zcat_desi_drq.fits --out ${inpath}/Correlations/xdmat.fits  --rej 0.99"
#     fi
#     /usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/do_xdmat.log
#     echo -e "Done. $(( SECONDS - start )) s\n\n"
# fi

# ### Run export.py
# if [ ${do_xexport} -gt 0 ]
# then
#     echo -e "Running export.py"
#     command="export.py --data ${inpath}/Correlations/xcf.fits --dmat ${inpath}/Correlations/xdmat.fits --out ${inpath}/Correlations/e_xcf_dmat.fits"
#     /usr/bin/time -f "%eReal %Uuser %Ssystem %PCPU %M " $command &> ${inpath}/logs/export_xdmat.log
#     echo -e "Done. $(( SECONDS - start )) s\n\n"
# fi

### Generate .ini files
echo -e "Generating .ini files"
# chi2.ini
echo "-- chi2.ini"
ini_files=""
if [ $fit_cf -gt 0 ]
then
    ini_files=${ini_files}"${inpath}/Fit/config_cf.ini "
fi
if [ $fit_xcf -gt 0 ]
then
    ini_files=${ini_files}"${inpath}/Fit/config_xcf.ini "
fi
if [ $fit_co -gt 0 ]
then
    ini_files=${ini_files}"${inpath}/Fit/config_co_${object}.ini "
fi
outfile="result"
if [ $fit_cf -gt 0 ]; then
    outfile=${outfile}"_cf"
fi
if [ $fit_xcf -gt 0 ]; then
    outfile=${outfile}"_xcf"
fi
if [ $fit_co -gt 0 ]; then
    outfile=${outfile}"_co_${object}"
fi
if [ $fit_pred -gt 0 ]; then
    outfile=${outfile}"_pred"
fi
if [ $zbins -gt 0 ]; then
    outfile=${outfile}_z_${zmin}_${zmax}
fi
outfile=${outfile}".h5"
cat > ${inpath}/Fit/chi2.ini << EOF
[data sets]
zeff = $zeff
ini files = ${ini_files}

[fiducial]
filename = PlanckDR12/PlanckDR12.fits

[cosmo-fit type]
cosmo fit func = ap_at

[output]
filename = ${inpath}/Fit/$outfile
EOF

# config_cf.ini
if [ $fit_cf -gt 0 ]
then
    echo "-- config_cf.ini"
    if [ $use_dmat -gt 0 ]; then
	filename=${inpath}/Correlations/e_cf_dmat.fits
    else
	if [ $fit_pred -gt 0 ]; then
	    filename=${inpath}/Correlations/e_cf_pred.fits
	else
	    filename=${inpath}/Correlations/e_cf.fits
	fi
	if [ $zbins -gt 0 ]; then
	    filename=${inpath}/Correlations/e_cf_z_${zmin}_${zmax}.fits
	fi
    fi
    cat > ${inpath}/Fit/config_cf.ini <<EOF
[data]
name = LYA(LYA)xLYA(LYA)
tracer1 = LYA
tracer2 = LYA
tracer1-type = continuous
tracer2-type = continuous
filename = $filename
ell-max = 6

[cuts]
rp-min = 0.
rp-max = 200.

rt-min = 0.
rt-max = 200.

r-min = 10.
r-max = 180.

mu-min = 0.
mu-max = 1.

[model]
model-pk = pk_kaiser
model-xi = xi
z evol LYA = bias_vs_z_std
growth function = growth_factor_no_de
# tranfer-func-mock = data/xi_g_to_xi_F.fits.gz
pk-gauss-smoothing = pk_gauss_smoothing

[parameters]

ap = 1. 0.1 0.5 1.5 fixed
at = 1. 0.1 0.5 1.5 fixed
bao_amp = 1. 0 None None fixed

sigmaNL_per = 0     0 None None fixed
sigmaNL_par = 0 0 None None fixed
1+f         = 1.966    0 None None fixed
growth_rate = 0.966    0 None None fixed

bias_eta_LYA  = -0.14  0.017 None None free
beta_LYA  = 1.8     0.1 None None free
alpha_LYA = 2.9    0   None None fixed

par binsize LYA(LYA)xLYA(LYA) = 4 0.4 0 None fixed
per binsize LYA(LYA)xLYA(LYA) = 4 0.4 0 None fixed

par_sigma_smooth = 3.1 0.1 0 None free
per_sigma_smooth = 3.1 0.1 0 None free
EOF
fi

# config_xcf.ini
if [ $fit_xcf -gt 0 ]
then
    echo "-- config_xcf.ini"
    if [ $use_xdmat -gt 0 ]; then
	filename=${inpath}/Correlations/e_xcf_dmat.fits
    else
	filename=${inpath}/Correlations/e_xcf.fits
	if [ $zbins -gt 0 ]; then
	    filename=${inpath}/Correlations/e_xcf_z_${zmin}_${zmax}.fits
	fi
    fi
    cat > ${inpath}/Fit/config_xcf.ini <<EOF
[data]
name = LYA(LYA)xQSO
tracer1 = QSO
tracer2 = LYA
tracer1-type = discrete
tracer2-type = continuous
filename = $filename
ell-max = 6

[cuts]
rp-min = -200.
rp-max = 200.

rt-min = 0.
rt-max = 200.

r-min = 10.
r-max = 180.

mu-min = -1.
mu-max = 1.

[model]
model-pk = pk_kaiser
model-xi = xi_drp
z evol LYA = bias_vs_z_std
z evol QSO = qso_bias_vs_z_croom
growth function = growth_factor_de
pk-gauss-smoothing = pk_gauss_smoothing
velocity dispersion = pk_velo_lorentz

[parameters]

bias_eta_QSO  = 1. 0. None None fixed
beta_QSO  = 0.3 0.1 None None fixed

croom_par0         = 0.53  0. None None fixed
croom_par1         = 0.289 0. None None fixed
drp_QSO            = 0. 0.1   None None free
sigma_velo_lorentz_QSO = 0. 0.    None None fixed

ap = 1. 0.1 0.5 1.5 fixed
at = 1. 0.1 0.5 1.5 fixed
bao_amp = 1. 0. None None fixed

sigmaNL_per = 0     0. None None fixed
sigmaNL_par = 0 0 None None fixed
# 1+f         = 1.966    0. None None fixed
growth_rate = 0.966 0. None None fixed

bias_eta_LYA  = -0.12  0.017 None None free
beta_LYA  = 1.8    0.1   None None free
alpha_LYA = 2.9     0.   None None fixed

par binsize LYA(LYA)xQSO = 4 0.1 None None fixed
per binsize LYA(LYA)xQSO = 4 0.1 None None fixed

par_sigma_smooth = 3.1 0.4 0 None free
per_sigma_smooth = 3.1 0.4 0 None free

EOF
fi

# config_co.ini
if [ $fit_co -gt 0 ]
then
    filename=${inpath}/Correlations/e_co_${object}.fits
    if [ $zbins -gt 0 ]; then
	filename=${inpath}/Correlations/e_co_z_${zmin}_${zmax}.fits
    fi
    echo -e "-- config_co.ini"
    cat > ${inpath}/Fit/config_co_${object}.ini <<EOF
[data]
name = QSOxQSO
tracer1 = QSO
tracer2 = QSO
tracer1-type = discrete
tracer2-type = discrete
filename = $filename
ell-max = 6

[cuts]
rp-min = -200.
rp-max = 200.

rt-min = 0.
rt-max = 200.

r-min = 20.
r-max = 180.

mu-min = -1.
mu-max = 1.

[model]
model-pk = pk_kaiser
model-xi = xi
z evol QSO = bias_vs_z_std
growth function = growth_factor_de
velocity dispersion = pk_velo_lorentz

[parameters]

bias_eta_QSO  = 1. 0. None None fixed
beta_QSO  = 0.5 0.1 None None free
alpha_QSO = 1.44 0 None None fixed

sigma_velo_lorentz_QSO = 0. 0.1    None None fixed

ap = 1. 0.1 0.5 1.5 free
at = 1. 0.1 0.5 1.5 free
bao_amp = 1. 0. None None fixed

sigmaNL_par = 6.36984     0. None None fixed
sigmaNL_per = 3.24     0. None None fixed
# 1+f         = 1.966    0. None None fixed
growth_rate = 0.962524 0. None None free

par binsize QSOxQSO = 4 0. None None fixed
per binsize QSOxQSO = 4 0. None None fixed
EOF
fi
echo -e "Done. $(( SECONDS - start )) s\n\n"

### Run the fit
echo -e "Running the fit"
if [ $sbatch -gt 0 ]; then
    command="srun -n 1 -c 64 fitter2 ${inpath}/Fit/chi2.ini"
else
    command="picca_fitter2.py ${inpath}/Fit/chi2.ini"
fi
logfile="fit"
if [ $fit_cf -gt 0 ]; then
    logfile=${logfile}"_cf"
fi
if [ $fit_xcf -gt 0 ]; then
    logfile=${logfile}"_xcf"
fi
if [ $fit_co -gt 0 ]; then
    logfile=${logfile}"_co_${object}"
fi
if [ $fit_pred -gt 0 ]; then
    logfile=${logfile}"_pred"
fi
if [ $zbins -gt 0 ]; then
    logfile=${logfile}_z_${zmin}_${zmax}
fi
logfile=${logfile}".log"

$command &> ${inpath}/Fit/$logfile
echo -e "Done. $(( SECONDS - start )) s\n\n"
