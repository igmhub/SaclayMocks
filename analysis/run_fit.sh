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
fit_pred=0
sbatch=0
do_deltas=1  # 1 is delta from do_delta, 0 is from transmission

# version=debug_v4.1_2
# version=debug_v4.4
# version=debug_z_dep_qso50
# version=debug_v4.6_38
# version=fit_z1.8
# version=dr16_paper
version=v4.7.22_masked_dla20.5_4

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
fit_metal=0


zeff=2.85

hesse=0  # set to 1 to print correlations between parameters
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

[hesse]
level = $hesse

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
    if [ $fit_metal -gt 0 ]; then
	metals="
[metals]
filename = ${inpath}/Correlations/metal_dmat.fits
model-pk-met = pk_kaiser
model-xi-met = xi
z evol = bias_vs_z_std
# in tracer1 = SiII(1260) NV(1243) NV(1239) SiIII(1207) NI(1200) SiII(1193) SiII(1190)
# in tracer2 = SiII(1260) NV(1243) NV(1239) SiIII(1207) NI(1200) SiII(1193) SiII(1190)
in tracer1 = CIV(eff) SiII(1260) SiIII(1207) SiII(1193) SiII(1190)
in tracer2 = CIV(eff) SiII(1260) SiIII(1207) SiII(1193) SiII(1190)
"
	metal_pars="
bias_eta_SiII(1260) = -0.002205 0.01 None 0. free
beta_SiII(1260) = 0.5 0. None None fixed
alpha_SiII(1260) = 1.0 0. None None fixed

bias_eta_SiIII(1207) = -0.004535 0.01 None 0. free
beta_SiIII(1207) = 0.5 0. None None fixed
alpha_SiIII(1207) = 1.0 0. None None fixed

bias_eta_SiII(1193) = -0.00209 0.01 None 0. free
beta_SiII(1193) = 0.5 0. None None fixed
alpha_SiII(1193) = 1.0 0. None None fixed

bias_eta_SiII(1190) = -0.00296 0.01 None 0. free
beta_SiII(1190) = 0.5 0. None None fixed
alpha_SiII(1190) = 1.0 0. None None fixed

# bias_eta_NV(1243) = -0.000 0.01 None 0. free
# beta_NV(1243) = 0.5 0. None None fixed
# alpha_NV(1243) = 1.0 0. None None fixed

# bias_eta_NV(1239) = -0.000 0.01 None 0. free
# beta_NV(1239) = 0.5 0. None None fixed
# alpha_NV(1239) = 1.0 0. None None fixed

# bias_eta_NI(1200) = -0.000 0.01 None 0. free
# beta_NI(1200) = 0.5 0. None None fixed
# alpha_NI(1200) = 1.0 0. None None fixed

bias_eta_CIV(eff) = -0.005156549561399881 0.001 None 0. free
beta_CIV(eff) = 0.27 0.01 None 1. fixed
alpha_CIV(eff) = 1. 0.01 None None fixed
"
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
rp-min = -200.
rp-max = 200.

rt-min = 0.
rt-max = 200.

r-min = 10.
r-max = 180.

mu-min = -1.
mu-max = 1.

[model]
# model-pk = pk_kaiser
model-pk = pk_hcd_Rogers2018
model-xi = xi
z evol LYA = bias_vs_z_std
# growth function = growth_factor_no_de
growth function = growth_factor_de
# tranfer-func-mock = data/xi_g_to_xi_F.fits.gz
pk-gauss-smoothing = pk_gauss_smoothing
# small scale nl = dnl_arinyo

${metals}

[parameters]
${metal_pars}
ap = 1. 0.1 0.5 1.5 fixed
at = 1. 0.1 0.5 1.5 fixed
bao_amp = 1. 0 None None fixed

sigmaNL_per = 0     0 None None fixed
sigmaNL_par = 0 0 None None fixed
# sigmaNL_per = 3.24     0 None None fixed
# sigmaNL_par = 6.36984 0 None None fixed
# growth_rate = 0.966    0 None None fixed
growth_rate = 0.970386193694752 0. None None fixed

bias_eta_LYA  = -0.20  0.02 None None free
beta_LYA  = 1.8     0.1 None None free
alpha_LYA = 2.9    0   None None fixed

bias_hcd = -0.05222070235530851 0.1 None 0. free
beta_hcd = 0.6098209629393987 0.1 None None free
L0_hcd = 10. 1. None None fixed

# dnl_arinyo_q1 = 0.8558 0.1 None None fixed
# dnl_arinyo_kv = 1.11454 0.1 None None fixed
# dnl_arinyo_av = 0.5378 0.1 None None fixed
# dnl_arinyo_bv = 1.607 0.1 None None fixed
# dnl_arinyo_kp = 19.47 0.1 None None fixed

# BB-LYA(LYA)xLYA(LYA)-0-broadband_sky-scale-sky = 0.009411087303413272 0.1 None None free
# BB-LYA(LYA)xLYA(LYA)-0-broadband_sky-sigma-sky = 31.41897384749371 0.1 None None free

par binsize LYA(LYA)xLYA(LYA) = 4 0.4 0 None fixed
per binsize LYA(LYA)xLYA(LYA) = 4 0.4 0 None fixed

par_sigma_smooth = 3.1 0.1 0 None free
per_sigma_smooth = 3.1 0.1 0 None free

# [broadband]
# bb1 = add pre rp,rt 0:0:1 0:0:1 broadband_sky

[priors]
beta_hcd = gaussian 0.5 0.09
# bias_eta_CIV(eff) = gaussian -0.005 0.0026

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
