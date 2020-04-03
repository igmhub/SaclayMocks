import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from SaclayMocks import util


SMALL_SIZE = 16
MEDIUM_SIZE = 20
BIGGER_SIZE = 22
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize

plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.rc('figure', figsize=(9,7))

### Options
plot_bias = True
plot_beta = True
plot_bias_eta = True
plot_beff = True
plot_b_hcd = True

# toplot = ['v4.7.22_raw', 'v4.7.22-0.0', 'v4.7.22-0.2_nhi20.3', 'v4.7.22-0.2_fvoigt_v4.7.22', 'redo_dr16', 'dr16_fvoigt_v4.7.22']
# toplot = ['v4.7.22-0.2_nhi20.3', 'v4.7.22-0.2_nhi20.5', 'v4.7.22-0.2_dndz3_nhi20.3', 'v4.7.22-0.2_fvoigt_v4.7.22', 'redo_dr16', 'dr16_fvoigt_v4.7.22']
# toplot = ['v4.7.22-0.2_nhi20.3', 'v4.7.22-0.2_nhi20.5', 'v4.7.22-0.2_dndz3_nhi20.3', 'redo_dr16']
# toplot = ['v4.7.22-0.2_nhi20.3', 'v4.7.22-0.2_nhi20.5', 'v4.7.22-0.2_dndz3_nhi20.3', 'v4.7.22-0.2_dndz3_nhi20.3_fixed_lya', 'redo_dr16']
# toplot = ['v4.7.22-0.2_nhi20.3', 'v4.7.22-0.2_fvoigt_v4.7.22', 'redo_dr16', 'dr16_fvoigt_v4.7.22']
# toplot = ['v4.7.22_raw', 'v4.7.22-0.0_bis', 'v4.7.22-0.2_nhi20.3', 'redo_dr16']
# toplot = ['redo_dr16', 'dr16_fvoigt_v4.7.22', 'dr16_mask_fvoigt_v4.7.22']
# toplot = ['v4.7.22-0.0', 'v4.7.22-0.0_bis']
toplot = ['redo_dr16']

# labels = toplot
# labels = ['mock_Rogers2018', 'mock_no_mask', 'DR16_Rogers2018', 'DR16_no_mask', 'dr16_fvoigt_v4.7.22']
# labels = ['nhi_20.3', 'nhi_20.5', '3*dndz_nhi_20.3', 'dr16']
# labels = ['nhi_20.3', 'nhi_20.5', '3*dndz_nhi_20.3', '3*dndz_nhi_20.3_fixed_lya', 'dr16']
# labels = ['raw mocks', 'mock-0.0', 'mock-0.2', 'DR16']
labels = ['DR16']
# labels = ['DR16_mask_Roger', 'DR16_Guy', 'DR16_mask_Guy']
# labels = ['eboss-0.0', 'eboss-0.0_seed126429']

mean_items = {}  # mocks used to compute means
mean_items['dla_mean'] = ['v4.7.22_dla', 'v4.7.27_dla']
mean_items['dla_coadd_mean'] = ['v4.7.22_dla_coadd', 'v4.7.27_dla_coadd']
mean_items['mock_mean'] = ['v4.7.22', 'v4.7.27']
mean_items['mock_coadd_mean'] = ['v4.7.22_coadd', 'v4.7.27_coadd']

plot_pred1 = False
plot_pred2 = False
plot_shades = False

correct_z_dep = False
gamma = {'bias_eta':2.2, 'beta':-1.8, 'bias':3.8, 'beff':3.2}

absolute_bias = True

def func(x, a, gamma):
    return a*(1+x)**gamma

z0 = 2.4

# create dictionnaries
files = {}
redshift = {}
bias_eta = {}
bias_eta_err = {}
beta = {}
beta_err = {}
bias = {}
bias_err = {}
beff = {}
beff_err = {}
b_hcd = {}
b_hcd_err = {}
cov = {}
colors = {}
p_bias = {}
p_beta = {}
p_bias_eta = {}
p_beff = {}

# Choose colors
colors['data_helion'] = 'black'
colors['redo_dr16'] = 'royalblue'
colors['dr16_fvoigt_v4.7.22'] = 'darkgreen'
colors['dr16_mask_fvoigt_v4.7.22'] = 'red'
colors['v4.7.22_raw'] = 'green'
colors['v4.7.22-0.0'] = 'royalblue'
colors['v4.7.22-0.0_bis'] = 'royalblue'
colors['v4.7.22-0.2_nhi20.3'] = 'r'
colors['v4.7.22-0.2_nhi20.5'] = 'darkgreen'
colors['v4.7.22-0.2_dndz3_nhi20.3'] = 'r'
colors['v4.7.22-0.2_dndz3_nhi20.3_fixed_lya'] = 'darkorange'
colors['v4.7.22-0.2_fvoigt_v4.7.22'] = 'darkorange'
colors['v4.7.22-0.2_fvoigt_exp'] = 'magenta'
# colors['v4.7.27_coadd'] = 'darkblue'
# colors['mock_mean'] = 'royalblue'
# colors['mock_coadd_mean'] = 'royalblue'
# colors['v4.7.22_dla'] = 'darkorange'
# colors['v4.7.22_dla_coadd'] = 'darkorange'
# colors['v4.7.27_dla'] = 'darkred'
# colors['v4.7.27_dla_coadd'] = 'darkred'
# colors['dla_mean'] = 'r'
# colors['dla_coadd_mean'] = 'r'
# colors['v4.7.22_dla_nomask'] = 'magenta'
# colors['v4.7.22_raw'] = 'green'
# colors['v4.7.22_raw_coadd'] = 'green'

### Inputs .h5 files :
# v4.7.22 raw mocks
files['v4.7.22_raw'] = ["/global/cscratch1/sd/tetourne/Out/v4.7.22_1/from_transmissions/Fit/result_cf.h5",
                        "/global/cscratch1/sd/tetourne/Out/v4.7.22_2/from_transmissions/Fit/result_cf.h5",
                        "/global/cscratch1/sd/tetourne/Out/v4.7.22_3/from_transmissions/Fit/result_cf.h5",
                        "/global/cscratch1/sd/tetourne/Out/v4.7.22_4/from_transmissions/Fit/result_cf.h5"]

# v4.7.22 eboss-0.0
files['v4.7.22-0.0'] = ["/global/cscratch1/sd/tetourne/Out/v4.7.22_1/from_quickquasars/Fit/result_cf.h5",
                        "/global/cscratch1/sd/tetourne/Out/v4.7.22_2/from_quickquasars/Fit/result_cf.h5",
                        "/global/cscratch1/sd/tetourne/Out/v4.7.22_3/from_quickquasars/Fit/result_cf.h5",
                        "/global/cscratch1/sd/tetourne/Out/v4.7.22_4/from_quickquasars/Fit/result_cf.h5"]

# v4.7.22 eboss-0.0
files['v4.7.22-0.0_bis'] = ["/global/projecta/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output/eboss-0.0_seed126429/cf_z_0_2.35-exp.h5",
                            "/global/projecta/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output/eboss-0.0_seed126429/cf_z_2.35_2.65-exp.h5",
                            "/global/projecta/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output/eboss-0.0_seed126429/cf_z_2.65_3.05-exp.h5",
                            "/global/projecta/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output/eboss-0.0_seed126429/cf_z_3.05_10-exp.h5"]

# v4.7.22 eboss-0.2 : DLAs masked with log(n_HI) > 20.3
files['v4.7.22-0.2_nhi20.3'] = ["/global/cscratch1/sd/tetourne/Out/v4.7.22_masked_dla20.3_1/from_quickquasars/Fit/result_cf.h5",
                                    "/global/cscratch1/sd/tetourne/Out/v4.7.22_masked_dla20.3_2/from_quickquasars/Fit/result_cf.h5",
                                    "/global/cscratch1/sd/tetourne/Out/v4.7.22_masked_dla20.3_3/from_quickquasars/Fit/result_cf.h5",
                                    "/global/cscratch1/sd/tetourne/Out/v4.7.22_masked_dla20.3_4/from_quickquasars/Fit/result_cf.h5"]

# v4.7.22 eboss-0.2 : DLAs masked with log(n_HI) > 20.5
files['v4.7.22-0.2_nhi20.5'] = ["/global/cscratch1/sd/tetourne/Out/v4.7.22_masked_dla20.5_1/from_quickquasars/Fit/result_cf.h5",
                                    "/global/cscratch1/sd/tetourne/Out/v4.7.22_masked_dla20.5_2/from_quickquasars/Fit/result_cf.h5",
                                    "/global/cscratch1/sd/tetourne/Out/v4.7.22_masked_dla20.5_3/from_quickquasars/Fit/result_cf.h5",
                                    "/global/cscratch1/sd/tetourne/Out/v4.7.22_masked_dla20.5_4/from_quickquasars/Fit/result_cf.h5"]

# v4.7.22 eboss-0.2 with 3 times more DLAs and cut with log(n_HI) > 20.3 :
files['v4.7.22-0.2_dndz3_nhi20.3'] = ["/global/cscratch1/sd/tetourne/Out/v4.7.22_dndz3_masked_dla20.3_1/from_quickquasars/Fit/result_cf.h5",
                                "/global/cscratch1/sd/tetourne/Out/v4.7.22_dndz3_masked_dla20.3_2/from_quickquasars/Fit/result_cf.h5",
                                "/global/cscratch1/sd/tetourne/Out/v4.7.22_dndz3_masked_dla20.3_3/from_quickquasars/Fit/result_cf.h5",
                                "/global/cscratch1/sd/tetourne/Out/v4.7.22_dndz3_masked_dla20.3_4/from_quickquasars/Fit/result_cf.h5"]

# v4.7.22 eboss-0.2 with 3 times more DLAs and cut with log(n_HI) > 20.3 and fixed LYA parameters :
files['v4.7.22-0.2_dndz3_nhi20.3_fixed_lya'] = ["/global/cscratch1/sd/tetourne/Out/v4.7.22_dndz3_masked_dla20.3_1/from_quickquasars/Fit/result_cf_fixed_lya.h5",
                                                "/global/cscratch1/sd/tetourne/Out/v4.7.22_dndz3_masked_dla20.3_2/from_quickquasars/Fit/result_cf_fixed_lya.h5",
                                                "/global/cscratch1/sd/tetourne/Out/v4.7.22_dndz3_masked_dla20.3_3/from_quickquasars/Fit/result_cf_fixed_lya.h5",
                                                "/global/cscratch1/sd/tetourne/Out/v4.7.22_dndz3_masked_dla20.3_4/from_quickquasars/Fit/result_cf_fixed_lya.h5"]

# v4.7.22 eboss-0.2 : no DLA masking, treaded with Fvoigt_v4.7.22.txt
files['v4.7.22-0.2_fvoigt_v4.7.22'] = ["/global/cscratch1/sd/tetourne/Out/v4.7.22_dla_1/from_quickquasars/Fit/result_cf.h5",
                                      "/global/cscratch1/sd/tetourne/Out/v4.7.22_dla_2/from_quickquasars/Fit/result_cf.h5",
                                      "/global/cscratch1/sd/tetourne/Out/v4.7.22_dla_3/from_quickquasars/Fit/result_cf.h5",
                                      "/global/cscratch1/sd/tetourne/Out/v4.7.22_dla_4/from_quickquasars/Fit/result_cf.h5"]

# v4.7.27 eboss-0.0
files['v4.7.27-0.0'] = ["/global/cscratch1/sd/tetourne/Out/v4.7.27_1/from_quickquasars/Fit/result_cf.h5",
                        "/global/cscratch1/sd/tetourne/Out/v4.7.27_2/from_quickquasars/Fit/result_cf.h5",
                        "/global/cscratch1/sd/tetourne/Out/v4.7.27_3/from_quickquasars/Fit/result_cf.h5",
                        "/global/cscratch1/sd/tetourne/Out/v4.7.27_4/from_quickquasars/Fit/result_cf.h5"]

# DR16 data analysis from Helion
files['data_helion'] = ["/global/u1/h/hdumasde/Run_programs/igmhub/picca/Data/eBOSS/picca_DR16_paper_analysis/Correlations/Fit/cf_z_0_2.35/result.h5",
                        "/global/u1/h/hdumasde/Run_programs/igmhub/picca/Data/eBOSS/picca_DR16_paper_analysis/Correlations/Fit/cf_z_2.35_2.65/result.h5",
                        "/global/u1/h/hdumasde/Run_programs/igmhub/picca/Data/eBOSS/picca_DR16_paper_analysis/Correlations/Fit/cf_z_2.65_3.05/result.h5",
                        "/global/u1/h/hdumasde/Run_programs/igmhub/picca/Data/eBOSS/picca_DR16_paper_analysis/Correlations/Fit/cf_z_3.05_10/result.h5"]

# DR16 data with standard analysis (as in DR16 paper)
files['redo_dr16'] = ["/global/cscratch1/sd/tetourne/Analysis/redo_dr16/Fits/cf/Kaiser_sky_met_hcd/z_0_2.35/result.h5",
                      "/global/cscratch1/sd/tetourne/Analysis/redo_dr16/Fits/cf/Kaiser_sky_met_hcd/z_2.35_2.65/result.h5",
                      "/global/cscratch1/sd/tetourne/Analysis/redo_dr16/Fits/cf/Kaiser_sky_met_hcd/z_2.65_3.05/result.h5",
                      "/global/cscratch1/sd/tetourne/Analysis/redo_dr16/Fits/cf/Kaiser_sky_met_hcd/z_3.05_10/result.h5"]

# DR16 data without DLA masking, treated with Fvoigt_v4.7.22.txt
files['dr16_fvoigt_v4.7.22'] = ["/global/cscratch1/sd/tetourne/Analysis/dr16_no_dla_masking/Fits/cf/Kaiser_sky_met_hcd/z_0_2.35/result_fvoigt_v4.7.22.h5",
                                "/global/cscratch1/sd/tetourne/Analysis/dr16_no_dla_masking/Fits/cf/Kaiser_sky_met_hcd/z_2.35_2.65/result_fvoigt_v4.7.22.h5",
                                "/global/cscratch1/sd/tetourne/Analysis/dr16_no_dla_masking/Fits/cf/Kaiser_sky_met_hcd/z_2.65_3.05/result_fvoigt_v4.7.22.h5",
                                "/global/cscratch1/sd/tetourne/Analysis/dr16_no_dla_masking/Fits/cf/Kaiser_sky_met_hcd/z_3.05_10/result_fvoigt_v4.7.22.h5"]

# DR16 data with DLA masked and fitted with HCD model from Julien Edmond
files['dr16_mask_fvoigt_v4.7.22'] = ["/global/cscratch1/sd/tetourne/Analysis/redo_dr16/Fits/cf/Kaiser_sky_met_hcd/z_0_2.35/result_fvoigt_v4.7.22_highcut_nhi20.3.h5",
                                     "/global/cscratch1/sd/tetourne/Analysis/redo_dr16/Fits/cf/Kaiser_sky_met_hcd/z_2.35_2.65/result_fvoigt_v4.7.22_highcut_nhi20.3.h5",
                                     "/global/cscratch1/sd/tetourne/Analysis/redo_dr16/Fits/cf/Kaiser_sky_met_hcd/z_2.65_3.05/result_fvoigt_v4.7.22_highcut_nhi20.3.h5",
                                     "/global/cscratch1/sd/tetourne/Analysis/redo_dr16/Fits/cf/Kaiser_sky_met_hcd/z_3.05_10/result_fvoigt_v4.7.22_highcut_nhi20.3.h5"]

# DR16 data without DLA masking, treated with Fvoigt_exp.txt
files['dr16_fvoigt_exp'] = ["/global/cscratch1/sd/tetourne/Analysis/dr16_no_dla_masking/Fits/cf/Kaiser_sky_met_hcd/z_0_2.35/result_fvoigt_exp.h5",
                            "/global/cscratch1/sd/tetourne/Analysis/dr16_no_dla_masking/Fits/cf/Kaiser_sky_met_hcd/z_2.35_2.65/result_fvoigt_exp.h5",
                            "/global/cscratch1/sd/tetourne/Analysis/dr16_no_dla_masking/Fits/cf/Kaiser_sky_met_hcd/z_2.65_3.05/result_fvoigt_exp.h5",
                            "/global/cscratch1/sd/tetourne/Analysis/dr16_no_dla_masking/Fits/cf/Kaiser_sky_met_hcd/z_3.05_10/result_fvoigt_exp.h5"]


### Input : all the various data sets
## les mocks avec quickquasars (distorsion matrix) + DLAs (v4.7.22_eboss-0.2_thomas)
# redshift['v4.7.22_dla_nomask'] = np.array([2.09, 2.21, 2.52, 2.85])
# bias_eta['v4.7.22_dla_nomask'] = np.array([-0.172, -0.195, -0.226, -0.267])
# bias_eta_err['v4.7.22_dla_nomask'] = np.array([0.003, 0.004, 0.007, 0.017])
# beta['v4.7.22_dla_nomask'] = np.array([1.55, 1.49, 1.07, 1.11])
# beta_err['v4.7.22_dla_nomask'] = np.array([0.07, 0.07, 0.07, 0.16])
# cor['v4.7.22_dla_nomask'] = -0.87

## les mocks avec quickquasars + masked DLAs log(n_HI) > 20 (v4.7.22_eboss-0.2)
# redshift['v4.7.22_dla'] = np.array([2.09, 2.21, 2.52, 2.85])
# bias_eta['v4.7.22_dla'] = np.array([-0.174, -0.195, -0.227, -0.279])
# bias_eta_err['v4.7.22_dla'] = np.array([0.003, 0.004, 0.006, 0.018])
# beta['v4.7.22_dla'] = np.array([1.69, 1.58, 1.16, 1.29])
# beta_err['v4.7.22_dla'] = np.array([0.08, 0.07, 0.07, 0.21])
# cor['v4.7.22_dla'] = -0.87
# redshift['v4.7.22_dla'] = np.array([2.09, 2.21, 2.52, 2.85])
# bias_eta['v4.7.22_dla'] = np.array([-0.17834206335604685, -0.19824135015363673, -0.23347086797901706, -0.26844055267731287])
# bias_eta_err['v4.7.22_dla'] = np.array([0.0023817769778462423, 0.002676011550310295, 0.0049626231393716065, 0.012809679373192641])
# beta['v4.7.22_dla'] = np.array([1.6972452265815168, 1.571389214992836, 1.1929697421092473, 1.0865762569014479])
# beta_err['v4.7.22_dla'] = np.array([0.06321656925591052, 0.05876429477516389, 0.042423648043276974, 0.11576246791024232])
# cor['v4.7.22_dla'] = -0.87

# # le coadd de v4.7.22_eboss-0.2
# redshift['v4.7.22_dla_coadd'] = 2.20
# bias_eta['v4.7.22_dla_coadd'] = -0.1895
# bias_eta_err['v4.7.22_dla_coadd'] = 0.0029
# beta['v4.7.22_dla_coadd'] = 1.57
# beta_err['v4.7.22_dla_coadd'] = 0.06
# cor['v4.7.22_dla_coadd'] = -0.90

# # les mocks qvec quickquasars + masked DLAs with log(n_HI) > 20 (v4.7.27_eboss-0.2)
# redshift['v4.7.27_dla'] = np.array([2.09, 2.21, 2.52, 2.85])
# bias_eta['v4.7.27_dla'] = np.array([-0.169, -0.196, -0.235, -0.277])
# bias_eta_err['v4.7.27_dla'] = np.array([0.003, 0.003, 0.007, 0.018])
# beta['v4.7.27_dla'] = np.array([1.56, 1.61, 1.32, 1.03])
# beta_err['v4.7.27_dla'] = np.array([0.07, 0.07, 0.09, 0.15])
# cor['v4.7.27_dla'] = -0.87

# # le coadd de v4.7.27_eboss-0.2
# redshift['v4.7.27_dla_coadd'] = 2.20
# bias_eta['v4.7.27_dla_coadd'] = -0.1876
# bias_eta_err['v4.7.27_dla_coadd'] = 0.0028
# beta['v4.7.27_dla_coadd'] = 1.54
# beta_err['v4.7.27_dla_coadd'] = 0.05
# cor['v4.7.27_dla_coadd'] = -0.89

# # les mocks avec quickquasars sans DLAs (v4.7.22_eboss-0.0)
# redshift['v4.7.22'] = np.array([2.09, 2.21, 2.52, 2.85])
# bias_eta['v4.7.22'] = np.array([-0.1802, -0.1987, -0.246, -0.270])
# bias_eta_err['v4.7.22'] = np.array([0.0019, 0.0021, 0.004, 0.010])
# beta['v4.7.22'] = np.array([1.76, 1.553, 1.44, 1.06])
# beta_err['v4.7.22'] = np.array([0.03, 0.031, 0.04, 0.06])
# cor['v4.7.22'] = -0.97

# # le coadd de v4.7.22_eboss-0.0
# redshift['v4.7.22_coadd'] = 2.20
# bias_eta['v4.7.22_coadd'] = -0.1942
# bias_eta_err['v4.7.22_coadd'] = 0.0015
# beta['v4.7.22_coadd'] = 1.564
# beta_err['v4.7.22_coadd'] = 0.022
# cor['v4.7.22_coadd'] = -0.96

# # les mocks avec quickquasars sans DLAs (v4.7.27_eboss-0.0)
# redshift['v4.7.27'] = np.array([2.09, 2.21, 2.52, 2.85])
# bias_eta['v4.7.27'] = np.array([-0.1777, -0.1956, -0.241, -0.282])
# bias_eta_err['v4.7.27'] = np.array([0.0018, 0.0021, 0.004, 0.010])
# beta['v4.7.27'] = np.array([1.69, 1.525, 1.38, 1.14])
# beta_err['v4.7.27'] = np.array([0.03, 0.030, 0.04, 0.07])
# cor['v4.7.27'] = -0.965

# # le coadd de v4.7.27_eboss-0.0
# redshift['v4.7.27_coadd'] = 2.20
# bias_eta['v4.7.27_coadd'] = -0.1903
# bias_eta_err['v4.7.27_coadd'] = 0.0015
# beta['v4.7.27_coadd'] = 1.502
# beta_err['v4.7.27_coadd'] = 0.021
# cor['v4.7.27_coadd'] = -0.96

# # les mocks directement sur les transmissions (v4.7.22)
# redshift['v4.7.22_raw'] = np.array([2.10, 2.24, 2.53, 2.87])
# bias_eta['v4.7.22_raw'] = np.array([-0.1786, -0.2054, -0.244, -0.277])
# bias_eta_err['v4.7.22_raw'] = np.array([0.002, 0.0022, 0.004, 0.013])
# beta['v4.7.22_raw'] = np.array([1.73, 1.58, 1.39, 1.07])
# beta_err['v4.7.22_raw'] = np.array([0.04, 0.03, 0.04, 0.09])
# cor['v4.7.22_raw'] = -0.96

# # le coadd de v4.7.22 raw mocks
# redshift['v4.7.22_raw_coadd'] = 2.25
# bias_eta['v4.7.22_raw_coadd'] = -0.2034
# bias_eta_err['v4.7.22_raw_coadd'] = 0.0016
# beta['v4.7.22_raw_coadd'] = 1.569
# beta_err['v4.7.22_raw_coadd'] = 0.023
# cor['v4.7.22_raw_coadd'] = -0.96

### # les donnees
# data_helion_bias_eta= np.array( [[ 2.13781251, -0.18558474,  0.00650885],
#                           [ 2.27259149, -0.20089904,  0.00566513],
#                           [ 2.54744584, -0.23754895,  0.00859699],
#                           [ 2.9263245,  -0.28984325,  0.01907289]] )
# data_helion_beta = np.array([[2.13781251, 1.71518188, 0.14670333],
#                       [2.27259149, 1.63599596, 0.12377548],
#                       [2.54744584, 1.34603853, 0.10685724],
#                       [2.9263245,  1.20116875, 0.16525412]])
# redshift['data_helion'] = data_helion_beta[:,0]
# bias_eta['data_helion'] = data_helion_bias_eta[:,1]
# bias_eta_err['data_helion'] = data_helion_bias_eta[:,2]
# beta['data_helion'] = data_helion_beta[:,1]
# beta_err['data_helion'] = data_helion_beta[:,2]
# cor['data_helion'] = -0.9

# Mon analyse dr16
# redshift['redo_dr16'] = np.array([2.136, 2.276, 2.551, 2.914])
# bias_eta['redo_dr16'] = np.array([-0.1795977211203897, -0.19377778657236358, -0.2236876829818197, -0.2928514012382003])
# bias_eta_err['redo_dr16'] = np.array([0.005769053847991895, 0.00526939132784233, 0.008421287322167362, 0.018740373333761263])
# beta['redo_dr16'] = np.array([2.093809306317262, 1.7112686242434043, 1.4274033068601393, 1.264422978848562])
# beta_err['redo_dr16'] = np.array([0.2104492929009926, 0.13322965498374356, 0.13841133835520955, 0.19412583530626168])
# cor['redo_dr16'] = -0.875

# DR16 without DLA masking and Fvoigt profile v4.7.22
# redshift['dr16_fvoigt_v4.7.22'] = np.array([2.135, 2.275, 2.551, 2.914])
# bias_eta['dr16_fvoigt_v4.7.22'] = np.array([-0.16747388922748288, -0.18802498580322063, -0.21458917298426639, -0.30164555366354134])
# bias_eta_err['dr16_fvoigt_v4.7.22'] = np.array([0.0070574051927350815, 0.006319810140388751, 0.009660827961796433, 0.02049020057648377])
# beta['dr16_fvoigt_v4.7.22'] = np.array([2.770406496234253, 2.189148675345433, 1.650081980537978, 1.4357382321776464])
# beta_err['dr16_fvoigt_v4.7.22'] = np.array([0.5202757717928898, 0.28690340773842793, 0.24132001880155676, 0.29638415622576236])
# cor['dr16_fvoigt_v4.7.22'] = np.array([-0.70, -0.75, -0.80, -0.82])

# # Le fit de la pred, faites dans chaque bin de v4.7.22_raw sur 10 < r < 180:
# zpred1 = np.array([2.10, 2.24, 2.53, 2.87])
# bias_eta_pred1 = np.array([-0.1753, -0.2019, -0.237, -0.284])
# beta_pred1 = np.array([1.68, 1.57, 1.37, 1.18])
# bias_pred1 = bias_eta_pred1 * 0.97 / beta_pred1
# beff_pred1 = bias_pred1 * np.sqrt(1+2/3*beta_pred1+1/5*beta_pred1**2)

# # Le fit de la pred, faites dans chaque bin de v4.7.22_raw, sur 40 < r < 180 :
# zpred2 = np.array([2.10, 2.24, 2.53, 2.87])
# bias_eta_pred2 = np.array([-0.180, -0.205, -0.240, -0.290])
# beta_pred2 = np.array([1.76, 1.63, 1.40, 1.19])
# bias_pred2 = bias_eta_pred2 * 0.97 / beta_pred2
# beff_pred2 = bias_pred2 * np.sqrt(1+2/3*beta_pred2+1/5*beta_pred2**2)

# # Compute means on mocks
# for h in toplot:
#     if 'mean' in h:
#         redshift[h] = np.array([redshift[m] for m in mean_items[h]]).mean(axis=0)
#         bias_eta[h] = np.array([bias_eta[m] for m in mean_items[h]]).mean(axis=0)  # assume same errors
#         bias_eta_err[h] = np.array([bias_eta_err[m] for m in mean_items[h]]).mean(axis=0)
#         bias_eta_err[h] /= np.sqrt(len(mean_items[h]))
#         beta[h] = np.array([beta[m] for m in mean_items[h]]).mean(axis=0)
#         beta_err[h] = np.array([beta_err[m] for m in mean_items[h]]).mean(axis=0)
#         beta_err[h] /= np.sqrt(len(mean_items[h]))
#         cor[h] = np.mean([cor[m] for m in mean_items[h]])

### Compute bias, b_eff, and fit with power law
for item in toplot:
    redshift[item] = np.zeros(len(files[item]))
    bias_eta[item] = np.zeros(len(files[item]))
    bias_eta_err[item] = np.zeros(len(files[item]))
    beta[item] = np.zeros(len(files[item]))
    beta_err[item] = np.zeros(len(files[item]))
    b_hcd[item] = np.zeros(len(files[item]))
    b_hcd_err[item] = np.zeros(len(files[item]))
    cov[item] = np.zeros(len(files[item]))
    for i, f in enumerate(files[item]):
        pars = util.extract_h5file(f)
        redshift[item][i] = pars[2]['zeff']
        if 'bias_eta_LYA' in pars[2]:
            bias_eta[item][i] = pars[2]['bias_eta_LYA']
            bias_eta_err[item][i] = pars[3]['bias_eta_LYA']
            cov[item][i] = pars[4]['cov[beta_LYA, bias_eta_LYA]']
        else:
            bias_eta[item][i] = pars[2]['bias_LYA']
            bias_eta_err[item][i] = pars[3]['bias_LYA']
            cov[item][i] = pars[4]['cov[beta_LYA, bias_LYA]']
        beta[item][i] = pars[2]['beta_LYA']
        beta_err[item][i] = pars[3]['beta_LYA']
        if 'bias_hcd' in pars[2]:
            b_hcd[item][i] = pars[2]['bias_hcd']
            b_hcd_err[item][i] = pars[3]['bias_hcd']
        else:
            b_hcd[item] = np.bool_(b_hcd[item])
        cov[item][i] = pars[4]['cov[beta_LYA, bias_eta_LYA]']

    bias[item] = bias_eta[item] * 0.97 / beta[item]
    bias_err[item] = util.bias_err(bias_eta[item], bias_eta_err[item], beta[item], beta_err[item], cov[item])
    beff[item] = bias[item] * np.sqrt(1+2/3*beta[item]+1/5*beta[item]**2)
    beff_err[item] = util.beff_err(bias_eta[item], bias_eta_err[item], beta[item], beta_err[item], cov[item], pars[2]['growth_rate'])

    # Correct redshift dependency
    if correct_z_dep:
        print("Correcting the redshift dependency of b_eff...")
        bias_eta[item] /= (1+redshift[item])**gamma['bias_eta']
        beta[item] /= (1+redshift[item])**gamma['beta']
        bias[item] /= (1+redshift[item])**gamma['bias']
        beff[item] /= (1+redshift[item])**gamma['beff']
        bias_eta_err[item] /= (1+redshift[item])**gamma['bias_eta']
        beta_err[item] /= (1+redshift[item])**gamma['beta']
        bias_err[item] /= (1+redshift[item])**gamma['bias']
        beff_err[item] /= (1+redshift[item])**gamma['beff']

    # Take absolute value on bias
    bias_eta[item] = np.abs(bias_eta[item])
    bias[item] = np.abs(bias[item])
    beff[item] = np.abs(beff[item])
    b_hcd[item] = np.abs(b_hcd[item])

    if 'coadd' in item: continue

    if np.array_equal([0,0,0,0],bias_eta_err[item]): continue

    print("Fits on {}:".format(item))
    p_bias[item] = sp.optimize.curve_fit(func, redshift[item], bias[item], sigma=bias_err[item])
    print("bias(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_bias[item][0][0], p_bias[item][1][0,0], p_bias[item][0][1], p_bias[item][1][1,1]))
    p_beta[item] = sp.optimize.curve_fit(func, redshift[item], beta[item], sigma=beta_err[item])
    print("beta(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_beta[item][0][0], p_beta[item][1][0,0], p_beta[item][0][1], p_beta[item][1][1,1]))
    p_bias_eta[item] = sp.optimize.curve_fit(func, redshift[item], bias_eta[item], sigma=bias_eta_err[item])
    print("bias_eta(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_bias_eta[item][0][0], p_bias_eta[item][1][0,0], p_bias_eta[item][0][1], p_bias_eta[item][1][1,1]))
    p_beff[item] = sp.optimize.curve_fit(func, redshift[item], beff[item], sigma=beff_err[item])
    print("beff(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_beff[item][0][0], p_beff[item][1][0,0], p_beff[item][0][1], p_beff[item][1][1,1]))


### Plots
if plot_bias_eta:
    fig1, ax1 = plt.subplots()
    for i, item in enumerate(toplot):
        label = labels[i]
        fmt = 'o'
        if 'coadd' in item:
            label = None
            fmt = 'x'
        ax1.errorbar(redshift[item], bias_eta[item], yerr=bias_eta_err[item], fmt=fmt, label=label, color=colors[item])
        if 'coadd' in item: continue
        z = np.linspace(redshift[item].min(), redshift[item].max(), 100)
        if item in p_bias_eta.keys():
            plt.plot(z, func(z, p_bias_eta[item][0][0], p_bias_eta[item][0][1]), linestyle='--', color=colors[item])
            if plot_shades:
                plt.fill_between(z, func(z, p_bias_eta[item][0][0]+p_bias_eta[item][1][0,0], p_bias_eta[item][0][1]-p_bias_eta[item][1][1,1]),
                                 func(z, p_bias_eta[item][0][0]-p_bias_eta[item][1][0,0], p_bias_eta[item][0][1]+p_bias_eta[item][1][1,1]),
                                 color=colors[item], alpha=0.2)
    if plot_pred1:
        ax1.plot(zpred1, bias_eta_pred1, 'x', color='darkgreen', label='pred 10<r<180')
    if plot_pred2:
        ax1.plot(zpred2, bias_eta_pred2, 'x', color='limegreen', label='pred 40<r<180')
    ax1.legend()
    ax1.grid()
    ax1.set_xlabel('z')
    if absolute_bias:
        ylabel = r'$|\mathrm{bias\_eta}_{LYA}|$'
    else:
        ylabel = r'$\mathrm{bias\_eta}_{LYA}$'
    if correct_z_dep:
        ylabel += ' / (1+z)^{}'.format(gamma['bias_eta'])
    ax1.set_ylabel(ylabel)
    plt.tight_layout()

if plot_beta:
    fig2, ax2 = plt.subplots()
    for i, item in enumerate(toplot):
        label = labels[i]
        fmt = 'o'
        if 'coadd' in item:
            label = None
            fmt = 'x'
        ax2.errorbar(redshift[item], beta[item], yerr=beta_err[item], fmt=fmt, label=label, color=colors[item])
        if 'coadd' in item: continue
        z = np.linspace(redshift[item].min(), redshift[item].max(), 100)
        if item in p_beta.keys():
            plt.plot(z, func(z, p_beta[item][0][0], p_beta[item][0][1]), linestyle='--', color=colors[item])
            if plot_shades:
                plt.fill_between(z, func(z, p_beta[item][0][0]+p_beta[item][1][0,0], p_beta[item][0][1]-p_beta[item][1][1,1]),
                                 func(z, p_beta[item][0][0]-p_beta[item][1][0,0], p_beta[item][0][1]+p_beta[item][1][1,1]),
                                 color=colors[item], alpha=0.2)
    if plot_pred1:
        ax2.plot(zpred1, beta_pred1, 'x', color='darkgreen', label='pred 10<r<180')
    if plot_pred2:
        ax2.plot(zpred2, beta_pred2, 'x', color='limegreen', label='pred 40<r<180')
    ax2.legend()
    ax2.grid()
    ax2.set_xlabel('z')
    ylabel = r'$\beta_{LYA}$'
    if correct_z_dep:
        ylabel += ' / (1+z)^{}'.format(gamma['beta'])
    ax2.set_ylabel(ylabel)
    plt.tight_layout()

if plot_bias:
    fig3, ax3 = plt.subplots()
    for i, item in enumerate(toplot):
        label = labels[i]
        fmt = 'o'
        if 'coadd' in item:
            label = None
            fmt = 'x'
        ax3.errorbar(redshift[item], bias[item], yerr=bias_err[item], fmt=fmt, label=label, color=colors[item])
        if 'coadd' in item: continue
        z = np.linspace(redshift[item].min(), redshift[item].max(), 100)
        if item in p_bias.keys():
            plt.plot(z, func(z, p_bias[item][0][0], p_bias[item][0][1]), linestyle='--', color=colors[item])
            if plot_shades:
                plt.fill_between(z, func(z, p_bias[item][0][0]+p_bias[item][1][0,0], p_bias[item][0][1]-p_bias[item][1][1,1]),
                                 func(z, p_bias[item][0][0]-p_bias[item][1][0,0], p_bias[item][0][1]+p_bias[item][1][1,1]),
                                 color=colors[item], alpha=0.2)
    if plot_pred1:
        ax3.plot(zpred1, bias_pred1, 'x', color='darkgreen', label='pred 10<r<180')
    if plot_pred2:
        ax3.plot(zpred2, bias_pred2, 'x', color='limegreen', label='pred 40<r<180')
    ax3.legend()
    ax3.grid()
    ax3.set_xlabel('z')
    if absolute_bias:
        ylabel = r'$|\mathrm{b}_{LYA}|$'
    else:
        ylabel = r'$\mathrm{b}_{LYA}$'
    if correct_z_dep:
        ylabel += ' / (1+z)^{}'.format(gamma['bias'])
    ax3.set_ylabel(ylabel)
    plt.tight_layout()

if plot_beff:
    fig4, ax4 = plt.subplots()
    for i, item in enumerate(toplot):
        label = labels[i]
        fmt = 'o'
        if 'coadd' in item:
            label = None
            fmt = 'x'
        ax4.errorbar(redshift[item], beff[item], yerr=beff_err[item], fmt=fmt, label=label, color=colors[item])
        if 'coadd' in item: continue
        z = np.linspace(redshift[item].min(), redshift[item].max(), 100)
        if item in p_beff.keys():
            plt.plot(z, func(z, p_beff[item][0][0], p_beff[item][0][1]), linestyle='--', color=colors[item])
            if plot_shades:
                plt.fill_between(z, func(z, p_beff[item][0][0]+p_beff[item][1][0,0], p_beff[item][0][1]-p_beff[item][1][1,1]),
                                 func(z, p_beff[item][0][0]-p_beff[item][1][0,0], p_beff[item][0][1]+p_beff[item][1][1,1]),
                                 color=colors[item], alpha=0.2)
    if plot_pred1:
        ax4.plot(zpred1, beff_pred1, 'x', color='darkgreen', label='pred 10<r<180')
    if plot_pred2:
        ax4.plot(zpred2, beff_pred2, 'x', color='limegreen', label='pred 40<r<180')
    ax4.legend()
    ax4.grid()
    ax4.set_xlabel('z')
    if absolute_bias:
        ylabel=r'$|\mathrm{beff}_{LYA}|$'
    else:
        ylabel=r'$\mathrm{beff}_{LYA}$'
    if correct_z_dep:
        ylabel += ' / (1+z)^{}'.format(gamma['beff'])
    ax4.set_ylabel(ylabel)
    plt.tight_layout()

if plot_beff:
    fig5, ax5 = plt.subplots()
    for i, item in enumerate(toplot):
        # beff[item] *= (1+z0) / (1+redshift[item])
        label = labels[i]
        fmt = 'o'
        if 'coadd' in item:
            label = None
            fmt = 'x'
        ax5.errorbar(redshift[item], beff[item]*(1+z0)/(1+redshift[item]), yerr=beff_err[item]*(1+z0)/(1+redshift[item]), fmt=fmt, label=label, color=colors[item])
        if 'coadd' in item: continue
        z = np.linspace(redshift[item].min(), redshift[item].max(), 100)
        if item in p_beff.keys():
            plt.plot(z, func(z, p_beff[item][0][0], p_beff[item][0][1])*(1+z0)/(1+z), linestyle='--', color=colors[item])
            if plot_shades:
                plt.fill_between(z, func(z, p_beff[item][0][0]+p_beff[item][1][0,0], p_beff[item][0][1]-p_beff[item][1][1,1])*(1+z0)/(1+z),
                                 func(z, p_beff[item][0][0]-p_beff[item][1][0,0], p_beff[item][0][1]+p_beff[item][1][1,1])*(1+z0)/(1+z),
                                 color=colors[item], alpha=0.2)
    if plot_pred1:
        ax5.plot(zpred1, beff_pred1*(1+z0)/(1+zpred1), 'x', color='darkgreen', label='pred 10<r<180')
    if plot_pred2:
        ax5.plot(zpred2, beff_pred2*(1+z0)/(1+zpred2), 'x', color='limegreen', label='pred 40<r<180')
    ax5.legend()
    ax5.grid()
    ax5.set_xlabel('z')
    if absolute_bias:
        ylabel = r'$|\mathrm{beff}_{LYA}|$'+' * G(z) / G({})'.format(z0)
    else:
        ylabel = r'$\mathrm{beff}_{LYA}$'+' * G(z) / G({})'.format(z0)
    if correct_z_dep:
        ylabel += ' / (1+z)^{}'.format(gamma['beff'])
    ax5.set_ylabel(ylabel)
    plt.tight_layout()

if plot_b_hcd:
    fig6, ax6 = plt.subplots()
    for i, item in enumerate(toplot):
        if b_hcd[item][0] is None: continue
        label = labels[i]
        # fmt = 'o'
        # if 'coadd' in item:
        #     label = None
        #     fmt = 'x'
        if np.array_equal(b_hcd_err[item], np.zeros(len(b_hcd[item]))): continue
        ax6.errorbar(redshift[item], b_hcd[item], yerr=b_hcd_err[item], marker='o', label=label, color=colors[item])
    ax6.legend()
    ax6.grid()
    ax6.set_xlabel('z')
    if absolute_bias:
        ylabel='|b_hcd|'
    else:
        ylabel='b_hcd'
    ax6.set_ylabel(ylabel)
    plt.tight_layout()

plt.show()
