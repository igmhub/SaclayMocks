import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from SaclayMocks import util
from iminuit import Minuit


SMALL_SIZE = 11
MEDIUM_SIZE = 14
BIGGER_SIZE = 14
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize

plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.rc('figure', figsize=(11,8))


toplot = [
    "/global/project/projectdirs/eboss/lya_forest/dr16/redo_4_zbins/Fits/cf/Kaiser_sky_met_hcd/z_{zmin}_{zmax}/result.h5"
    # "/global/project/projectdirs/eboss/lya_forest/dr16/redo_4_zbins/Fits/cf/Kaiser_sky_met_hcd/z_{zmin}_{zmax}/result_Rogers2.8_rmin10.h5"
    # ,"/global/project/projectdirs/eboss/lya_forest/dr16/redo_4_zbins/Fits/cf/Kaiser_sky_met_hcd/z_{zmin}_{zmax}/result_Rogers2.8_rmin20.h5"
    ,"/global/project/projectdirs/eboss/lya_forest/dr16/redo_4_zbins/Fits/cf/Kaiser_sky_met_hcd/z_{zmin}_{zmax}/result_C-G.h5"
    # ,"/global/project/projectdirs/eboss/lya_forest/dr16/redo_4_zbins/Fits/cf/Kaiser_sky_met_hcd/z_{zmin}_{zmax}/result_C-G_rmin20.h5"
    ,"/global/cscratch1/sd/tetourne/Analysis/dr16_no_dla_masking/Fits/cf/Kaiser_sky_met_hcd/z_{zmin}_{zmax}/result.h5"

    # "/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-raw/cf_z_{zmin}_{zmax}-exp.h5"

    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.0/cf_z_{zmin}_{zmax}-exp.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.0/cf_z_{zmin}_{zmax}-exp_distorsion_corrected.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.0/xcf_z_{zmin}_{zmax}-exp.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.0/xcf_z_{zmin}_{zmax}-exp_fixed_drp.h5"
    #,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.0/xcf_z_{zmin}_{zmax}-exp_rmin40.h5"

    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.0/cf_z_{zmin}_{zmax}-exp_rmin10.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.0/xcf_z_{zmin}_{zmax}-exp_rmin30.h5"

    # "/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.0/xcf_z_{zmin}_{zmax}-exp_fixed_drp.h5"

    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2/cf_z_{zmin}_{zmax}-exp.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2/xcf_z_{zmin}_{zmax}-exp.h5"

    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2/cf_z_{zmin}_{zmax}-exp_Rogers2.8.h5"
    # "/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2/cf_z_{zmin}_{zmax}-exp_C-G.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2/xcf_z_{zmin}_{zmax}-exp_C-G.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2/xcf_z_{zmin}_{zmax}-exp_C-G_rmin10.h5"
    #,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2/xcf_z_{zmin}_{zmax}-exp_C-G_rmin40.h5"

    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2/cf_z_{zmin}_{zmax}-exp_rmin10.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2/xcf_z_{zmin}_{zmax}-exp_rmin10.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2/xcf_z_{zmin}_{zmax}-exp_fixed_lya.h5"

    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_10/eboss-0.2_no_hcd_masking/cf_z_{zmin}_{zmax}-exp.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_10/eboss-0.2_no_hcd_masking/xcf_z_{zmin}_{zmax}-exp.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_10/eboss-0.2_no_hcd_masking/xcf_z_{zmin}_{zmax}-exp_rmin10.h5"

    #,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2_no_hcd_masking/cf_z_{zmin}_{zmax}-exp.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2_no_hcd_masking/xcf_z_{zmin}_{zmax}-exp.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2_no_hcd_masking/xcf_z_{zmin}_{zmax}-exp_fixed_drp.h5"
    #,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2_no_hcd_masking/xcf_z_{zmin}_{zmax}-exp_rmin40.h5"
    ]

refs = [
    # "/global/project/projectdirs/eboss/lya_forest/dr16/redo_4_zbins/Fits/cf/Kaiser_sky_met_hcd/z_0_10/result.h5"
    # "/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-raw/cf_z_0_10-exp.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-raw/cf_z_0_2.65-exp.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-raw/cf_z_2.35_3.05-exp.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-raw/cf_z_2.65_10-exp.h5"
    # "/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.0/cf_z_0_10-exp.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.0/cf_z_0_2.65-exp.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.0/cf_z_2.35_3.05-exp.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.0/cf_z_2.65_10-exp.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.0/xcf_z_0_10-exp.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.0/xcf_z_0_10-exp_fixed_drp.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.0/xcf_z_0_10-exp_rmin40.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2/cf_z_0_10-exp_Rogers2.8.h5"
    # "/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2/cf_z_0_10-exp_C-G.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2/cf_z_0_2.65-exp_C-G.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2/cf_z_2.35_3.05-exp_C-G.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2/cf_z_2.65_10-exp_C-G.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2/xcf_z_0_10-exp_C-G.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2/xcf_z_0_10-exp_C-G_rmin40.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_10/eboss-0.2_no_hcd_masking/cf_z_0_10-exp.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2_no_hcd_masking/cf_z_0_10-exp.h5"
    #,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2_no_hcd_masking/xcf_z_0_10-exp.h5"
    #,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2_no_hcd_masking/xcf_z_0_10-exp_fixed_drp.h5"
    # ,"/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2_no_hcd_masking/xcf_z_0_10-exp_rmin40.h5"
]

# labels = ['DR16', 'auto-raw', 'auto-0.0', 'cross-0.0', 'auto-0.2', 'cross-0.2']
# labels = ['DR16', 'auto-0.0_rmin20', 'cross-0.0_rmin20','auto-0.0_rmin10', 'cross-0.0_rmin10']
# labels = ['DR16', 'auto-0.0', 'cross-0.0', 'cross-0.2_rmin30']
# labels = ['DR16_Rogers', 'DR16_C-G', 'eboss-0.2_Rogers', 'eboss-0.2_C-G']
# labels = ['DR16_Rogers', 'DR16_C-G', 'eboss-0.0', 'eboss-0.2_Rogers', 'eboss-0.2_C-G']
# labels = ['Rogers_rmin10', 'Rogers_rmin20', 'CG_rmin_10(nhi<20.3)', 'CG_rmin20(nhi<20.3)', 'CG_rmin10']
# labels = ['eboss-0.0', 'Rogers', 'C-G(nhi<20.3)', 'C-G']
# labels = ['DR16', 'raw mocks', 'eboss-0.0', 'eboss-0.2']
# labels = ['eboss-raw', 'eboss-0-0', 'eboss-0.2 (mask)', 'eboss-0.2 (no mask)']
# labels = ['eboss-0-0', 'eboss-0.2 (mask)', 'eboss-0.2 (no mask)']
# labels = ['cf_rmin20', 'xcf_rmin20', 'xcf_rmin40']
# labels = ['raw mocks', 'cf-0.0', 'xcf-0.0', 'cf-0.2(mask)', 'xcf-0.2(mask)', 'cf-0.2', 'xcf-0.2']
# labels = ['Rogers', 'CG(nhi<20.3)', 'CG']
# labels = ['DR16', 'auto-0.2', 'cross-0.2', 'auto-0.2_C-G_nhi<20.3', 'cross-0.2_C-G_nhi<20.3', 'auto-0.2_C-G', 'cross-0.2_C-G']
# labels = ['DR16_Rogers', 'DR16_C-G', 'eboss-0.0', 'eboss-0.2_Rogers', 'eboss-0.2_C-G', 'eboss-0.2_C-G']
# labels = ['DR16', 'cross-0.2', 'cross-0.2_C-G_nhi<20.3', 'cross-0.2_C-G_nhi<20.3_rmin10', 'cross-0.2_C-G', 'cross-0.2_C-G_rmin10']
# labels = ['DR16', 'cross-0.2_Rogers', 'cross-0.2_C-G_nhi<20.3', 'cross-0.2_C-G']
labels = ['eboss-raw', 'eboss-0.0', 'eboss-0.0_corrected']

colors = ['b', 'green', 'orange', 'r']
# colors = ['green', 'orange', 'r']
# colors = ['b', 'green', 'green', 'orange', 'orange', 'r', 'r']
# colors = ['b', 'green', 'red']
# colors = ['b', 'green', 'grey', 'orange', 'r']
# colors = ['b', 'green', 'orange', 'red', 'deeppink']
# colors = ['black', 'b', 'green', 'r', 'r',  'g', 'r']
# fmts = ['o', 'o', 'o', 'o', 'o', 'o', 'o', 'o']
fmts = ['o', 'o', 'x', 'o', 'x', 'o', 'x']
# colors_ref = colors
colors_ref = ['r', 'green', 'green', 'green']
fmts_ref = ['x', 'x', 'x', 'x', 'x', 'x', 'x']

z_bins = np.array(['0', '2.35', '2.65', '3.05', '10'])

mu_bins = np.array([0, 0.5, 0.8, 0.95, 1])
# mu_bins = np.array([0,1])
# mu_bins = np.array([0, 0.2, 0.5, 1])
# colors = ['b', 'g', 'orange', 'red']
# colors = ['b', 'g', 'red']

plot_4_panels = True
legends = True
plot_refs = False
correct_z_dep = False
xmin=2.05
xmax=3.1
z_plot = np.linspace(xmin-0.05,xmax+0.05,1000)
delta_z = np.zeros(10)
# delta_z = np.array([-0.01, 0, 0.01])
# delta_z = np.array([-0.02, -0.01, 0.0, 0.01, 0.02, 0, 0, 0, 0, 0, 0, 0])  # small shift to better see the errorbars when the z_eff is the same
gamma_beff = 3.55
gamma_beta = -2.5
zeff = 2.4
### End of config


if plot_4_panels:
    f7, ax7 = plt.subplots(nrows=2, ncols=2)
else:
    f1, ax1 = plt.subplots()
    f2, ax2 = plt.subplots()
    f3, ax3 = plt.subplots()
    f4, ax4 = plt.subplots()
    f5, ax5 = plt.subplots()
    f6, ax6 = plt.subplots()

def func(x, a, gamma):
    return a*(1+x)**gamma

class Fitter(object):
    def __init__(self, x, y, y_err):
        self.x = x
        self.y = y
        self.y_err = y_err
    def chi2(self, a, gamma):
        chi2 = ((self.y - func(self.x, a, gamma)) / self.y_err)**2
        chi2 = chi2.sum()
        return chi2

for i, item in enumerate(toplot):
    redshift = np.zeros(len(z_bins)-1)
    beff = np.zeros(len(z_bins)-1)
    beff_err = np.zeros(len(z_bins)-1)
    bias = np.zeros(len(z_bins)-1)
    bias_err = np.zeros(len(z_bins)-1)
    beta = np.zeros(len(z_bins)-1)
    beta_err = np.zeros(len(z_bins)-1)
    b_hcd = np.zeros(len(z_bins)-1)
    b_hcd_err = np.zeros(len(z_bins)-1)
    beta_hcd = np.zeros(len(z_bins)-1)
    beta_hcd_err = np.zeros(len(z_bins)-1)
    growth_rate = np.zeros(len(z_bins)-1)
    for j in range(len(z_bins)-1):
        print("Reading {}".format(item.format(zmin=z_bins[j], zmax=z_bins[j+1])))
        pars = util.extract_h5file(item.format(zmin=z_bins[j], zmax=z_bins[j+1]))
        redshift[j] = pars[2]['zeff']
        beff[j] = pars[2]['beff_LYA']
        beff_err[j] = pars[3]['beff_LYA']
        bias[j] = pars[2]['bias_LYA']
        bias_err[j] = pars[3]['bias_LYA']
        beta[j] = pars[2]['beta_LYA']
        beta_err[j] = pars[3]['beta_LYA']
        growth_rate[j] = pars[2]['growth_rate']
        if 'bias_hcd' in pars[2]:
            b_hcd[j] = pars[2]['bias_hcd']
            b_hcd_err[j] = pars[3]['bias_hcd']
            beta_hcd[j] = pars[2]['beta_hcd']
            beta_hcd_err[j] = pars[3]['beta_hcd']

    # Take absolute value on bias
    beff = np.abs(beff)
    b_hcd = np.abs(b_hcd)
    bias = np.abs(bias)
    bias_sum = bias + b_hcd
    bias_sum_err = np.sqrt(bias_err**2 + b_hcd_err**2)
    bias_prod = bias*beta + b_hcd*beta_hcd
    bias_prod_err = np.sqrt(beta**2*bias_err**2 + bias**2*beta_err**2 + b_hcd**2*beta_hcd_err**2 + beta_hcd**2*b_hcd_err**2)
    if correct_z_dep:
        beff *= ((1+zeff)/(1+redshift))**gamma_beff
        beff_err *= ((1+zeff)/(1+redshift))**gamma_beff
        beta *= ((1+zeff)/(1+redshift))**gamma_beta
        beta_err *= ((1+zeff)/(1+redshift))**gamma_beta

    # Fitting with iminuit
    f_beta = Fitter(redshift, beta, beta_err)
    m_beta = Minuit(f_beta.chi2, a=28, error_a=5, gamma=-2, error_gamma=0.5, errordef=1)
    m_beta.migrad()
    m_beta.hesse()
    f_beff = Fitter(redshift, beff, beff_err)
    m_beff = Minuit(f_beff.chi2, a=0, error_a=1e-3, gamma=3.4, error_gamma=0.5, errordef=1)
    m_beff.migrad()
    m_beff.hesse()
    f_bias = Fitter(redshift, bias, bias_err)
    m_bias = Minuit(f_bias.chi2, a=0, error_a=1e-3, gamma=3.4, error_gamma=0.5, errordef=1)
    m_bias.migrad()
    m_bias.hesse()
    f_bias_sum = Fitter(redshift, bias_sum, bias_err)
    m_bias_sum = Minuit(f_bias_sum.chi2, a=0, error_a=1e-3, gamma=3.4, error_gamma=0.5, errordef=1)
    m_bias_sum.migrad()
    m_bias_sum.hesse()
    print("beta(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(m_beta.values[0], m_beta.errors[0], m_beta.values[1], m_beta.errors[1]))
    print("beff(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(m_beff.values[0], m_beff.errors[0], m_beff.values[1], m_beff.errors[1]))
    print("bias(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(m_bias.values[0], m_bias.errors[0], m_bias.values[1], m_bias.errors[1]))
    print("bias_sum(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(m_bias_sum.values[0], m_bias_sum.errors[0], m_bias_sum.values[1], m_bias_sum.errors[1]))

    # shift redshift to distinguish data sets in plots
    redshift += delta_z[i]

    # Plots
    if plot_4_panels:
        ax7[0,0].errorbar(redshift, beff, yerr=beff_err, fmt=fmts[i], color=colors[i], label=labels[i])
        ax7[0,0].plot(z_plot, func(z_plot, m_beff.values[0], m_beff.values[1]), linestyle='--', color=colors[i])
        ax7[0,1].errorbar(redshift, beta, yerr=beta_err, fmt=fmts[i], color=colors[i], label=labels[i])
        ax7[0,1].plot(z_plot, func(z_plot, m_beta.values[0], m_beta.values[1]), linestyle='--', color=colors[i])
        if 0 not in b_hcd:
            ax7[1,0].errorbar(redshift, b_hcd, yerr=b_hcd_err, fmt=fmts[i], color=colors[i], label=labels[i])
            # ax7[1,1].errorbar(redshift, bias_sum, yerr=bias_sum_err, fmt=fmts[i], color=colors[i], label=labels[i])
            ax7[1,1].errorbar(redshift, beta_hcd, yerr=beta_hcd_err, fmt=fmts[i], color=colors[i], label=labels[i])
    else:
        ax1.errorbar(redshift, beff, yerr=beff_err, fmt=fmts[i], color=colors[i], label=labels[i])
        ax1.plot(z_plot, func(z_plot, m_beff.values[0], m_beff.values[1]), linestyle='--', color=colors[i])
        ax2.errorbar(redshift, beta, yerr=beta_err, fmt=fmts[i], color=colors[i], label=labels[i])
        ax2.plot(z_plot, func(z_plot, m_beta.values[0], m_beta.values[1]), linestyle='--', color=colors[i])
        ax3.errorbar(redshift, bias, yerr=bias_err, fmt=fmts[i], color=colors[i], label=labels[i])
        ax3.plot(z_plot, func(z_plot, m_bias.values[0], m_bias.values[1]), linestyle='--', color=colors[i])
        if 0 not in b_hcd:
            ax4.errorbar(redshift, b_hcd, yerr=b_hcd_err, fmt=fmts[i], color=colors[i], label=labels[i])
            ax5.errorbar(redshift, bias_sum, yerr=bias_sum_err, fmt=fmts[i], color=colors[i], label=labels[i])
            # ax5.plot(z_plot, func(z_plot, m_bias_sum.values[0], m_bias_sum.values[1]), linestyle='--', color=colors[i])
            ax6.errorbar(redshift, bias_prod, yerr=bias_prod_err, fmt=fmts[i], color=colors[i], label=labels[i])

### Plots refs
if plot_refs:
    for i, item in enumerate(refs):
        print("Reading {}".format(item))
        pars = util.extract_h5file(item)
        redshift = pars[2]['zeff']
        beff = pars[2]['beff_LYA']
        beff_err = pars[3]['beff_LYA']
        bias = pars[2]['bias_LYA']
        bias_err = pars[3]['bias_LYA']
        beta = pars[2]['beta_LYA']
        beta_err = pars[3]['beta_LYA']
        growth_rate = pars[2]['growth_rate']
        if 'bias_hcd' in pars[2]:
            b_hcd = pars[2]['bias_hcd']
            b_hcd_err = pars[3]['bias_hcd']
            beta_hcd = pars[2]['beta_hcd']
            beta_hcd_err = pars[3]['beta_hcd']
        # Absolute values
        beff = np.abs(beff)
        b_hcd = np.abs(b_hcd)
        bias = np.abs(bias)
        bias_sum = bias + b_hcd
        bias_sum_err = np.sqrt(bias_err**2 + b_hcd_err**2)
        bias_prod = bias*beta + b_hcd*beta_hcd
        bias_prod_err = np.sqrt(beta**2*bias_err**2 + bias**2*beta_err**2 + b_hcd**2*beta_hcd_err**2 + beta_hcd**2*b_hcd_err**2)
        if correct_z_dep:
            beff *= ((1+zeff)/(1+redshift))**gamma_beff
            beff_err *= ((1+zeff)/(1+redshift))**gamma_beff
            beta *= ((1+zeff)/(1+redshift))**gamma_beta
            beta_err *= ((1+zeff)/(1+redshift))**gamma_beta
        # shift redshift to distinguish data sets in plots
        redshift += delta_z[i]
        # Plots
        if plot_4_panels:
            ax7[0,0].errorbar(redshift, beff, yerr=beff_err, fmt=fmts_ref[i], color=colors_ref[i])
            ax7[0,1].errorbar(redshift, beta, yerr=beta_err, fmt=fmts_ref[i], color=colors_ref[i])
            if 'bias_hcd' in pars[2]:
                ax7[1,0].errorbar(redshift, b_hcd, yerr=b_hcd_err, fmt=fmts_ref[i], color=colors_ref[i])
                ax7[1,1].errorbar(redshift, beta_hcd, yerr=beta_hcd_err, fmt=fmts_ref[i], color=colors_ref[i])
        else:
            ax1.errorbar(redshift, beff, yerr=beff_err, fmt=fmts_ref[i], color=colors_ref[i])
            ax2.errorbar(redshift, beta, yerr=beta_err, fmt=fmts_ref[i], color=colors_ref[i])
            ax3.errorbar(redshift, bias, yerr=bias_err, fmt=fmts_ref[i], color=colors_ref[i])
            if 'bias_hcd' in pars[2]:
                ax4.errorbar(redshift, b_hcd, yerr=b_hcd_err, fmt=fmts_ref[i], color=colors_ref[i])
                ax5.errorbar(redshift, bias_sum, yerr=bias_sum_err, fmt=fmts_ref[i], color=colors_ref[i])
                ax6.errorbar(redshift, bias_prod, yerr=bias_prod_err, fmt=fmts_ref[i], color=colors_ref[i])


if plot_4_panels:
    for ax in ax7.ravel():
        ax.grid()
        ax.set_xlabel('z')
        ax.set_xlim(xmin, xmax)
    if legends:
        ax7[0,0].legend()
    ax7[1,0].set_xlabel('z')
    ax7[1,1].set_xlabel('z')
    ax7[0,0].set_ylabel(r'$|b_{\mathrm{eff},\mathrm{Ly}\alpha}|$')
    ax7[0,1].set_ylabel(r'$\beta_{\mathrm{Ly}\alpha}$')
    ax7[1,0].set_ylabel(r'$|b_{\mathrm{HCD}}|$')
    # ax7[1,1].set_ylabel(r'$|b_{\mathrm{Ly}\alpha} + b_{\mathrm{HCD}}|$')
    ax7[1,1].set_ylabel(r'$\beta_{\mathrm{HCD}}$')
    f7.tight_layout()
else:
    for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
        ax.grid()
        if legends:
            ax.legend()
        ax.set_xlabel('z')
        ax.set_xlim(2.05, 3.5)

    ax1.set_ylabel(r'$|b_{\mathrm{eff},\mathrm{Ly}\alpha}|$')
    ax2.set_ylabel(r'$\beta_{\mathrm{Ly}\alpha}$')
    ax3.set_ylabel(r'$|b_{\mathrm{Ly}\alpha}|$')
    ax4.set_ylabel(r'$|b_{\mathrm{HCD}}|$')
    ax5.set_ylabel(r'$|b_{\mathrm{Ly}\alpha} + b_{\mathrm{HCD}}|$')
    ax6.set_ylabel(r'$|b_{\mathrm{Ly}\alpha}|\beta_{\mathrm{Ly}\alpha} + |b_{\mathrm{HCD}}|\beta_{\mathrm{HCD}}$')

    for f in [f1, f2, f3, f4, f5, f6]:
        f.tight_layout()

plt.show()
