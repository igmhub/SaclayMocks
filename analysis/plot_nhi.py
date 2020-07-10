import numpy as np
import glob
from SaclayMocks import util
import matplotlib.pyplot as plt


# pyplot config
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

files = {}
colors = {}
nhi = {}
chi2 = {}
beff_lya = {}
beff_err_lya = {}
beta_lya = {}
beta_err_lya = {}
b_hcd = {}
b_err_hcd = {}
beta_hcd = {}
beta_err_hcd = {}
L0 = {}
L0_err = {}

### Parameters and config:
## References
# refs = ['dr16', 'mock-0.0', 'mock-0.2']
refs = ['mock-0.0']

# files['dr16'] = "/global/cscratch1/sd/tetourne/Analysis/dr16_no_dla_masking/Fits/cf/Kaiser_sky_met_hcd/z_0_10/result_fvoigt_v4.7.22.h5"
files['dr16'] = "/global/project/projectdirs/eboss/lya_forest/dr16/redo_4_zbins/Fits/cf/Kaiser_sky_met_hcd/z_0_10/result.h5"

# files['mock-0.0']  = "/global/cscratch1/sd/tetourne/Out/v4.7.22/from_quickquasars/Fit/result_cf.h5"
files['mock-0.0'] = "/global/project/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output_old/eboss-0.0_seed126429/cf_z_0_10-exp.h5"

# files['mock-0.2']  = "/global/cscratch1/sd/tetourne/Out/v4.7.22_dla/from_quickquasars/Fit/result_cf.h5"
files['mock-0.2']  = "/global/project/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output_old/eboss-0.2_mask_nhi20.3/cf_z_0_10-exp_Rogers.h5"

## version of mocks to plot
# mocks = ['v4.7.22', 'downsampling0.5', 'lowcut_nhi20.3', 'highcut_nhi20.3', 'nhi_bin']
# mocks = ['v4.7.22', 'lowcut_nhi20.3', 'highcut_nhi20.3', 'nhi_bin']
# mocks = ['v4.7.22', 'nhi_bin', 'nhi_bin_fixed_lya-0.0', 'nhi_bin_fixed_lya-0.2']
# mocks = ['v4.7.22']
# mocks = ['v4.7.22', 'nhi_bin']
# mocks = ['c-g', 'rogers_1', 'rogers_2']
# mocks = ['rogers_1', 'rogers_2', 'rogers_3']
mocks = ['rogers_3']
# mocks = ['rogers_1', 'rogers_2', 'rogers_fixed_lya']
# mocks = ['rogers', 'rogers_rmin10', 'c-g', 'c-g_rmin10']
# mocks = ['rogers', 'rogers_rmin10']
# mocks = ['c-g', 'c-g_rmin10']
# mocks = ['rogers_fixed_lya']

labels_mocks = mocks
# labels_mocks = ['C-G', 'Rogers 1', 'Rogers 2']
# labels_mocks = ['Rogers 1', 'Rogers 2', 'Rogers 3']
# labels_mocks = ['all', 'nhi>20.3', 'nhi<20.3', 'nhi_bin']
# labels_mocks = ['all', 'nhi_bin', 'nhi_bin_fixed_lya-0.0', 'nhi_bin_fixed_lya-0.2']
# labels_mocks = ['v4.7.22_nhi_bins']

# filename_mock = "/global/projecta/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output_nhi*/eboss-0.2/Fits/result_fvoigt_v4.7.22.h5"
# files['v4.7.22'] = "/global/projecta/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output_nhi*/eboss-0.2/Fits/result_fvoigt_v4.7.22.h5"
# files['downsampling0.5'] = "/global/projecta/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output_nhi*/eboss-0.2/Fits/result_fvoigt_v4.7.22_downsampling0.5.h5"
# files['lowcut_nhi20.3'] = "/global/projecta/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output_nhi*/eboss-0.2/Fits/result_fvoigt_v4.7.22_lowcut_nhi20.3.h5"
# files['highcut_nhi20.3'] = "/global/projecta/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output_nhi*/eboss-0.2/Fits/result_fvoigt_v4.7.22_highcut_nhi20.3.h5"
# files['nhi_bin'] = "/global/projecta/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output_nhi*/eboss-0.2/Fits/result_fvoigt_v4.7.22_nhi_??.?_??.?.h5"
# files['nhi_bin_fixed_lya-0.0'] = "/global/projecta/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output_nhi*/eboss-0.2/Fits/result_fvoigt_v4.7.22_nhi_*_fixed_lya-0.0.h5"
# files['nhi_bin_fixed_lya-0.2'] = "/global/projecta/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output_nhi*/eboss-0.2/Fits/result_fvoigt_v4.7.22_nhi_*_fixed_lya-0.2.h5"
files['rogers_1'] = "/global/project/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output_nhi_*/eboss-0.2/result_rogers.h5"
files['rogers_2'] = "/global/project/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output_nhi_*/eboss-0.2/result_rogers_method2.h5"
files['rogers_3'] = "/global/project/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output_nhi_*/eboss-0.2/result_rogers_fixed_L0.h5"
files['rogers_rmin10'] = "/global/project/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output_nhi_*/eboss-0.2/result_rogers_rmin10.h5"
files['rogers_fixed_lya'] = "/global/project/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output_nhi_*/eboss-0.2/result_rogers_fixed_lya.h5"
files['c-g'] = "/global/project/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output_nhi_*/eboss-0.2/result_c-g.h5"
files['c-g_rmin10'] = "/global/project/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output_nhi_*/eboss-0.2/result_c-g_rmin10.h5"

delta_nhi = np.array([0,0,0,0,0,0,0,0,0])

# L0
nhi_for_L0 = np.array([17.6, 18.2, 18.8, 19.4, 20])
L0['method_1'] = np.array([0.28, 0.60, 1.2, 2.4, 4.8])  # using C-G
L0['method_2'] = np.array([2, 3, 4, 8, 16])  # Using Jim's model


# colors
colors['dr16'] = 'black'
colors['mock-0.0'] = 'royalblue'
colors['mock-0.2'] = 'darkgreen'
colors['rogers_1'] = 'orange'
colors['rogers_2'] = 'darkviolet'
colors['rogers_3'] = 'deeppink'
colors['rogers_rmin10'] = 'darkorange'
colors['rogers_fixed_lya'] = 'deeppink'
colors['c-g'] = 'red'
colors['c-g_rmin10'] = 'darkorange'
# colors['v4.7.22'] = 'red'
# colors['downsampling0.5'] = 'darkorange'
# colors['lowcut_nhi20.3'] = 'darkblue'
# colors['highcut_nhi20.3'] = 'darkmagenta'
# colors['nhi_bin'] = 'darkorange'
# colors['nhi_bin_fixed_lya-0.0'] = 'indigo'
# colors['nhi_bin_fixed_lya-0.2'] = 'deeppink'

# title = '0 < z < 10 - Fvoigt_v4.7.22.txt'
title = ''
absolute_bias = True
legends = False

# Read refs
for item in refs:
    print("Reading from {}".format(files[item]))
    pars = util.extract_h5file(files[item])
    chi2[item] = pars[2]['chi2']
    beff_lya[item] = pars[2]['beff_LYA']
    beff_err_lya[item] = pars[3]['beff_LYA']
    beta_lya[item] = pars[2]['beta_LYA']
    beta_err_lya[item] = pars[3]['beta_LYA']
    if 'bias_hcd' in pars[2].keys():
        b_hcd[item] = pars[2]['bias_hcd']
        b_err_hcd[item] = pars[3]['bias_hcd']
        beta_hcd[item] = pars[2]['beta_hcd']
        beta_err_hcd[item] = pars[3]['beta_hcd']
    else:
        b_hcd[item] = None
        b_err_hcd[item] = None
        beta_hcd[item] = None
        beta_err_hcd[item] = None
    if 'L0_hcd' in pars[2].keys():
        L0[item] = pars[2]['L0_hcd']
        L0_err[item] = pars[3]['L0_hcd']
    else:
        L0[item] = None
        L0_err[item] = None

# Read bins
for i, item in enumerate(mocks):
    files_mock = glob.glob(files[item])
    print("Reading mocks from {}".format(files_mock))
    nhi[item] = []
    chi2[item] = []
    beff_lya[item] = []
    beff_err_lya[item] = []
    beta_lya[item] = []
    beta_err_lya[item] = []
    b_hcd[item] = []
    b_err_hcd[item] = []
    beta_hcd[item] = []
    beta_err_hcd[item] = []
    L0[item] = []
    L0_err[item] = []
    
    for f in files_mock:
        idx = f.find('output_nhi') + 11
        nhi[item].append(np.mean([float(f[idx:idx+4]), float(f[idx+5:idx+9])]))
        pars = util.extract_h5file(f)
        chi2[item].append(pars[2]['chi2'])
        beff_lya[item].append(pars[2]['beff_LYA'])
        beff_err_lya[item].append(pars[3]['beff_LYA'])
        beta_lya[item].append(pars[2]['beta_LYA'])
        beta_err_lya[item].append(pars[3]['beta_LYA'])
        b_hcd[item].append(pars[2]['bias_hcd'])
        b_err_hcd[item].append(pars[3]['bias_hcd'])
        beta_hcd[item].append(pars[2]['beta_hcd'])
        beta_err_hcd[item].append(pars[3]['beta_hcd'])
        if 'L0_hcd' in pars[2].keys():
            L0[item].append(pars[2]['L0_hcd'])
            L0_err[item].append(pars[3]['L0_hcd'])

    nhi[item] += delta_nhi[i]

# Take absolute values of biases
if absolute_bias:
    for item in mocks+refs:
        beff_lya[item] = np.abs(beff_lya[item])
        if b_hcd[item] is not None:
            b_hcd[item] = np.abs(b_hcd[item])

### Plots
print("Plotting...")
# beff_lya
f1, ax1 = plt.subplots()
nhi_min = np.min([nhi[item] for item in mocks])-0.5
nhi_max = np.max([nhi[item] for item in mocks])+0.5
for i, item in enumerate(mocks):
    label = labels_mocks[i]
    ax1.errorbar(nhi[item], beff_lya[item], yerr=beff_err_lya[item], fmt='o', label=label, color=colors[item])
for item in refs:
    ax1.axhline(beff_lya[item], label=item, color=colors[item])
    ax1.axhline(beff_lya[item] - beff_err_lya[item], linestyle='--', color=colors[item])
    ax1.axhline(beff_lya[item] + beff_err_lya[item], linestyle='--', color=colors[item])
    ax1.fill_between([nhi_min, nhi_max], [beff_lya[item] - beff_err_lya[item], beff_lya[item] - beff_err_lya[item]],
                     y2=[beff_lya[item] + beff_err_lya[item], beff_lya[item] + beff_err_lya[item]], color=colors[item], alpha=0.2)
ax1.set_xlim(nhi_min, nhi_max)
ax1.set_xlabel(r'$\log n_{\mathrm{HI}}$')
if absolute_bias:
    ax1.set_ylabel(r'$|b_{\mathrm{eff},\mathrm{Ly}\alpha}|$')
else:
    ax1.set_ylabel(r'$b_{\mathrm{eff},\mathrm{Ly}\alpha}$')
ax1.set_title(title)
if legends:
    ax1.legend()
ax1.grid()
plt.tight_layout()

# beta_lya
f1, ax1 = plt.subplots()
nhi_min = np.min([nhi[item] for item in mocks])-0.5
nhi_max = np.max([nhi[item] for item in mocks])+0.5
for i, item in enumerate(mocks):
    label = labels_mocks[i]
    ax1.errorbar(nhi[item], beta_lya[item], yerr=beta_err_lya[item], fmt='o', label=label, color=colors[item])
for item in refs:
    ax1.axhline(beta_lya[item], label=item, color=colors[item])
    ax1.axhline(beta_lya[item] - beta_err_lya[item], linestyle='--', color=colors[item])
    ax1.axhline(beta_lya[item] + beta_err_lya[item], linestyle='--', color=colors[item])
    ax1.fill_between([nhi_min, nhi_max], [beta_lya[item] - beta_err_lya[item], beta_lya[item] - beta_err_lya[item]],
                     y2=[beta_lya[item] + beta_err_lya[item], beta_lya[item] + beta_err_lya[item]], color=colors[item], alpha=0.2)
ax1.set_xlim(nhi_min, nhi_max)
ax1.set_xlabel(r'$\log n_{\mathrm{HI}}$')
ax1.set_ylabel(r'$\beta_{\mathrm{Ly}\alpha}}$')
ax1.set_title(title)
if legends:
    ax1.legend()

ax1.grid()
plt.tight_layout()

# b_hcd
f1, ax1 = plt.subplots()
nhi_min = np.min([nhi[item] for item in mocks])-0.5
nhi_max = np.max([nhi[item] for item in mocks])+0.5
for i, item in enumerate(mocks):
    label = labels_mocks[i]
    ax1.errorbar(nhi[item], b_hcd[item], yerr=b_err_hcd[item], fmt='o', label=label, color=colors[item])
for item in refs:
    if b_hcd[item] is not None:
        ax1.axhline(b_hcd[item], label=item, color=colors[item])
        ax1.axhline(b_hcd[item] - b_err_hcd[item], linestyle='--', color=colors[item])
        ax1.axhline(b_hcd[item] + b_err_hcd[item], linestyle='--', color=colors[item])
        ax1.fill_between([nhi_min, nhi_max], [b_hcd[item] - b_err_hcd[item], b_hcd[item] - b_err_hcd[item]],
                         y2=[b_hcd[item] + b_err_hcd[item], b_hcd[item] + b_err_hcd[item]], color=colors[item], alpha=0.2)
ax1.set_xlim(nhi_min, nhi_max)
ax1.set_xlabel(r'$\log n_{\mathrm{HI}}$')
if absolute_bias:
    ax1.set_ylabel(r'$|b_{\mathrm{HCD}}|$')
else:
    ax1.set_ylabel(r'$b_{\mathrm{HCD}}$')
ax1.set_title(title)
if legends:
    ax1.legend()

ax1.grid()
plt.tight_layout()

# beta_hcd
f1, ax1 = plt.subplots()
nhi_min = np.min([nhi[item] for item in mocks])-0.5
nhi_max = np.max([nhi[item] for item in mocks])+0.5
for i, item in enumerate(mocks):
    label = labels_mocks[i]
    ax1.errorbar(nhi[item], beta_hcd[item], yerr=beta_err_hcd[item], fmt='o', label=label, color=colors[item])
for item in refs:
    if beta_hcd[item] is not None:
        ax1.axhline(beta_hcd[item], label=item, color=colors[item])
        ax1.axhline(beta_hcd[item] - beta_err_hcd[item], linestyle='--', color=colors[item])
        ax1.axhline(beta_hcd[item] + beta_err_hcd[item], linestyle='--', color=colors[item])
        ax1.fill_between([nhi_min, nhi_max], [beta_hcd[item] - beta_err_hcd[item], beta_hcd[item] - beta_err_hcd[item]],
                         y2=[beta_hcd[item] + beta_err_hcd[item], beta_hcd[item] + beta_err_hcd[item]], color=colors[item], alpha=0.2)
ax1.set_xlim(nhi_min, nhi_max)
ax1.set_xlabel(r'$\log n_{\mathrm{HI}}$')
ax1.set_ylabel(r'$\beta_{\mathrm{HCD}}$')
ax1.set_title(title)
if legends:
    ax1.legend()

ax1.grid()
plt.tight_layout()

# chi2
f1, ax1 = plt.subplots()
nhi_min = np.min([nhi[item] for item in mocks])-0.5
nhi_max = np.max([nhi[item] for item in mocks])+0.5
for i, item in enumerate(mocks):
    label = labels_mocks[i]
    ax1.plot(nhi[item], chi2[item], 'o', label=label, color=colors[item])
for item in refs:
    ax1.axhline(chi2[item], label=item, color=colors[item])
ax1.set_xlim(nhi_min, nhi_max)
ax1.set_xlabel(r'$\log n_{\mathrm{HI}}$')
ax1.set_ylabel(r'$\chi^2$')
ax1.set_title(title)
if legends:
    ax1.legend()

ax1.grid()
plt.tight_layout()

# L0_hcd
f1, ax1 = plt.subplots()
nhi_min = np.min([nhi[item] for item in mocks])-0.5
nhi_max = np.max([nhi[item] for item in mocks])+0.5
for i, item in enumerate(mocks):
    label = labels_mocks[i]
    ax1.errorbar(nhi[item], L0[item], yerr=L0_err[item], fmt='o', label=label, color=colors[item])
# for item in refs:
#     if L0[item] is not None:
#         ax1.axhline(L0[item], label=item, color=colors[item])
#         ax1.axhline(L0[item] - L0_err[item], linestyle='--', color=colors[item])
#         ax1.axhline(L0[item] + L0_err[item], linestyle='--', color=colors[item])
#         ax1.fill_between([nhi_min, nhi_max], [L0[item] - L0_err[item], L0[item] - L0_err[item]],
#                          y2=[L0[item] + L0_err[item], L0[item] + L0_err[item]], color=colors[item], alpha=0.2)
ax1.axhline(10, linestyle='--', color='black')
# ax1.plot(nhi_for_L0, L0['method_1'], 'o', color='royalblue', label='method 1')
# ax1.plot(nhi_for_L0, L0['method_2'], 'o', color='darkblue', label='method 2')
ax1.set_xlim(nhi_min, nhi_max)
ax1.set_xlabel(r'$\log n_{\mathrm{HI}}$')
ax1.set_ylabel(r'$L_{\mathrm{HCD}}$')
ax1.set_title(title)
if legends:
    ax1.legend()

ax1.grid()
plt.tight_layout()


plt.show()
