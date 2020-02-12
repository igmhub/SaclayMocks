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
beff_lya = {}
beff_err_lya = {}
beta_lya = {}
beta_err_lya = {}
b_hcd = {}
b_err_hcd = {}
beta_hcd = {}
beta_err_hcd = {}

### Parameters and config:
## References
refs = ['dr16', 'mock-0.0', 'mock-0.2']

files['dr16'] = "/global/cscratch1/sd/tetourne/Analysis/dr16_no_dla_masking/Fits/cf/Kaiser_sky_met_hcd/z_0_10/result_fvoigt_v4.7.22.h5"

files['mock-0.0']  = "/global/cscratch1/sd/tetourne/Out/v4.7.22/from_quickquasars/Fit/result_cf.h5"

files['mock-0.2']  = "/global/cscratch1/sd/tetourne/Out/v4.7.22_dla/from_quickquasars/Fit/result_cf.h5"

## mocks
# mocks = ['v4.7.22', 'downsampling0.5', 'lowcut_nhi20.3', 'highcut_nhi20.3', 'nhi_bin']
mocks = ['v4.7.22', 'lowcut_nhi20.3', 'highcut_nhi20.3', 'nhi_bin']
# mocks = ['v4.7.22', 'nhi_bin']
labels_mocks = mocks
labels_mocks = ['all', 'nhi>20.3', 'nhi<20.3', 'nhi_bin']

# filename_mock = "/global/projecta/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output_nhi*/eboss-0.2/Fits/result_fvoigt_v4.7.22.h5"
files['v4.7.22'] = "/global/projecta/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output_nhi*/eboss-0.2/Fits/result_fvoigt_v4.7.22.h5"
files['downsampling0.5'] = "/global/projecta/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output_nhi*/eboss-0.2/Fits/result_fvoigt_v4.7.22_downsampling0.5.h5"
files['lowcut_nhi20.3'] = "/global/projecta/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output_nhi*/eboss-0.2/Fits/result_fvoigt_v4.7.22_lowcut_nhi20.3.h5"
files['highcut_nhi20.3'] = "/global/projecta/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output_nhi*/eboss-0.2/Fits/result_fvoigt_v4.7.22_highcut_nhi20.3.h5"
files['nhi_bin'] = "/global/projecta/projectdirs/desi/mocks/lya_forest/develop/saclay/v4.7/intermediate_files/mock_22/output_nhi*/eboss-0.2/Fits/result_fvoigt_v4.7.22_nhi_*.h5"

# colors
colors['dr16'] = 'black'
colors['mock-0.0'] = 'royalblue'
colors['mock-0.2'] = 'darkgreen'
colors['v4.7.22'] = 'red'
colors['downsampling0.5'] = 'darkorange'
colors['lowcut_nhi20.3'] = 'darkblue'
colors['highcut_nhi20.3'] = 'darkmagenta'
colors['nhi_bin'] = 'darkorange'

title = '0 < z < 10 - Fvoigt_v4.7.22.txt'
absolute_bias = True

# Read refs
for item in refs:
    print("Reading from {}".format(files[item]))
    pars = util.extract_h5file(files[item])
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

# Read bins
for item in mocks:
    files_mock = glob.glob(files[item])
    print("Reading mocks from {}".format(files_mock))
    nhi[item] = []
    beff_lya[item] = []
    beff_err_lya[item] = []
    beta_lya[item] = []
    beta_err_lya[item] = []
    b_hcd[item] = []
    b_err_hcd[item] = []
    beta_hcd[item] = []
    beta_err_hcd[item] = []
    
    for f in files_mock:
        idx = f.find('output_nhi') + 11
        nhi[item].append(np.mean([float(f[idx:idx+4]), float(f[idx+5:idx+9])]))
        pars = util.extract_h5file(f)
        beff_lya[item].append(pars[2]['beff_LYA'])
        beff_err_lya[item].append(pars[3]['beff_LYA'])
        beta_lya[item].append(pars[2]['beta_LYA'])
        beta_err_lya[item].append(pars[3]['beta_LYA'])
        b_hcd[item].append(pars[2]['bias_hcd'])
        b_err_hcd[item].append(pars[3]['bias_hcd'])
        beta_hcd[item].append(pars[2]['beta_hcd'])
        beta_err_hcd[item].append(pars[3]['beta_hcd'])

# Take absolute values of biases
if absolute_bias:
    for item in mocks+refs:
        beff_lya[item] = np.abs(beff_lya)
        b_hcd[item] = np.abs(b_hcd)

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
ax1.set_xlabel('log(n_HI)')
if absolute_bias:
    ax1.set_ylabel(r'$|beff_{LYA}|$')
else:
    ax1.set_ylabel(r'$beff_{LYA}$')
ax1.set_title(title)
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
ax1.set_xlabel('log(n_HI)')
ax1.set_ylabel(r'$\beta_{LYA}$')
ax1.set_title(title)
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
ax1.set_xlabel('log(n_HI)')
if absolute_bias:
    ax1.set_ylabel(r'$|b_{HCD}|$')
else:
    ax1.set_ylabel(r'$b_{HCD}$')
ax1.set_title(title)
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
ax1.set_xlabel('log(n_HI)')
ax1.set_ylabel(r'$\beta_{HCD}$')
ax1.set_title(title)
ax1.legend()
ax1.grid()
plt.tight_layout()

plt.show()
