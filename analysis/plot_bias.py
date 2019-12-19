import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from SaclayMocks import util


### Options
plot_bias = True
plot_beta = True
plot_bias_eta = True
plot_beff = True

# toplot = ['mock1', 'mock2', 'mock3', 'mock_raw', 'data', 'mock_mean']
toplot = ['mock_mean', 'mock_raw', 'data']
mocks = ['mock2', 'mock3']  # mocks used to compute means

# plot_data = True
# plot_mock1 = True
# plot_mock2 = True
# plot_mock3 = True
# plot_mockmean = True
# plot_mock_raw = True
plot_pred1 = True
plot_pred2 = True
plot_shades = False

correct_beff = False

def f(x, a, gamma):
    return a*(1+x)**gamma

if correct_beff:
    z0 = 2.4
    print("beff is corrected by G(zz)/G(z0) with z0 = {}".format(z0))

# create dictionnaries
redshift = {}
bias_eta = {}
bias_eta_err = {}
beta = {}
beta_err = {}
bias = {}
bias_err = {}
cor = {}
colors = {}
beff = {}
beff_err = {}
p_bias = {}
p_beta = {}
p_bias_eta = {}
p_beff = {}

### Input : all the various data sets
# les mocks avec quickquasars (distorsion matrix) + DLAs (v4.7.22)
redshift['mock1'] = np.array([2.09, 2.21, 2.52, 2.85])
bias_eta['mock1'] = np.array([-0.172, -0.195, -0.226, -0.267])
bias_eta_err['mock1'] = np.array([0.003, 0.004, 0.007, 0.017])
beta['mock1'] = np.array([1.55, 1.49, 1.07, 1.11])
beta_err['mock1'] = np.array([0.07, 0.07, 0.07, 0.16])
cor['mock1'] = -0.87
colors['mock1'] = 'darkorange'

# les mocks avec quickquasars + masked DLAs log(n_HI) > 20 (v4.7.22)
redshift['mock2'] = np.array([2.09, 2.21, 2.52, 2.85])
bias_eta['mock2'] = np.array([-0.174, -0.195, -0.227, -0.279])
bias_eta_err['mock2'] = np.array([0.003, 0.004, 0.006, 0.018])
beta['mock2'] = np.array([1.69, 1.58, 1.16, 1.29])
beta_err['mock2'] = np.array([0.08, 0.07, 0.07, 0.21])
cor['mock2'] = -0.87
colors['mock2'] = 'darkviolet'

# les mocks qvec quickquasars + masked DLAs with log(n_HI) > 20 (v4.7.27)
redshift['mock3'] = np.array([2.09, 2.21, 2.52, 2.85])
bias_eta['mock3'] = np.array([-0.169, -0.196, -0.235, -0.277])
bias_eta_err['mock3'] = np.array([0.003, 0.003, 0.007, 0.018])
beta['mock3'] = np.array([1.56, 1.61, 1.32, 1.03])
beta_err['mock3'] = np.array([0.07, 0.07, 0.09, 0.15])
cor['mock3'] = -0.87
colors['mock3'] = 'magenta'

# les mocks directement sur les transmissions (v4.7.22)
redshift['mock_raw'] = np.array([2.10, 2.24, 2.53, 2.87])
bias_eta['mock_raw'] = np.array([-0.1786, -0.2054, -0.244, -0.277])
bias_eta_err['mock_raw'] = np.array([0.002, 0.0022, 0.004, 0.013])
beta['mock_raw'] = np.array([1.73, 1.58, 1.39, 1.07])
beta_err['mock_raw'] = np.array([0.04, 0.03, 0.04, 0.09])
cor['mock_raw'] = -0.96
colors['mock_raw'] = 'green'

# les donnees
data_bias_eta= np.array( [[ 2.13781251, -0.18558474,  0.00650885],
                          [ 2.27259149, -0.20089904,  0.00566513],
                          [ 2.54744584, -0.23754895,  0.00859699],
                          [ 2.9263245,  -0.28984325,  0.01907289]] )
data_beta = np.array([[2.13781251, 1.71518188, 0.14670333],
                      [2.27259149, 1.63599596, 0.12377548],
                      [2.54744584, 1.34603853, 0.10685724],
                      [2.9263245,  1.20116875, 0.16525412]])
redshift['data'] = data_beta[:,0]
bias_eta['data'] = data_bias_eta[:,1]
bias_eta_err['data'] = data_bias_eta[:,2]
beta['data'] = data_beta[:,1]
beta_err['data'] = data_beta[:,2]
cor['data'] = -0.9
colors['data'] = 'royalblue'

# Le fit de la pred, faites dans chaque bin de mock_raw sur 10 < r < 180:
zpred1 = np.array([2.10, 2.24, 2.53, 2.87])
bias_eta_pred1 = np.array([-0.1753, -0.2019, -0.237, -0.284])
beta_pred1 = np.array([1.68, 1.57, 1.37, 1.18])
bias_pred1 = bias_eta_pred1 * 0.97 / beta_pred1
beff_pred1 = bias_pred1 * np.sqrt(1+2/3*beta_pred1+1/5*beta_pred1**2)

# Le fit de la pred, faites dans chaque bin de mock_raw, sur 40 < r < 180 :
zpred2 = np.array([2.10, 2.24, 2.53, 2.87])
bias_eta_pred2 = np.array([-0.180, -0.205, -0.240, -0.290])
beta_pred2 = np.array([1.76, 1.63, 1.40, 1.19])
bias_pred2 = bias_eta_pred2 * 0.97 / beta_pred2
beff_pred2 = bias_pred2 * np.sqrt(1+2/3*beta_pred2+1/5*beta_pred2**2)

# Compute means on mocks
if 'mock_mean' in toplot:
    redshift['mock_mean'] = np.array([redshift[m] for m in mocks]).mean(axis=0)
    bias_eta['mock_mean'] = np.array([bias_eta[m] for m in mocks]).mean(axis=0)  # assume same errors
    bias_eta_err['mock_mean'] = np.array([bias_eta_err[m] for m in mocks]).mean(axis=0)
    bias_eta_err['mock_mean'] /= np.sqrt(len(mocks))
    beta['mock_mean'] = np.array([beta[m] for m in mocks]).mean(axis=0)
    beta_err['mock_mean'] = np.array([beta_err[m] for m in mocks]).mean(axis=0)
    beta_err['mock_mean'] /= np.sqrt(len(mocks))
    cor['mock_mean'] = np.mean([cor[m] for m in mocks])
    colors['mock_mean'] = 'r'

### Compute bias, b_eff, and fit with power law
for item in toplot:
    bias[item] = bias_eta[item] * 0.97 / beta[item]
    bias_err[item] = util.bias_err(bias_eta[item], bias_eta_err[item], beta[item], beta_err[item], cor[item])
    beff[item] = bias[item] * np.sqrt(1+2/3*beta[item]+1/5*beta[item]**2)
    if correct_beff:
        beff[item] *= (1+z0) / (1+redshift[item])
    beff_err[item] = bias_err[item]

    print("Fits on {}:".format(item))
    p_bias[item] = sp.optimize.curve_fit(f, redshift[item], bias[item], sigma=bias_err[item])
    print("bias(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_bias[item][0][0], p_bias[item][1][0,0], p_bias[item][0][1], p_bias[item][1][1,1]))
    p_beta[item] = sp.optimize.curve_fit(f, redshift[item], beta[item], sigma=beta_err[item])
    print("beta(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_beta[item][0][0], p_beta[item][1][0,0], p_beta[item][0][1], p_beta[item][1][1,1]))
    p_bias_eta[item] = sp.optimize.curve_fit(f, redshift[item], bias_eta[item], sigma=bias_eta_err[item])
    print("bias_eta(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_bias_eta[item][0][0], p_bias_eta[item][1][0,0], p_bias_eta[item][0][1], p_bias_eta[item][1][1,1]))
    p_beff[item] = sp.optimize.curve_fit(f, redshift[item], beff[item], sigma=beff_err[item])
    print("beff(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_beff[item][0][0], p_beff[item][1][0,0], p_beff[item][0][1], p_beff[item][1][1,1]))

### Plots
if plot_bias_eta:
    fig, ax = plt.subplots()
    for i, item in enumerate(toplot):
        ax.errorbar(redshift[item], bias_eta[item], yerr=bias_eta_err[item], fmt='.', label=item, color=colors[item])
        z = np.linspace(redshift[item].min(), redshift[item].max(), 100)
        plt.plot(z, f(z, p_bias_eta[item][0][0], p_bias_eta[item][0][1]), linestyle='--', color=colors[item])
        if plot_shades:
            plt.fill_between(z, f(z, p_bias_eta[item][0][0]+p_bias_eta[item][1][0,0], p_bias_eta[item][0][1]-p_bias_eta[item][1][1,1]),
                         f(z, p_bias_eta[item][0][0]-p_bias_eta[item][1][0,0], p_bias_eta[item][0][1]+p_bias_eta[item][1][1,1]),
                         color=colors[item], alpha=0.2)
    if plot_pred1:
        ax.plot(zpred1, bias_eta_pred1, 'x', color='darkgreen', label='pred 10<r<180')
    if plot_pred2:
        ax.plot(zpred2, bias_eta_pred2, 'x', color='limegreen', label='pred 40<r<180')
    ax.legend()
    ax.grid()
    ax.set_xlabel('z')
    ax.set_ylabel('bias_eta')

if plot_beta:
    fig, ax = plt.subplots()
    for i, item in enumerate(toplot):
        ax.errorbar(redshift[item], beta[item], yerr=beta_err[item], fmt='.', label=item, color=colors[item])
        z = np.linspace(redshift[item].min(), redshift[item].max(), 100)
        plt.plot(z, f(z, p_beta[item][0][0], p_beta[item][0][1]), linestyle='--', color=colors[item])
        if plot_shades:
            plt.fill_between(z, f(z, p_beta[item][0][0]+p_beta[item][1][0,0], p_beta[item][0][1]-p_beta[item][1][1,1]),
                         f(z, p_beta[item][0][0]-p_beta[item][1][0,0], p_beta[item][0][1]+p_beta[item][1][1,1]),
                         color=colors[item], alpha=0.2)
    if plot_pred1:
        ax.plot(zpred1, beta_pred1, 'x', color='darkgreen', label='pred 10<r<180')
    if plot_pred2:
        ax.plot(zpred2, beta_pred2, 'x', color='limegreen', label='pred 40<r<180')
    ax.legend()
    ax.grid()
    ax.set_xlabel('z')
    ax.set_ylabel('beta')

if plot_bias:
    fig, ax = plt.subplots()
    for i, item in enumerate(toplot):
        ax.errorbar(redshift[item], bias[item], yerr=bias_err[item], fmt='.', label=item, color=colors[item])
        z = np.linspace(redshift[item].min(), redshift[item].max(), 100)
        plt.plot(z, f(z, p_bias[item][0][0], p_bias[item][0][1]), linestyle='--', color=colors[item])
        if plot_shades:
            plt.fill_between(z, f(z, p_bias[item][0][0]+p_bias[item][1][0,0], p_bias[item][0][1]-p_bias[item][1][1,1]),
                         f(z, p_bias[item][0][0]-p_bias[item][1][0,0], p_bias[item][0][1]+p_bias[item][1][1,1]),
                         color=colors[item], alpha=0.2)
    if plot_pred1:
        ax.plot(zpred1, bias_pred1, 'x', color='darkgreen', label='pred 10<r<180')
    if plot_pred2:
        ax.plot(zpred2, bias_pred2, 'x', color='limegreen', label='pred 40<r<180')
    ax.legend()
    ax.grid()
    ax.set_xlabel('z')
    ax.set_ylabel('bias')

if plot_beff:
    fig, ax = plt.subplots()
    for i, item in enumerate(toplot):
        ax.errorbar(redshift[item], beff[item], yerr=beff_err[item], fmt='.', label=item, color=colors[item])
        z = np.linspace(redshift[item].min(), redshift[item].max(), 100)
        plt.plot(z, f(z, p_beff[item][0][0], p_beff[item][0][1]), linestyle='--', color=colors[item])
        if plot_shades:
            plt.fill_between(z, f(z, p_beff[item][0][0]+p_beff[item][1][0,0], p_beff[item][0][1]-p_beff[item][1][1,1]),
                         f(z, p_beff[item][0][0]-p_beff[item][1][0,0], p_beff[item][0][1]+p_beff[item][1][1,1]),
                         color=colors[item], alpha=0.2)
    if plot_pred1:
        ax.plot(zpred1, beff_pred1, 'x', color='darkgreen', label='pred 10<r<180')
    if plot_pred2:
        ax.plot(zpred2, beff_pred2, 'x', color='limegreen', label='pred 40<r<180')
    ax.legend()
    ax.grid()
    ax.set_xlabel('z')
    ylabel='b_eff'
    if correct_beff:
        ylabel += ' G(z) / G({})'.format(z0)
    ax.set_ylabel(ylabel)


plt.show()

# if plot_bias:
#     fig, ax = plt.subplots()
#     if plot_data:
#         ax.errorbar(redshift['data'], bias['data'], yerr=bias_err['data'], fmt='.', label='data', color='royalblue')
#         z = np.linspace(redshift['data'].min(), redshift['data'].max(), 100)
#         plt.plot(z, f(z, p_bias['data'][0][0], p_bias['data'][0][1]), linestyle='--', color='royalblue')
#         if plot_shades:
#             plt.fill_between(z, f(z, p_bias['data'][0][0]+p_bias['data'][1][0,0], p_bias['data'][0][1]-p_bias['data'][1][1,1]),
#                          f(z, p_bias['data'][0][0]-p_bias['data'][1][0,0], p_bias['data'][0][1]+p_bias['data'][1][1,1]),
#                          color='royalblue', alpha=0.2)
#     if plot_mock1:
#         ax.errorbar(redshift['mock1'], bias['mock1'], yerr=bias_err['mock1'], fmt='.', label='mock1', color='darkorange')
#         z = np.linspace(redshift['mock1'].min(), redshift['mock1'].max(), 100)
#         plt.plot(z, f(z, p_bias['mock1'][0][0], p_bias['mock1'][0][1]), linestyle='--', color='darkorange')
#         if plot_shades:
#             plt.fill_between(z, f(z, p_bias['mock1'][0][0]+p_bias['mock1'][1][0,0], p_bias['mock1'][0][1]-p_bias['mock1'][1][1,1]),
#                          f(z, p_bias['mock1'][0][0]-p_bias['mock1'][1][0,0], p_bias['mock1'][0][1]+p_bias['mock1'][1][1,1]),
#                          color='darkorange', alpha=0.2)
#     if plot_mock2:
#         ax.errorbar(redshift['mock2'], bias['mock2'], yerr=bias_err['mock2'], fmt='.', label='mock2', color='red')
#         z = np.linspace(redshift['mock2'].min(), redshift['mock2'].max(), 100)
#         plt.plot(z, f(z, p_bias['mock2'][0][0], p_bias['mock2'][0][1]), linestyle='--', color='red')
#         if plot_shades:
#             plt.fill_between(z, f(z, p_bias['mock2'][0][0]+p_bias['mock2'][1][0,0], p_bias['mock2'][0][1]-p_bias['mock2'][1][1,1]),
#                          f(z, p_bias['mock2'][0][0]-p_bias['mock2'][1][0,0], p_bias['mock2'][0][1]+p_bias['mock2'][1][1,1]),
#                          color='red', alpha=0.2)
#     if plot_mock3:
#         ax.errorbar(redshift['mock3'], bias['mock3'], yerr=bias_err['mock3'], fmt='.', label='mock3', color='magenta')
#         z = np.linspace(redshift['mock3'].min(), redshift['mock3'].max(), 100)
#         plt.plot(z, f(z, p_bias['mock3'][0][0], p_bias['mock3'][0][1]), linestyle='--', color='magenta')
#         if plot_shades:
#             plt.fill_between(z, f(z, p_bias['mock3'][0][0]+p_bias['mock3'][1][0,0], p_bias['mock3'][0][1]-p_bias['mock3'][1][1,1]),
#                          f(z, p_bias['mock3'][0][0]-p_bias['mock3'][1][0,0], p_bias['mock3'][0][1]+p_bias['mock3'][1][1,1]),
#                          color='magenta', alpha=0.2)
#     if plot_mock_raw:
#         ax.errorbar(redshift['mock_raw'], bias['mock_raw'], yerr=bias_err['mock_raw'], fmt='.', label='mock_raw', color='green')
#         z = np.linspace(redshift['mock_raw'].min(), redshift['mock_raw'].max(), 100)
#         plt.plot(z, f(z, p_bias['mock_raw'][0][0], p_bias['mock_raw'][0][1]), linestyle='--', color='green')
#         if plot_shades:
#             plt.fill_between(z, f(z, p_bias['mock_raw'][0][0]+p_bias['mock_raw'][1][0,0], p_bias['mock_raw'][0][1]-p_bias['mock_raw'][1][1,1]),
#                          f(z, p_bias['mock_raw'][0][0]-p_bias['mock_raw'][1][0,0], p_bias['mock_raw'][0][1]+p_bias['mock_raw'][1][1,1]),
#                          color='green', alpha=0.2)
#     if plot_pred1:
#         ax.plot(zpred1, bias_pred1, 'x', color='darkgreen', label='pred 10<r<180')
#     if plot_pred2:
#         ax.plot(zpred2, bias_pred2, 'x', color='limegreen', label='pred 40<r<180')

#     ax.legend()
#     ax.grid()
#     ax.set_xlabel('z')
#     ax.set_ylabel('bias')

# if plot_beta:
#     fig, ax = plt.subplots()
#     if plot_data:
#         ax.errorbar(redshift['data'], beta['data'], yerr=beta_err['data'], fmt='.', label='data', color='royalblue')
#         z = np.linspace(redshift['data'].min(), redshift['data'].max(), 100)
#         plt.plot(z, f(z, p_beta['data'][0][0], p_beta['data'][0][1]), linestyle='--', color='royalblue')
#         if plot_shades:
#             plt.fill_between(z, f(z, p_beta['data'][0][0]+p_beta['data'][1][0,0], p_beta['data'][0][1]-p_beta['data'][1][1,1]),
#                          f(z, p_beta['data'][0][0]-p_beta['data'][1][0,0], p_beta['data'][0][1]+p_beta['data'][1][1,1]),
#                          color='royalblue', alpha=0.2)
#     if plot_mock1:
#         ax.errorbar(redshift['mock1'], beta['mock1'], yerr=beta_err['mock1'], fmt='.', label='mock1', color='darkorange')
#         z = np.linspace(redshift['mock1'].min(), redshift['mock1'].max(), 100)
#         plt.plot(z, f(z, p_beta['mock1'][0][0], p_beta['mock1'][0][1]), linestyle='--', color='darkorange')
#         if plot_shades:
#             plt.fill_between(z, f(z, p_beta['mock1'][0][0]+p_beta['mock1'][1][0,0], p_beta['mock1'][0][1]-p_beta['mock1'][1][1,1]),
#                          f(z, p_beta['mock1'][0][0]-p_beta['mock1'][1][0,0], p_beta['mock1'][0][1]+p_beta['mock1'][1][1,1]),
#                          color='darkorange', alpha=0.2)
#     if plot_mock2:
#         ax.errorbar(redshift['mock2'], beta['mock2'], yerr=beta_err['mock2'], fmt='.', label='mock2', color='red')
#         z = np.linspace(redshift['mock2'].min(), redshift['mock2'].max(), 100)
#         plt.plot(z, f(z, p_beta['mock2'][0][0], p_beta['mock2'][0][1]), linestyle='--', color='red')
#         if plot_shades:
#             plt.fill_between(z, f(z, p_beta['mock2'][0][0]+p_beta['mock2'][1][0,0], p_beta['mock2'][0][1]-p_beta['mock2'][1][1,1]),
#                          f(z, p_beta['mock2'][0][0]-p_beta['mock2'][1][0,0], p_beta['mock2'][0][1]+p_beta['mock2'][1][1,1]),
#                          color='red', alpha=0.2)
#     if plot_mock3:
#         ax.errorbar(redshift['mock3'], beta['mock3'], yerr=beta_err['mock3'], fmt='.', label='mock3', color='magenta')
#         z = np.linspace(redshift['mock3'].min(), redshift['mock3'].max(), 100)
#         plt.plot(z, f(z, p_beta['mock3'][0][0], p_beta['mock3'][0][1]), linestyle='--', color='magenta')
#         if plot_shades:
#             plt.fill_between(z, f(z, p_beta['mock3'][0][0]+p_beta['mock3'][1][0,0], p_beta['mock3'][0][1]-p_beta['mock3'][1][1,1]),
#                          f(z, p_beta['mock3'][0][0]-p_beta['mock3'][1][0,0], p_beta['mock3'][0][1]+p_beta['mock3'][1][1,1]),
#                          color='magenta', alpha=0.2)
#     if plot_mock_raw:
#         ax.errorbar(redshift['mock_raw'], beta['mock_raw'], yerr=beta_err['mock_raw'], fmt='.', label='mock_raw', color='green')
#         z = np.linspace(redshift['mock_raw'].min(), redshift['mock_raw'].max(), 100)
#         plt.plot(z, f(z, p_beta['mock_raw'][0][0], p_beta['mock_raw'][0][1]), linestyle='--', color='green')
#         if plot_shades:
#             plt.fill_between(z, f(z, p_beta['mock_raw'][0][0]+p_beta['mock_raw'][1][0,0], p_beta['mock_raw'][0][1]-p_beta['mock_raw'][1][1,1]),
#                          f(z, p_beta['mock_raw'][0][0]-p_beta['mock_raw'][1][0,0], p_beta['mock_raw'][0][1]+p_beta['mock_raw'][1][1,1]),
#                          color='green', alpha=0.2)
#     if plot_pred1:
#         ax.plot(zpred1, beta_pred1, 'x', color='darkgreen', label='pred 10<r<180')
#     if plot_pred2:
#         ax.plot(zpred2, beta_pred2, 'x', color='limegreen', label='pred 40<r<180')

#     ax.legend()
#     ax.grid()
#     ax.set_xlabel('z')
#     ax.set_ylabel('beta')

# if plot_bias_eta:
#     fig, ax = plt.subplots()
#     if plot_data:
#         ax.errorbar(redshift['data'], bias_eta['data'], yerr=bias_eta_err['data'], fmt='.', label='data', color='royalblue')
#         z = np.linspace(redshift['data'].min(), redshift['data'].max(), 100)
#         plt.plot(z, f(z, p_bias_eta['data'][0][0], p_bias_eta['data'][0][1]), linestyle='--', color='royalblue')
#         if plot_shades:
#             plt.fill_between(z, f(z, p_bias_eta['data'][0][0]+p_bias_eta['data'][1][0,0], p_bias_eta['data'][0][1]-p_bias_eta['data'][1][1,1]),
#                          f(z, p_bias_eta['data'][0][0]-p_bias_eta['data'][1][0,0], p_bias_eta['data'][0][1]+p_bias_eta['data'][1][1,1]),
#                          color='royalblue', alpha=0.2)
#     if plot_mock1:
#         ax.errorbar(redshift['mock1'], bias_eta['mock1'], yerr=bias_eta_err['mock1'], fmt='.', label='mock1', color='darkorange')
#         z = np.linspace(redshift['mock1'].min(), redshift['mock1'].max(), 100)
#         plt.plot(z, f(z, p_bias_eta['mock1'][0][0], p_bias_eta['mock1'][0][1]), linestyle='--', color='darkorange')
#         if plot_shades:
#             plt.fill_between(z, f(z, p_bias_eta['mock1'][0][0]+p_bias_eta['mock1'][1][0,0], p_bias_eta['mock1'][0][1]-p_bias_eta['mock1'][1][1,1]),
#                          f(z, p_bias_eta['mock1'][0][0]-p_bias_eta['mock1'][1][0,0], p_bias_eta['mock1'][0][1]+p_bias_eta['mock1'][1][1,1]),
#                          color='darkorange', alpha=0.2)
#     if plot_mock2:
#         ax.errorbar(redshift['mock2'], bias_eta['mock2'], yerr=bias_eta_err['mock2'], fmt='.', label='mock2', color='red')
#         z = np.linspace(redshift['mock2'].min(), redshift['mock2'].max(), 100)
#         plt.plot(z, f(z, p_bias_eta['mock2'][0][0], p_bias_eta['mock2'][0][1]), linestyle='--', color='red')
#         if plot_shades:
#             plt.fill_between(z, f(z, p_bias_eta['mock2'][0][0]+p_bias_eta['mock2'][1][0,0], p_bias_eta['mock2'][0][1]-p_bias_eta['mock2'][1][1,1]),
#                          f(z, p_bias_eta['mock2'][0][0]-p_bias_eta['mock2'][1][0,0], p_bias_eta['mock2'][0][1]+p_bias_eta['mock2'][1][1,1]),
#                          color='red', alpha=0.2)
#     if plot_mock3:
#         ax.errorbar(redshift['mock3'], bias_eta['mock3'], yerr=bias_eta_err['mock3'], fmt='.', label='mock3', color='magenta')
#         z = np.linspace(redshift['mock3'].min(), redshift['mock3'].max(), 100)
#         plt.plot(z, f(z, p_bias_eta['mock3'][0][0], p_bias_eta['mock3'][0][1]), linestyle='--', color='magenta')
#         if plot_shades:
#             plt.fill_between(z, f(z, p_bias_eta['mock3'][0][0]+p_bias_eta['mock3'][1][0,0], p_bias_eta['mock3'][0][1]-p_bias_eta['mock3'][1][1,1]),
#                          f(z, p_bias_eta['mock3'][0][0]-p_bias_eta['mock3'][1][0,0], p_bias_eta['mock3'][0][1]+p_bias_eta['mock3'][1][1,1]),
#                          color='magenta', alpha=0.2)
#     if plot_mock_raw:
#         ax.errorbar(redshift['mock_raw'], bias_eta['mock_raw'], yerr=bias_eta_err['mock_raw'], fmt='.', label='mock_raw', color='green')
#         z = np.linspace(redshift['mock_raw'].min(), redshift['mock_raw'].max(), 100)
#         plt.plot(z, f(z, p_bias_eta['mock_raw'][0][0], p_bias_eta['mock_raw'][0][1]), linestyle='--', color='green')
#         if plot_shades:
#             plt.fill_between(z, f(z, p_bias_eta['mock_raw'][0][0]+p_bias_eta['mock_raw'][1][0,0], p_bias_eta['mock_raw'][0][1]-p_bias_eta['mock_raw'][1][1,1]),
#                          f(z, p_bias_eta['mock_raw'][0][0]-p_bias_eta['mock_raw'][1][0,0], p_bias_eta['mock_raw'][0][1]+p_bias_eta['mock_raw'][1][1,1]),
#                          color='green', alpha=0.2)
#     if plot_pred1:
#         ax.plot(zpred1, bias_eta_pred1, 'x', color='darkgreen', label='pred 10<r<180')
#     if plot_pred2:
#         ax.plot(zpred2, bias_eta_pred2, 'x', color='limegreen', label='pred 40<r<180')

#     ax.legend()
#     ax.grid()
#     ax.set_xlabel('z')
#     ax.set_ylabel('bias_eta')

# if plot_beff:
#     fig, ax = plt.subplots()
#     if plot_data:
#         ax.errorbar(redshift['data'], beff['data'], yerr=beff_err['data'], fmt='.', label='data', color='royalblue')
#         z = np.linspace(redshift['data'].min(), redshift['data'].max(), 100)
#         plt.plot(z, f(z, p_beff['data'][0][0], p_beff['data'][0][1]), linestyle='--', color='royalblue')
#         if plot_shades:
#             plt.fill_between(z, f(z, p_beff['data'][0][0]+p_beff['data'][1][0,0], p_beff['data'][0][1]-p_beff['data'][1][1,1]),
#                          f(z, p_beff['data'][0][0]-p_beff['data'][1][0,0], p_beff['data'][0][1]+p_beff['data'][1][1,1]),
#                          color='royalblue', alpha=0.2)
#     if plot_mock1:
#         ax.errorbar(redshift['mock1'], beff['mock1'], yerr=beff_err['mock1'], fmt='.', label='mock1', color='darkorange')
#         z = np.linspace(redshift['mock1'].min(), redshift['mock1'].max(), 100)
#         plt.plot(z, f(z, p_beff['mock1'][0][0], p_beff['mock1'][0][1]), linestyle='--', color='darkorange')
#         if plot_shades:
#             plt.fill_between(z, f(z, p_beff['mock1'][0][0]+p_beff['mock1'][1][0,0], p_beff['mock1'][0][1]-p_beff['mock1'][1][1,1]),
#                          f(z, p_beff['mock1'][0][0]-p_beff['mock1'][1][0,0], p_beff['mock1'][0][1]+p_beff['mock1'][1][1,1]),
#                          color='darkorange', alpha=0.2)
#     if plot_mock2:
#         ax.errorbar(redshift['mock2'], beff['mock2'], yerr=beff_err['mock2'], fmt='.', label='mock2', color='red')
#         z = np.linspace(redshift['mock2'].min(), redshift['mock2'].max(), 100)
#         plt.plot(z, f(z, p_beff['mock2'][0][0], p_beff['mock2'][0][1]), linestyle='--', color='red')
#         if plot_shades:
#             plt.fill_between(z, f(z, p_beff['mock2'][0][0]+p_beff['mock2'][1][0,0], p_beff['mock2'][0][1]-p_beff['mock2'][1][1,1]),
#                          f(z, p_beff['mock2'][0][0]-p_beff['mock2'][1][0,0], p_beff['mock2'][0][1]+p_beff['mock2'][1][1,1]),
#                          color='red', alpha=0.2)
#     if plot_mock3:
#         ax.errorbar(redshift['mock3'], beff['mock3'], yerr=beff_err['mock3'], fmt='.', label='mock3', color='magenta')
#         z = np.linspace(redshift['mock3'].min(), redshift['mock3'].max(), 100)
#         plt.plot(z, f(z, p_beff['mock3'][0][0], p_beff['mock3'][0][1]), linestyle='--', color='magenta')
#         if plot_shades:
#             plt.fill_between(z, f(z, p_beff['mock3'][0][0]+p_beff['mock3'][1][0,0], p_beff['mock3'][0][1]-p_beff['mock3'][1][1,1]),
#                          f(z, p_beff['mock3'][0][0]-p_beff['mock3'][1][0,0], p_beff['mock3'][0][1]+p_beff['mock3'][1][1,1]),
#                          color='magenta', alpha=0.2)
#     if plot_mock_raw:
#         ax.errorbar(redshift['mock_raw'], beff['mock_raw'], yerr=beff_err['mock_raw'], fmt='.', label='mock_raw', color='green')
#         z = np.linspace(redshift['mock_raw'].min(), redshift['mock_raw'].max(), 100)
#         plt.plot(z, f(z, p_beff['mock_raw'][0][0], p_beff['mock_raw'][0][1]), linestyle='--', color='green')
#         if plot_shades:
#             plt.fill_between(z, f(z, p_beff['mock_raw'][0][0]+p_beff['mock_raw'][1][0,0], p_beff['mock_raw'][0][1]-p_beff['mock_raw'][1][1,1]),
#                          f(z, p_beff['mock_raw'][0][0]-p_beff['mock_raw'][1][0,0], p_beff['mock_raw'][0][1]+p_beff['mock_raw'][1][1,1]),
#                          color='green', alpha=0.2)
#     if plot_pred1:
#         ax.plot(zpred1, beff_pred1, 'x', color='darkgreen', label='pred 10<r<180')
#     if plot_pred2:
#         ax.plot(zpred2, beff_pred2, 'x', color='limegreen', label='pred 40<r<180')

#     ax.legend()
#     ax.grid()
#     ax.set_xlabel('z')
#     ylabel='b_eff'
#     if correct_beff:
#         ylabel += ' G(z) / G({})'.format(z0)
#     ax.set_ylabel(ylabel)

# plt.show()
