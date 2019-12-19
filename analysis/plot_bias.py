import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from SaclayMocks import util


# Options
plot_bias = True
plot_beta = True
plot_bias_eta = True
plot_beff = True

toplot = ['mock1', 'mock2', 'mock3', 'mock_raw', 'data']
colors = ['darkorange', 'red', 'magenta', 'green', 'royalblue']

plot_data = True
plot_mock1 = True
plot_mock2 = True
plot_mock3 = True
plot_mockmean = True
plot_mock_raw = True
plot_pred1 = True
plot_pred2 = True
plot_shades = False

correct_beff = True

def f(x, a, gamma):
    return a*(1+x)**gamma

if correct_beff:
    z0 = 2.4
    print("beff is corrected by G(zz)/G(z0) with z0 = {}".format(z0))

redshift = {}
bias_eta = {}
bias_eta_err = {}
beta = {}
beta_err = {}
bias = {}
bias_err = {}
cor = {}
beff = {}
beff_err = {}
p_bias = {}
p_beta = {}
p_bias_eta = {}
p_beff = {}

# les mocks avec quickquasars (distorsion matrix) + DLAs (v4.7.22)
redshift['mock1'] = np.array([2.09, 2.21, 2.52, 2.85])
bias_eta['mock1'] = np.array([-0.172, -0.195, -0.226, -0.267])
bias_eta_err['mock1'] = np.array([0.003, 0.004, 0.007, 0.017])
beta['mock1'] = np.array([1.55, 1.49, 1.07, 1.11])
beta_err['mock1'] = np.array([0.07, 0.07, 0.07, 0.16])
cor['mock1'] = -0.87
# bias['mock1'] = np.array([-0.107, -0.127, -0.204, -0.233])
# bias_err['mock1'] = np.array([0.005, 0.006, 0.009, 0.025])
bias['mock1'] = bias_eta['mock1'] * 0.97 / beta['mock1']
bias_err['mock1'] = util.bias_err(bias_eta['mock1'], bias_eta_err['mock1'], beta['mock1'], beta_err['mock1'], cor['mock1'])
beff['mock1'] = bias['mock1'] * np.sqrt(1+2/3*beta['mock1']+1/5*beta['mock1']**2)
if correct_beff:
    beff['mock1'] *= (1+z0) / (1+redshift['mock1'])
beff_err['mock1'] = bias_err['mock1'] 

print("Fits on mock1:")
p_bias['mock1'] = sp.optimize.curve_fit(f, redshift['mock1'], bias['mock1'], sigma=bias_err['mock1'])
print("bias(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_bias['mock1'][0][0], p_bias['mock1'][1][0,0], p_bias['mock1'][0][1], p_bias['mock1'][1][1,1]))
p_beta['mock1'] = sp.optimize.curve_fit(f, redshift['mock1'], beta['mock1'], sigma=beta_err['mock1'])
print("beta(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_beta['mock1'][0][0], p_beta['mock1'][1][0,0], p_beta['mock1'][0][1], p_beta['mock1'][1][1,1]))
p_bias_eta['mock1'] = sp.optimize.curve_fit(f, redshift['mock1'], bias_eta['mock1'], sigma=bias_eta_err['mock1'])
print("bias_eta(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_bias_eta['mock1'][0][0], p_bias_eta['mock1'][1][0,0], p_bias_eta['mock1'][0][1], p_bias_eta['mock1'][1][1,1]))
p_beff['mock1'] = sp.optimize.curve_fit(f, redshift['mock1'], beff['mock1'], sigma=beff_err['mock1'])
print("beff(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_beff['mock1'][0][0], p_beff['mock1'][1][0,0], p_beff['mock1'][0][1], p_beff['mock1'][1][1,1]))


# les mocks avec quickquasars + masked DLAs log(n_HI) > 20 (v4.7.22)
redshift['zmock2'] = np.array([2.09, 2.21, 2.52, 2.85])
bias_eta['mock2'] = np.array([-0.174, -0.195, -0.227, -0.279])
bias_eta_err['mock2'] = np.array([0.003, 0.004, 0.006, 0.018])
beta['mock2'] = np.array([1.69, 1.58, 1.16, 1.29])
beta_err['mock2'] = np.array([0.08, 0.07, 0.07, 0.21])
bias['mock2'] = bias_eta['mock2'] * 0.97 / beta['mock2']
cor['mock2'] = -0.87
bias_err['mock2'] = util.bias_err(bias_eta['mock2'], bias_eta_err['mock2'], beta['mock2'], beta_err['mock2'], cor['mock2'])
beff['mock2'] = bias['mock2'] * np.sqrt(1+2/3*beta['mock2']+1/5*beta['mock2']**2)
if correct_beff:
    beff['mock2'] *= (1+z0) / (1+redshift['zmock2'])
beff_err['mock2'] = bias_err['mock2']

print("Fits on mock2:")
p_bias['mock2'] = sp.optimize.curve_fit(f, redshift['zmock2'], bias['mock2'], sigma=bias_err['mock2'])
print("bias(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_bias['mock2'][0][0], p_bias['mock2'][1][0,0], p_bias['mock2'][0][1], p_bias['mock2'][1][1,1]))
p_beta['mock2'] = sp.optimize.curve_fit(f, redshift['zmock2'], beta['mock2'], sigma=beta_err['mock2'])
print("beta(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_beta['mock2'][0][0], p_beta['mock2'][1][0,0], p_beta['mock2'][0][1], p_beta['mock2'][1][1,1]))
p_bias_eta['mock2'] = sp.optimize.curve_fit(f, redshift['zmock2'], bias_eta['mock2'], sigma=bias_eta_err['mock2'])
print("bias_eta(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_bias_eta['mock2'][0][0], p_bias_eta['mock2'][1][0,0], p_bias_eta['mock2'][0][1], p_bias_eta['mock2'][1][1,1]))
p_beff['mock2'] = sp.optimize.curve_fit(f, redshift['zmock2'], beff['mock2'], sigma=beff_err['mock2'])
print("beff(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_beff['mock2'][0][0], p_beff['mock2'][1][0,0], p_beff['mock2'][0][1], p_beff['mock2'][1][1,1]))


# les mocks qvec quickquasars + masked DLAs with log(n_HI) > 20 (v4.7.27)
redshift['zmock3'] = np.array([2.09, 2.21, 2.52, 2.85])
bias_eta['mock3'] = np.array([-0.169, -0.196, -0.235, -0.277])
bias_eta_err['mock3'] = np.array([0.003, 0.003, 0.007, 0.018])
beta['mock3'] = np.array([1.56, 1.61, 1.32, 1.03])
beta_err['mock3'] = np.array([0.07, 0.07, 0.09, 0.15])
bias['mock3'] = bias_eta['mock3'] * 0.97 / beta['mock3']
cor['mock3'] = -0.87
bias_err['mock3'] = util.bias_err(bias_eta['mock3'], bias_eta_err['mock3'], beta['mock3'], beta_err['mock3'], cor['mock3'])
beff['mock3'] = bias['mock3'] * np.sqrt(1+2/3*beta['mock3']+1/5*beta['mock3']**2)
if correct_beff:
    beff['mock3'] *= (1+z0) / (1+redshift['zmock3'])
beff_err['mock3'] = bias_err['mock3']

print("Fits on mock3:")
p_bias['mock3'] = sp.optimize.curve_fit(f, redshift['zmock3'], bias['mock3'], sigma=bias_err['mock3'])
print("bias(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_bias['mock3'][0][0], p_bias['mock3'][1][0,0], p_bias['mock3'][0][1], p_bias['mock3'][1][1,1]))
p_beta['mock3'] = sp.optimize.curve_fit(f, redshift['zmock3'], beta['mock3'], sigma=beta_err['mock3'])
print("beta(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_beta['mock3'][0][0], p_beta['mock3'][1][0,0], p_beta['mock3'][0][1], p_beta['mock3'][1][1,1]))
p_bias_eta['mock3'] = sp.optimize.curve_fit(f, redshift['zmock3'], bias_eta['mock3'], sigma=bias_eta_err['mock3'])
print("bias_eta(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_bias_eta['mock3'][0][0], p_bias_eta['mock3'][1][0,0], p_bias_eta['mock3'][0][1], p_bias_eta['mock3'][1][1,1]))
p_beff['mock3'] = sp.optimize.curve_fit(f, redshift['zmock3'], beff['mock3'], sigma=beff_err['mock3'])
print("beff(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_beff['mock3'][0][0], p_beff['mock3'][1][0,0], p_beff['mock3'][0][1], p_beff['mock3'][1][1,1]))


# # la moyenne des mocks
# mocks = ['mock2', 'mock3']
# zmockmean = 


# les mocks directement sur les transmissions (v4.7.22)
redshift['mock_raw'] = np.array([2.10, 2.24, 2.53, 2.87])
bias_eta['mock_raw'] = np.array([-0.1786, -0.2054, -0.244, -0.277])
bias_eta_err['mock_raw'] = np.array([0.002, 0.0022, 0.004, 0.013])
beta['mock_raw'] = np.array([1.73, 1.58, 1.39, 1.07])
beta_err['mock_raw'] = np.array([0.04, 0.03, 0.04, 0.09])
# bias['mock_raw'] = np.array([-0.10013988, -0.1261    , -0.17027338, -0.25111215])
# bias_err['mock_raw'] = np.array([-0.00131754, -0.00119692, -0.00242659, -0.01066821])
bias['mock_raw'] = bias_eta['mock_raw'] * 0.97 / beta['mock_raw']
cor['mock_raw'] = -0.96
bias_err['mock_raw'] = util.bias_err(bias_eta['mock_raw'], bias_eta_err['mock_raw'], beta['mock_raw'], beta_err['mock_raw'], cor['mock_raw'])
beff['mock_raw'] = bias['mock_raw'] * np.sqrt(1+2/3*beta['mock_raw']+1/5*beta['mock_raw']**2)
if correct_beff:
    beff['mock_raw'] *= (1+z0) / (1+redshift['mock_raw'])
beff_err['mock_raw'] = bias_err['mock_raw']

print("Fits on mock_raw:")
p_bias['mock_raw'] = sp.optimize.curve_fit(f, redshift['mock_raw'], bias['mock_raw'], sigma=bias_err['mock_raw'])
print("bias(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_bias['mock_raw'][0][0], p_bias['mock_raw'][1][0,0], p_bias['mock_raw'][0][1], p_bias['mock_raw'][1][1,1]))
p_beta['mock_raw'] = sp.optimize.curve_fit(f, redshift['mock_raw'], beta['mock_raw'], sigma=beta_err['mock_raw'])
print("beta(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_beta['mock_raw'][0][0], p_beta['mock_raw'][1][0,0], p_beta['mock_raw'][0][1], p_beta['mock_raw'][1][1,1]))
p_bias_eta['mock_raw'] = sp.optimize.curve_fit(f, redshift['mock_raw'], bias_eta['mock_raw'], sigma=bias_eta_err['mock_raw'])
print("bias_eta(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_bias_eta['mock_raw'][0][0], p_bias_eta['mock_raw'][1][0,0], p_bias_eta['mock_raw'][0][1], p_bias_eta['mock_raw'][1][1,1]))
p_beff['mock_raw'] = sp.optimize.curve_fit(f, redshift['mock_raw'], beff['mock_raw'], sigma=beff_err['mock_raw'])
print("beff(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_beff['mock_raw'][0][0], p_beff['mock_raw'][1][0,0], p_beff['mock_raw'][0][1], p_beff['mock_raw'][1][1,1]))


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
bias['data'] = bias_eta['data'] * 0.97 / beta['data']
cor['data'] = -0.9
bias_err['data'] = util.bias_err(bias_eta['data'], bias_eta_err['data'], beta['data'], beta_err['data'], cor['data'])
beff['data'] = bias['data'] * np.sqrt(1+2/3*beta['data']+1/5*beta['data']**2)
if correct_beff:
    beff['data'] *= (1+z0) / (1+redshift['data'])
beff_err['data'] = bias_err['data']

print("Fits on data:")
p_bias['data'] = sp.optimize.curve_fit(f, redshift['data'], bias['data'], sigma=bias_err['data'])
print("bias(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_bias['data'][0][0], p_bias['data'][1][0,0], p_bias['data'][0][1], p_bias['data'][1][1,1]))
p_beta['data'] = sp.optimize.curve_fit(f, redshift['data'], beta['data'], sigma=beta_err['data'])
print("beta(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_beta['data'][0][0], p_beta['data'][1][0,0], p_beta['data'][0][1], p_beta['data'][1][1,1]))
p_bias_eta['data'] = sp.optimize.curve_fit(f, redshift['data'], bias_eta['data'], sigma=bias_eta_err['data'])
print("bias_eta(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_bias_eta['data'][0][0], p_bias_eta['data'][1][0,0], p_bias_eta['data'][0][1], p_bias_eta['data'][1][1,1]))
p_beff['data'] = sp.optimize.curve_fit(f, redshift['data'], beff['data'], sigma=beff_err['data'])
print("beff(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_beff['data'][0][0], p_beff['data'][1][0,0], p_beff['data'][0][1], p_beff['data'][1][1,1]))


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


# Plots
if plot_bias_eta:
    fig, ax = plt.subplots()
    for i, item in enumerate(toplot):
        ax.errorbar(redshift[item], bias_eta[item], yerr=bias_eta_err[item], fmt='.', label=item, color=colors[i])
        z = np.linspace(redshift[item].min(), redshift[item].max(), 100)
        plt.plot(z, f(z, p_bias_eta[item][0][0], p_bias_eta[item][0][1]), linestyle='--', color=colors[i])
        if plot_shades:
            plt.fill_between(z, f(z, p_bias_eta[item][0][0]+p_bias_eta[item][1][0,0], p_bias_eta[item][0][1]-p_bias_eta[item][1][1,1]),
                         f(z, p_bias_eta[item][0][0]-p_bias_eta[item][1][0,0], p_bias_eta[item][0][1]+p_bias_eta[item][1][1,1]),
                         color=colors[i], alpha=0.2)
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
        ax.errorbar(redshift[item], beta[item], yerr=beta_err[item], fmt='.', label=item, color=colors[i])
        z = np.linspace(redshift[item].min(), redshift[item].max(), 100)
        plt.plot(z, f(z, p_beta[item][0][0], p_beta[item][0][1]), linestyle='--', color=colors[i])
        if plot_shades:
            plt.fill_between(z, f(z, p_beta[item][0][0]+p_beta[item][1][0,0], p_beta[item][0][1]-p_beta[item][1][1,1]),
                         f(z, p_beta[item][0][0]-p_beta[item][1][0,0], p_beta[item][0][1]+p_beta[item][1][1,1]),
                         color=colors[i], alpha=0.2)
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
        ax.errorbar(redshift[item], bias[item], yerr=bias_err[item], fmt='.', label=item, color=colors[i])
        z = np.linspace(redshift[item].min(), redshift[item].max(), 100)
        plt.plot(z, f(z, p_bias[item][0][0], p_bias[item][0][1]), linestyle='--', color=colors[i])
        if plot_shades:
            plt.fill_between(z, f(z, p_bias[item][0][0]+p_bias[item][1][0,0], p_bias[item][0][1]-p_bias[item][1][1,1]),
                         f(z, p_bias[item][0][0]-p_bias[item][1][0,0], p_bias[item][0][1]+p_bias[item][1][1,1]),
                         color=colors[i], alpha=0.2)
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
        ax.errorbar(redshift[item], beff[item], yerr=beff_err[item], fmt='.', label=item, color=colors[i])
        z = np.linspace(redshift[item].min(), redshift[item].max(), 100)
        plt.plot(z, f(z, p_beff[item][0][0], p_beff[item][0][1]), linestyle='--', color=colors[i])
        if plot_shades:
            plt.fill_between(z, f(z, p_beff[item][0][0]+p_beff[item][1][0,0], p_beff[item][0][1]-p_beff[item][1][1,1]),
                         f(z, p_beff[item][0][0]-p_beff[item][1][0,0], p_beff[item][0][1]+p_beff[item][1][1,1]),
                         color=colors[i], alpha=0.2)
    if plot_pred1:
        ax.plot(zpred1, beff_pred1, 'x', color='darkgreen', label='pred 10<r<180')
    if plot_pred2:
        ax.plot(zpred2, beff_pred2, 'x', color='limegreen', label='pred 40<r<180')
    ax.legend()
    ax.grid()
    ax.set_xlabel('z')
    ax.set_ylabel('beff')


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
#         ax.errorbar(redshift['zmock2'], bias['mock2'], yerr=bias_err['mock2'], fmt='.', label='mock2', color='red')
#         z = np.linspace(redshift['zmock2'].min(), redshift['zmock2'].max(), 100)
#         plt.plot(z, f(z, p_bias['mock2'][0][0], p_bias['mock2'][0][1]), linestyle='--', color='red')
#         if plot_shades:
#             plt.fill_between(z, f(z, p_bias['mock2'][0][0]+p_bias['mock2'][1][0,0], p_bias['mock2'][0][1]-p_bias['mock2'][1][1,1]),
#                          f(z, p_bias['mock2'][0][0]-p_bias['mock2'][1][0,0], p_bias['mock2'][0][1]+p_bias['mock2'][1][1,1]),
#                          color='red', alpha=0.2)
#     if plot_mock3:
#         ax.errorbar(redshift['zmock3'], bias['mock3'], yerr=bias_err['mock3'], fmt='.', label='mock3', color='magenta')
#         z = np.linspace(redshift['zmock3'].min(), redshift['zmock3'].max(), 100)
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
#         ax.errorbar(redshift['zmock2'], beta['mock2'], yerr=beta_err['mock2'], fmt='.', label='mock2', color='red')
#         z = np.linspace(redshift['zmock2'].min(), redshift['zmock2'].max(), 100)
#         plt.plot(z, f(z, p_beta['mock2'][0][0], p_beta['mock2'][0][1]), linestyle='--', color='red')
#         if plot_shades:
#             plt.fill_between(z, f(z, p_beta['mock2'][0][0]+p_beta['mock2'][1][0,0], p_beta['mock2'][0][1]-p_beta['mock2'][1][1,1]),
#                          f(z, p_beta['mock2'][0][0]-p_beta['mock2'][1][0,0], p_beta['mock2'][0][1]+p_beta['mock2'][1][1,1]),
#                          color='red', alpha=0.2)
#     if plot_mock3:
#         ax.errorbar(redshift['zmock3'], beta['mock3'], yerr=beta_err['mock3'], fmt='.', label='mock3', color='magenta')
#         z = np.linspace(redshift['zmock3'].min(), redshift['zmock3'].max(), 100)
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
#         ax.errorbar(redshift['zmock2'], bias_eta['mock2'], yerr=bias_eta_err['mock2'], fmt='.', label='mock2', color='red')
#         z = np.linspace(redshift['zmock2'].min(), redshift['zmock2'].max(), 100)
#         plt.plot(z, f(z, p_bias_eta['mock2'][0][0], p_bias_eta['mock2'][0][1]), linestyle='--', color='red')
#         if plot_shades:
#             plt.fill_between(z, f(z, p_bias_eta['mock2'][0][0]+p_bias_eta['mock2'][1][0,0], p_bias_eta['mock2'][0][1]-p_bias_eta['mock2'][1][1,1]),
#                          f(z, p_bias_eta['mock2'][0][0]-p_bias_eta['mock2'][1][0,0], p_bias_eta['mock2'][0][1]+p_bias_eta['mock2'][1][1,1]),
#                          color='red', alpha=0.2)
#     if plot_mock3:
#         ax.errorbar(redshift['zmock3'], bias_eta['mock3'], yerr=bias_eta_err['mock3'], fmt='.', label='mock3', color='magenta')
#         z = np.linspace(redshift['zmock3'].min(), redshift['zmock3'].max(), 100)
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
#         ax.errorbar(redshift['zmock2'], beff['mock2'], yerr=beff_err['mock2'], fmt='.', label='mock2', color='red')
#         z = np.linspace(redshift['zmock2'].min(), redshift['zmock2'].max(), 100)
#         plt.plot(z, f(z, p_beff['mock2'][0][0], p_beff['mock2'][0][1]), linestyle='--', color='red')
#         if plot_shades:
#             plt.fill_between(z, f(z, p_beff['mock2'][0][0]+p_beff['mock2'][1][0,0], p_beff['mock2'][0][1]-p_beff['mock2'][1][1,1]),
#                          f(z, p_beff['mock2'][0][0]-p_beff['mock2'][1][0,0], p_beff['mock2'][0][1]+p_beff['mock2'][1][1,1]),
#                          color='red', alpha=0.2)
#     if plot_mock3:
#         ax.errorbar(redshift['zmock3'], beff['mock3'], yerr=beff_err['mock3'], fmt='.', label='mock3', color='magenta')
#         z = np.linspace(redshift['zmock3'].min(), redshift['zmock3'].max(), 100)
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
