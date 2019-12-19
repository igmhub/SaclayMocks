import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from SaclayMocks import util


# Options
plot_bias = True
plot_beta = True
plot_bias_eta = True
plot_beff = True

plot_data = True
plot_mock1 = True
plot_mock2 = True
plot_mock_raw = True
plot_pred1 = True
plot_pred2 = True
plot_shades = False

correct_beff = False

def f(x, a, gamma):
    return a*(1+x)**gamma

if correct_beff:
    z0 = 2.4
    print("beff is corrected by G(zz)/G(z0) with z0 = {}".format(z0))

# les mocks avec quickquasars (distorsion matrix) + DLAs
zmock1 = np.array([2.09, 2.21, 2.52, 2.85])
bias_eta_mock1 = np.array([-0.172, -0.195, -0.226, -0.267])
bias_eta_err_mock1 = np.array([0.003, 0.004, 0.007, 0.017])
beta_mock1 = np.array([1.55, 1.49, 1.07, 1.11])
beta_err_mock1 = np.array([0.07, 0.07, 0.07, 0.16])
# bias_mock1 = np.array([-0.107, -0.127, -0.204, -0.233])
# bias_err_mock1 = np.array([0.005, 0.006, 0.009, 0.025])
bias_mock1 = bias_eta_mock1 * 0.97 / beta_mock1
cor_mock1 = -0.87
bias_err_mock1 = util.bias_err(bias_eta_mock1, bias_eta_err_mock1, beta_mock1, beta_err_mock1, cor_mock1)
beff_mock1 = bias_mock1 * np.sqrt(1+2/3*beta_mock1+1/5*beta_mock1**2)
if correct_beff:
    beff_mock1 *= (1+z0) / (1+zmock1)
beff_err_mock1 = bias_err_mock1 

print("Fits on mock1:")
p_bias_mock1 = sp.optimize.curve_fit(f, zmock1, bias_mock1, sigma=bias_err_mock1)
print("bias(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_bias_mock1[0][0], p_bias_mock1[1][0,0], p_bias_mock1[0][1], p_bias_mock1[1][1,1]))
p_beta_mock1 = sp.optimize.curve_fit(f, zmock1, beta_mock1, sigma=beta_err_mock1)
print("beta(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_beta_mock1[0][0], p_beta_mock1[1][0,0], p_beta_mock1[0][1], p_beta_mock1[1][1,1]))
p_bias_eta_mock1 = sp.optimize.curve_fit(f, zmock1, bias_eta_mock1, sigma=bias_eta_err_mock1)
print("bias_eta(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_bias_eta_mock1[0][0], p_bias_eta_mock1[1][0,0], p_bias_eta_mock1[0][1], p_bias_eta_mock1[1][1,1]))
p_beff_mock1 = sp.optimize.curve_fit(f, zmock1, beff_mock1, sigma=beff_err_mock1)
print("beff(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_beff_mock1[0][0], p_beff_mock1[1][0,0], p_beff_mock1[0][1], p_beff_mock1[1][1,1]))


# les mocks avec quickquasars + masked DLAs log(n_HI) > 20
zmock2 = np.array([2.09, 2.21, 2.52, 2.85])
bias_eta_mock2 = np.array([-0.174, -0.195, -0.227, -0.279])
bias_eta_err_mock2 = np.array([0.003, 0.004, 0.006, 0.018])
beta_mock2 = np.array([1.69, 1.58, 1.16, 1.29])
beta_err_mock2 = np.array([0.08, 0.07, 0.07, 0.21])
bias_mock2 = bias_eta_mock2 * 0.97 / beta_mock2
cor_mock2 = -0.87
bias_err_mock2 = util.bias_err(bias_eta_mock2, bias_eta_err_mock2, beta_mock2, beta_err_mock2, cor_mock2)
beff_mock2 = bias_mock2 * np.sqrt(1+2/3*beta_mock2+1/5*beta_mock2**2)
if correct_beff:
    beff_mock2 *= (1+z0) / (1+zmock2)
beff_err_mock2 = bias_err_mock2

print("Fits on mock2:")
p_bias_mock2 = sp.optimize.curve_fit(f, zmock2, bias_mock2, sigma=bias_err_mock2)
print("bias(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_bias_mock2[0][0], p_bias_mock2[1][0,0], p_bias_mock2[0][1], p_bias_mock2[1][1,1]))
p_beta_mock2 = sp.optimize.curve_fit(f, zmock2, beta_mock2, sigma=beta_err_mock2)
print("beta(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_beta_mock2[0][0], p_beta_mock2[1][0,0], p_beta_mock2[0][1], p_beta_mock2[1][1,1]))
p_bias_eta_mock2 = sp.optimize.curve_fit(f, zmock2, bias_eta_mock2, sigma=bias_eta_err_mock2)
print("bias_eta(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_bias_eta_mock2[0][0], p_bias_eta_mock2[1][0,0], p_bias_eta_mock2[0][1], p_bias_eta_mock2[1][1,1]))
p_beff_mock2 = sp.optimize.curve_fit(f, zmock2, beff_mock2, sigma=beff_err_mock2)
print("beff(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_beff_mock2[0][0], p_beff_mock2[1][0,0], p_beff_mock2[0][1], p_beff_mock2[1][1,1]))


# les mocks directement sur les transmissions
zmock_raw = np.array([2.10, 2.24, 2.53, 2.87])
bias_eta_mock_raw = np.array([-0.1786, -0.2054, -0.244, -0.277])
bias_eta_err_mock_raw = np.array([0.002, 0.0022, 0.004, 0.013])
beta_mock_raw = np.array([1.73, 1.58, 1.39, 1.07])
beta_err_mock_raw = np.array([0.04, 0.03, 0.04, 0.09])
# bias_mock_raw = np.array([-0.10013988, -0.1261    , -0.17027338, -0.25111215])
# bias_err_mock_raw = np.array([-0.00131754, -0.00119692, -0.00242659, -0.01066821])
bias_mock_raw = bias_eta_mock_raw * 0.97 / beta_mock_raw
cor_mock_raw = -0.96
bias_err_mock_raw = util.bias_err(bias_eta_mock_raw, bias_eta_err_mock_raw, beta_mock_raw, beta_err_mock_raw, cor_mock_raw)
beff_mock_raw = bias_mock_raw * np.sqrt(1+2/3*beta_mock_raw+1/5*beta_mock_raw**2)
if correct_beff:
    beff_mock_raw *= (1+z0) / (1+zmock_raw)
beff_err_mock_raw = bias_err_mock_raw

print("Fits on mock_raw:")
p_bias_mock_raw = sp.optimize.curve_fit(f, zmock_raw, bias_mock_raw, sigma=bias_err_mock_raw)
print("bias(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_bias_mock_raw[0][0], p_bias_mock_raw[1][0,0], p_bias_mock_raw[0][1], p_bias_mock_raw[1][1,1]))
p_beta_mock_raw = sp.optimize.curve_fit(f, zmock_raw, beta_mock_raw, sigma=beta_err_mock_raw)
print("beta(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_beta_mock_raw[0][0], p_beta_mock_raw[1][0,0], p_beta_mock_raw[0][1], p_beta_mock_raw[1][1,1]))
p_bias_eta_mock_raw = sp.optimize.curve_fit(f, zmock_raw, bias_eta_mock_raw, sigma=bias_eta_err_mock_raw)
print("bias_eta(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_bias_eta_mock_raw[0][0], p_bias_eta_mock_raw[1][0,0], p_bias_eta_mock_raw[0][1], p_bias_eta_mock_raw[1][1,1]))
p_beff_mock_raw = sp.optimize.curve_fit(f, zmock_raw, beff_mock_raw, sigma=beff_err_mock_raw)
print("beff(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_beff_mock_raw[0][0], p_beff_mock_raw[1][0,0], p_beff_mock_raw[0][1], p_beff_mock_raw[1][1,1]))


# les donnees
data_bias_eta= np.array( [[ 2.13781251, -0.18558474,  0.00650885],
                          [ 2.27259149, -0.20089904,  0.00566513],
                          [ 2.54744584, -0.23754895,  0.00859699],
                          [ 2.9263245,  -0.28984325,  0.01907289]] )
data_beta = np.array([[2.13781251, 1.71518188, 0.14670333],
                      [2.27259149, 1.63599596, 0.12377548],
                      [2.54744584, 1.34603853, 0.10685724],
                      [2.9263245,  1.20116875, 0.16525412]])
zdata = data_beta[:,0]
bias_eta_data = data_bias_eta[:,1]
bias_eta_err_data = data_bias_eta[:,2]
beta_data = data_beta[:,1]
beta_err_data = data_beta[:,2]
bias_data = bias_eta_data * 0.97 / beta_data
cor_data = -0.9
bias_err_data = util.bias_err(bias_eta_data, bias_eta_err_data, beta_data, beta_err_data, cor_data)
beff_data = bias_data * np.sqrt(1+2/3*beta_data+1/5*beta_data**2)
if correct_beff:
    beff_data *= (1+z0) / (1+zdata)
beff_err_data = bias_err_data

print("Fits on data:")
p_bias_data = sp.optimize.curve_fit(f, zdata, bias_data, sigma=bias_err_data)
print("bias(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_bias_data[0][0], p_bias_data[1][0,0], p_bias_data[0][1], p_bias_data[1][1,1]))
p_beta_data = sp.optimize.curve_fit(f, zdata, beta_data, sigma=beta_err_data)
print("beta(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_beta_data[0][0], p_beta_data[1][0,0], p_beta_data[0][1], p_beta_data[1][1,1]))
p_bias_eta_data = sp.optimize.curve_fit(f, zdata, bias_eta_data, sigma=bias_eta_err_data)
print("bias_eta(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_bias_eta_data[0][0], p_bias_eta_data[1][0,0], p_bias_eta_data[0][1], p_bias_eta_data[1][1,1]))
p_beff_data = sp.optimize.curve_fit(f, zdata, beff_data, sigma=beff_err_data)
print("beff(z) = ({:.4} +/- {:.4})*(1+z)^[{:.4} +/- {:.4}]".format(p_beff_data[0][0], p_beff_data[1][0,0], p_beff_data[0][1], p_beff_data[1][1,1]))


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
if plot_bias:
    fig, ax = plt.subplots()
    if plot_data:
        ax.errorbar(zdata, bias_data, yerr=bias_err_data, fmt='.', label='data', color='royalblue')
        z = np.linspace(zdata.min(), zdata.max(), 100)
        plt.plot(z, f(z, p_bias_data[0][0], p_bias_data[0][1]), linestyle='--', color='royalblue')
        if plot_shades:
            plt.fill_between(z, f(z, p_bias_data[0][0]+p_bias_data[1][0,0], p_bias_data[0][1]-p_bias_data[1][1,1]),
                         f(z, p_bias_data[0][0]-p_bias_data[1][0,0], p_bias_data[0][1]+p_bias_data[1][1,1]),
                         color='royalblue', alpha=0.2)
    if plot_mock1:
        ax.errorbar(zmock1, bias_mock1, yerr=bias_err_mock1, fmt='.', label='mock1', color='darkorange')
        z = np.linspace(zmock1.min(), zmock1.max(), 100)
        plt.plot(z, f(z, p_bias_mock1[0][0], p_bias_mock1[0][1]), linestyle='--', color='darkorange')
        if plot_shades:
            plt.fill_between(z, f(z, p_bias_mock1[0][0]+p_bias_mock1[1][0,0], p_bias_mock1[0][1]-p_bias_mock1[1][1,1]),
                         f(z, p_bias_mock1[0][0]-p_bias_mock1[1][0,0], p_bias_mock1[0][1]+p_bias_mock1[1][1,1]),
                         color='darkorange', alpha=0.2)
    if plot_mock2:
        ax.errorbar(zmock2, bias_mock2, yerr=bias_err_mock2, fmt='.', label='mock2', color='red')
        z = np.linspace(zmock2.min(), zmock2.max(), 100)
        plt.plot(z, f(z, p_bias_mock2[0][0], p_bias_mock2[0][1]), linestyle='--', color='red')
        if plot_shades:
            plt.fill_between(z, f(z, p_bias_mock2[0][0]+p_bias_mock2[1][0,0], p_bias_mock2[0][1]-p_bias_mock2[1][1,1]),
                         f(z, p_bias_mock2[0][0]-p_bias_mock2[1][0,0], p_bias_mock2[0][1]+p_bias_mock2[1][1,1]),
                         color='red', alpha=0.2)
    if plot_mock_raw:
        ax.errorbar(zmock_raw, bias_mock_raw, yerr=bias_err_mock_raw, fmt='.', label='mock_raw', color='green')
        z = np.linspace(zmock_raw.min(), zmock_raw.max(), 100)
        plt.plot(z, f(z, p_bias_mock_raw[0][0], p_bias_mock_raw[0][1]), linestyle='--', color='green')
        if plot_shades:
            plt.fill_between(z, f(z, p_bias_mock_raw[0][0]+p_bias_mock_raw[1][0,0], p_bias_mock_raw[0][1]-p_bias_mock_raw[1][1,1]),
                         f(z, p_bias_mock_raw[0][0]-p_bias_mock_raw[1][0,0], p_bias_mock_raw[0][1]+p_bias_mock_raw[1][1,1]),
                         color='green', alpha=0.2)
    if plot_pred1:
        ax.plot(zpred1, bias_pred1, 'x', color='darkgreen', label='pred 10<r<180')
    if plot_pred2:
        ax.plot(zpred2, bias_pred2, 'x', color='limegreen', label='pred 40<r<180')

    ax.legend()
    ax.grid()
    ax.set_xlabel('z')
    ax.set_ylabel('bias')

if plot_beta:
    fig, ax = plt.subplots()
    if plot_data:
        ax.errorbar(zdata, beta_data, yerr=beta_err_data, fmt='.', label='data', color='royalblue')
        z = np.linspace(zdata.min(), zdata.max(), 100)
        plt.plot(z, f(z, p_beta_data[0][0], p_beta_data[0][1]), linestyle='--', color='royalblue')
        if plot_shades:
            plt.fill_between(z, f(z, p_beta_data[0][0]+p_beta_data[1][0,0], p_beta_data[0][1]-p_beta_data[1][1,1]),
                         f(z, p_beta_data[0][0]-p_beta_data[1][0,0], p_beta_data[0][1]+p_beta_data[1][1,1]),
                         color='royalblue', alpha=0.2)
    if plot_mock1:
        ax.errorbar(zmock1, beta_mock1, yerr=beta_err_mock1, fmt='.', label='mock1', color='darkorange')
        z = np.linspace(zmock1.min(), zmock1.max(), 100)
        plt.plot(z, f(z, p_beta_mock1[0][0], p_beta_mock1[0][1]), linestyle='--', color='darkorange')
        if plot_shades:
            plt.fill_between(z, f(z, p_beta_mock1[0][0]+p_beta_mock1[1][0,0], p_beta_mock1[0][1]-p_beta_mock1[1][1,1]),
                         f(z, p_beta_mock1[0][0]-p_beta_mock1[1][0,0], p_beta_mock1[0][1]+p_beta_mock1[1][1,1]),
                         color='darkorange', alpha=0.2)
    if plot_mock2:
        ax.errorbar(zmock2, beta_mock2, yerr=beta_err_mock2, fmt='.', label='mock2', color='red')
        z = np.linspace(zmock2.min(), zmock2.max(), 100)
        plt.plot(z, f(z, p_beta_mock2[0][0], p_beta_mock2[0][1]), linestyle='--', color='red')
        if plot_shades:
            plt.fill_between(z, f(z, p_beta_mock2[0][0]+p_beta_mock2[1][0,0], p_beta_mock2[0][1]-p_beta_mock2[1][1,1]),
                         f(z, p_beta_mock2[0][0]-p_beta_mock2[1][0,0], p_beta_mock2[0][1]+p_beta_mock2[1][1,1]),
                         color='red', alpha=0.2)
    if plot_mock_raw:
        ax.errorbar(zmock_raw, beta_mock_raw, yerr=beta_err_mock_raw, fmt='.', label='mock_raw', color='green')
        z = np.linspace(zmock_raw.min(), zmock_raw.max(), 100)
        plt.plot(z, f(z, p_beta_mock_raw[0][0], p_beta_mock_raw[0][1]), linestyle='--', color='green')
        if plot_shades:
            plt.fill_between(z, f(z, p_beta_mock_raw[0][0]+p_beta_mock_raw[1][0,0], p_beta_mock_raw[0][1]-p_beta_mock_raw[1][1,1]),
                         f(z, p_beta_mock_raw[0][0]-p_beta_mock_raw[1][0,0], p_beta_mock_raw[0][1]+p_beta_mock_raw[1][1,1]),
                         color='green', alpha=0.2)
    if plot_pred1:
        ax.plot(zpred1, beta_pred1, 'x', color='darkgreen', label='pred 10<r<180')
    if plot_pred2:
        ax.plot(zpred2, beta_pred2, 'x', color='limegreen', label='pred 40<r<180')

    ax.legend()
    ax.grid()
    ax.set_xlabel('z')
    ax.set_ylabel('beta')

if plot_bias_eta:
    fig, ax = plt.subplots()
    if plot_data:
        ax.errorbar(zdata, bias_eta_data, yerr=bias_eta_err_data, fmt='.', label='data', color='royalblue')
        z = np.linspace(zdata.min(), zdata.max(), 100)
        plt.plot(z, f(z, p_bias_eta_data[0][0], p_bias_eta_data[0][1]), linestyle='--', color='royalblue')
        if plot_shades:
            plt.fill_between(z, f(z, p_bias_eta_data[0][0]+p_bias_eta_data[1][0,0], p_bias_eta_data[0][1]-p_bias_eta_data[1][1,1]),
                         f(z, p_bias_eta_data[0][0]-p_bias_eta_data[1][0,0], p_bias_eta_data[0][1]+p_bias_eta_data[1][1,1]),
                         color='royalblue', alpha=0.2)
    if plot_mock1:
        ax.errorbar(zmock1, bias_eta_mock1, yerr=bias_eta_err_mock1, fmt='.', label='mock1', color='darkorange')
        z = np.linspace(zmock1.min(), zmock1.max(), 100)
        plt.plot(z, f(z, p_bias_eta_mock1[0][0], p_bias_eta_mock1[0][1]), linestyle='--', color='darkorange')
        if plot_shades:
            plt.fill_between(z, f(z, p_bias_eta_mock1[0][0]+p_bias_eta_mock1[1][0,0], p_bias_eta_mock1[0][1]-p_bias_eta_mock1[1][1,1]),
                         f(z, p_bias_eta_mock1[0][0]-p_bias_eta_mock1[1][0,0], p_bias_eta_mock1[0][1]+p_bias_eta_mock1[1][1,1]),
                         color='darkorange', alpha=0.2)
    if plot_mock2:
        ax.errorbar(zmock2, bias_eta_mock2, yerr=bias_eta_err_mock2, fmt='.', label='mock2', color='red')
        z = np.linspace(zmock2.min(), zmock2.max(), 100)
        plt.plot(z, f(z, p_bias_eta_mock2[0][0], p_bias_eta_mock2[0][1]), linestyle='--', color='red')
        if plot_shades:
            plt.fill_between(z, f(z, p_bias_eta_mock2[0][0]+p_bias_eta_mock2[1][0,0], p_bias_eta_mock2[0][1]-p_bias_eta_mock2[1][1,1]),
                         f(z, p_bias_eta_mock2[0][0]-p_bias_eta_mock2[1][0,0], p_bias_eta_mock2[0][1]+p_bias_eta_mock2[1][1,1]),
                         color='red', alpha=0.2)
    if plot_mock_raw:
        ax.errorbar(zmock_raw, bias_eta_mock_raw, yerr=bias_eta_err_mock_raw, fmt='.', label='mock_raw', color='green')
        z = np.linspace(zmock_raw.min(), zmock_raw.max(), 100)
        plt.plot(z, f(z, p_bias_eta_mock_raw[0][0], p_bias_eta_mock_raw[0][1]), linestyle='--', color='green')
        if plot_shades:
            plt.fill_between(z, f(z, p_bias_eta_mock_raw[0][0]+p_bias_eta_mock_raw[1][0,0], p_bias_eta_mock_raw[0][1]-p_bias_eta_mock_raw[1][1,1]),
                         f(z, p_bias_eta_mock_raw[0][0]-p_bias_eta_mock_raw[1][0,0], p_bias_eta_mock_raw[0][1]+p_bias_eta_mock_raw[1][1,1]),
                         color='green', alpha=0.2)
    if plot_pred1:
        ax.plot(zpred1, bias_eta_pred1, 'x', color='darkgreen', label='pred 10<r<180')
    if plot_pred2:
        ax.plot(zpred2, bias_eta_pred2, 'x', color='limegreen', label='pred 40<r<180')

    ax.legend()
    ax.grid()
    ax.set_xlabel('z')
    ax.set_ylabel('bias_eta')

if plot_beff:
    fig, ax = plt.subplots()
    if plot_data:
        ax.errorbar(zdata, beff_data, yerr=beff_err_data, fmt='.', label='data', color='royalblue')
        z = np.linspace(zdata.min(), zdata.max(), 100)
        plt.plot(z, f(z, p_beff_data[0][0], p_beff_data[0][1]), linestyle='--', color='royalblue')
        if plot_shades:
            plt.fill_between(z, f(z, p_beff_data[0][0]+p_beff_data[1][0,0], p_beff_data[0][1]-p_beff_data[1][1,1]),
                         f(z, p_beff_data[0][0]-p_beff_data[1][0,0], p_beff_data[0][1]+p_beff_data[1][1,1]),
                         color='royalblue', alpha=0.2)
    if plot_mock1:
        ax.errorbar(zmock1, beff_mock1, yerr=beff_err_mock1, fmt='.', label='mock1', color='darkorange')
        z = np.linspace(zmock1.min(), zmock1.max(), 100)
        plt.plot(z, f(z, p_beff_mock1[0][0], p_beff_mock1[0][1]), linestyle='--', color='darkorange')
        if plot_shades:
            plt.fill_between(z, f(z, p_beff_mock1[0][0]+p_beff_mock1[1][0,0], p_beff_mock1[0][1]-p_beff_mock1[1][1,1]),
                         f(z, p_beff_mock1[0][0]-p_beff_mock1[1][0,0], p_beff_mock1[0][1]+p_beff_mock1[1][1,1]),
                         color='darkorange', alpha=0.2)
    if plot_mock2:
        ax.errorbar(zmock2, beff_mock2, yerr=beff_err_mock2, fmt='.', label='mock2', color='red')
        z = np.linspace(zmock2.min(), zmock2.max(), 100)
        plt.plot(z, f(z, p_beff_mock2[0][0], p_beff_mock2[0][1]), linestyle='--', color='red')
        if plot_shades:
            plt.fill_between(z, f(z, p_beff_mock2[0][0]+p_beff_mock2[1][0,0], p_beff_mock2[0][1]-p_beff_mock2[1][1,1]),
                         f(z, p_beff_mock2[0][0]-p_beff_mock2[1][0,0], p_beff_mock2[0][1]+p_beff_mock2[1][1,1]),
                         color='red', alpha=0.2)
    if plot_mock_raw:
        ax.errorbar(zmock_raw, beff_mock_raw, yerr=beff_err_mock_raw, fmt='.', label='mock_raw', color='green')
        z = np.linspace(zmock_raw.min(), zmock_raw.max(), 100)
        plt.plot(z, f(z, p_beff_mock_raw[0][0], p_beff_mock_raw[0][1]), linestyle='--', color='green')
        if plot_shades:
            plt.fill_between(z, f(z, p_beff_mock_raw[0][0]+p_beff_mock_raw[1][0,0], p_beff_mock_raw[0][1]-p_beff_mock_raw[1][1,1]),
                         f(z, p_beff_mock_raw[0][0]-p_beff_mock_raw[1][0,0], p_beff_mock_raw[0][1]+p_beff_mock_raw[1][1,1]),
                         color='green', alpha=0.2)
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
