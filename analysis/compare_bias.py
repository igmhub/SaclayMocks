import numpy as np
from SaclayMocks import util

zref = 2.3

pars0 = util.extract_h5file("/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-raw/cf_z_0_10-exp.h5")
pars1 = util.extract_h5file("/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.0/cf_z_0_10-exp.h5")
pars2 = util.extract_h5file("/global/cfs/cdirs/desi/science/lya/picca_on_mocks/saclay/v4.7/global-01_30/eboss-0.2/cf_z_0_10-exp.h5")
pars3 = util.extract_h5file("/global/project/projectdirs/eboss/lya_forest/dr16/redo_4_zbins/Fits/cf/Kaiser_sky_met_hcd/z_0_10/result.h5")

bias_gamma0 = 3.46
bias_gamma1 = 3.256
bias_gamma2 = 3.183
bias_gamma3 = 3.473

beta_gamma0 = -2.497
beta_gamma1 = -2.371
beta_gamma2 = -2.116
beta_gamma3 = -2.316

bias0 = pars0[2]['beff_LYA']*((1+zref)/(1+pars0[2]['zeff']))**bias_gamma0
bias1 = pars1[2]['beff_LYA']*((1+zref)/(1+pars1[2]['zeff']))**bias_gamma1
bias2 = pars2[2]['beff_LYA']*((1+zref)/(1+pars2[2]['zeff']))**bias_gamma2
bias3 = pars3[2]['beff_LYA']*((1+zref)/(1+pars3[2]['zeff']))**bias_gamma3
bias_err0 = pars0[3]['beff_LYA']*((1+zref)/(1+pars0[2]['zeff']))**bias_gamma0
bias_err1 = pars1[3]['beff_LYA']*((1+zref)/(1+pars1[2]['zeff']))**bias_gamma1
bias_err2 = pars2[3]['beff_LYA']*((1+zref)/(1+pars2[2]['zeff']))**bias_gamma2
bias_err3 = pars3[3]['beff_LYA']*((1+zref)/(1+pars3[2]['zeff']))**bias_gamma3

beta0 = pars0[2]['beta_LYA']*((1+zref)/(1+pars0[2]['zeff']))**beta_gamma0
beta1 = pars1[2]['beta_LYA']*((1+zref)/(1+pars1[2]['zeff']))**beta_gamma1
beta2 = pars2[2]['beta_LYA']*((1+zref)/(1+pars2[2]['zeff']))**beta_gamma2
beta3 = pars3[2]['beta_LYA']*((1+zref)/(1+pars3[2]['zeff']))**beta_gamma3
beta_err0 = pars0[3]['beta_LYA']*((1+zref)/(1+pars0[2]['zeff']))**beta_gamma0
beta_err1 = pars1[3]['beta_LYA']*((1+zref)/(1+pars1[2]['zeff']))**beta_gamma1
beta_err2 = pars2[3]['beta_LYA']*((1+zref)/(1+pars2[2]['zeff']))**beta_gamma2
beta_err3 = pars3[3]['beta_LYA']*((1+zref)/(1+pars3[2]['zeff']))**beta_gamma3

# # horizontal
# print("\\toprule")
# print("Param\\`etre  & raw mocks & mocks eboss-0.0 & mocks eboss-0.2 & donn\\'ees DR16 \\\\")
# print("\\midrule")
# row = "$b_{\\mathrm{eff},\\mathrm{Ly}\\alpha}$"
# row += " & ${} \\pm {}$".format(np.round(bias0,4), np.round(bias_err0,4))
# row += " & ${} \\pm {}$".format(np.round(bias1,4), np.round(bias_err1,4))
# row += " & ${} \\pm {}$".format(np.round(bias2,4), np.round(bias_err2,4))
# row += " & ${} \\pm {}$".format(np.round(bias3,4), np.round(bias_err3,4))
# row += " \\\\"
# print(row)
# row = "$\\beta_{\\mathrm{Ly}\\alpha}$"
# row += " & ${} \\pm {}$".format(np.round(beta0,3), np.round(beta_err0,3))
# row += " & ${} \\pm {}$".format(np.round(beta1,3), np.round(beta_err1,3))
# row += " & ${} \\pm {}$".format(np.round(beta2,3), np.round(beta_err2,3))
# row += " & ${} \\pm {}$".format(np.round(beta3,3), np.round(beta_err3,3))
# row += " \\\\"
# print(row)
# print("\\bottomrule")

# vertical
print("\\toprule")
print(" & $b_{\\mathrm{eff},\\mathrm{Ly}\\alpha}$ & $\\beta_{\\mathrm{Ly}\\alpha}$ \\\\")
print("\\midrule")
row = "donn\\'ees DR16"
row += " & ${} \\pm {}$".format(np.round(bias3,4), np.round(bias_err3,4))
row += " & ${} \\pm {}$".format(np.round(beta3,3), np.round(beta_err3,3))
row += " \\\\"
print(row)
row = "raw mocks"
row += " & ${} \\pm {}$".format(np.round(bias0,4), np.round(bias_err0,4))
row += " & ${} \\pm {}$".format(np.round(beta0,3), np.round(beta_err0,3))
row += " \\\\"
print(row)
row = "mocks eboss-0.0"
row += " & ${} \\pm {}$".format(np.round(bias1,4), np.round(bias_err1,4))
row += " & ${} \\pm {}$".format(np.round(beta1,3), np.round(beta_err1,3))
row += " \\\\"
print(row)
row = "mocks eboss-0.2"
row += " & ${} \\pm {}$".format(np.round(bias2,4), np.round(bias_err2,4))
row += " & ${} \\pm {}$".format(np.round(beta2,3), np.round(beta_err2,3))
row += " \\\\"
print(row)
print("\\bottomrule")
