

r_pow = 2
rtmax = 200
rpmax = 200

f, ax = plt.subplots(figsize=(12,8))

data = f1
da = data[1]['DA'][:]
co = data[1]['CO'][:]
mu0, mu1, mu2, mu3 = 0, 0.5, 0.8, 1
w1 = picca.wedgize.wedge(mumin=mu0,mumax=mu1, rtmax=rtmax, rpmax=rpmax)
w2 = picca.wedgize.wedge(mumin=mu1,mumax=mu2, rtmax=rtmax, rpmax=rpmax)
w3 = picca.wedgize.wedge(mumin=mu2,mumax=mu3, rtmax=rtmax, rpmax=rpmax)
data_wedge1 = w1.wedge(da,co)
coef1 = data_wedge1[0]**r_pow
data_wedge2 = w2.wedge(da,co)
coef2 = data_wedge2[0]**r_pow
data_wedge3 = w3.wedge(da,co)
coef3 = data_wedge3[0]**r_pow

ax.errorbar(data_wedge1[0],coef1*data_wedge1[1],yerr=coef1*sp.sqrt(sp.diag(data_wedge1[2])),marker='.', label="1", color='r')
ax.errorbar(data_wedge2[0],coef2*data_wedge2[1],yerr=coef2*sp.sqrt(sp.diag(data_wedge2[2])),marker='.', color='r')
ax.errorbar(data_wedge3[0],coef3*data_wedge3[1],yerr=coef3*sp.sqrt(sp.diag(data_wedge3[2])),marker='.', color='r')

data = f2
da = data[1]['DA'][:]
co = data[1]['CO'][:]
mu0, mu1, mu2, mu3 = 0, 0.5, 0.8, 1
w1 = picca.wedgize.wedge(mumin=mu0,mumax=mu1, rtmax=rtmax, rpmax=rpmax)
w2 = picca.wedgize.wedge(mumin=mu1,mumax=mu2, rtmax=rtmax, rpmax=rpmax)
w3 = picca.wedgize.wedge(mumin=mu2,mumax=mu3, rtmax=rtmax, rpmax=rpmax)
data_wedge1 = w1.wedge(da,co)
coef1 = data_wedge1[0]**r_pow
data_wedge2 = w2.wedge(da,co)
coef2 = data_wedge2[0]**r_pow
data_wedge3 = w3.wedge(da,co)
coef3 = data_wedge3[0]**r_pow

ax.errorbar(data_wedge1[0],coef1*data_wedge1[1],yerr=coef1*sp.sqrt(sp.diag(data_wedge1[2])),marker='.', label="2", color='b')
ax.errorbar(data_wedge2[0],coef2*data_wedge2[1],yerr=coef2*sp.sqrt(sp.diag(data_wedge2[2])),marker='.', color='b')
ax.errorbar(data_wedge3[0],coef3*data_wedge3[1],yerr=coef3*sp.sqrt(sp.diag(data_wedge3[2])),marker='.', color='b')

data = f3
da = data[1]['DA'][:]
co = data[1]['CO'][:]
mu0, mu1, mu2, mu3 = 0, 0.5, 0.8, 1
w1 = picca.wedgize.wedge(mumin=mu0,mumax=mu1, rtmax=rtmax, rpmax=rpmax)
w2 = picca.wedgize.wedge(mumin=mu1,mumax=mu2, rtmax=rtmax, rpmax=rpmax)
w3 = picca.wedgize.wedge(mumin=mu2,mumax=mu3, rtmax=rtmax, rpmax=rpmax)
data_wedge1 = w1.wedge(da,co)
coef1 = data_wedge1[0]**r_pow
data_wedge2 = w2.wedge(da,co)
coef2 = data_wedge2[0]**r_pow
data_wedge3 = w3.wedge(da,co)
coef3 = data_wedge3[0]**r_pow

ax.errorbar(data_wedge1[0],coef1*data_wedge1[1],yerr=coef1*sp.sqrt(sp.diag(data_wedge1[2])),marker='.', label="3", color='g')
ax.errorbar(data_wedge2[0],coef2*data_wedge2[1],yerr=coef2*sp.sqrt(sp.diag(data_wedge2[2])),marker='.', color='g')
ax.errorbar(data_wedge3[0],coef3*data_wedge3[1],yerr=coef3*sp.sqrt(sp.diag(data_wedge3[2])),marker='.', color='g')

ax.set_xlabel(r"$r \, [\mathrm{Mpc \, h^{-1}}]$",fontsize=20)
if r_pow == 2:
    ax.set_ylabel(r"$r^{2}\xi(r) \, [\mathrm{Mpc \, h^{-1}}]$",fontsize=20)
if r_pow == 1:
    ax.set_ylabel(r"$r\xi(r) \, [\mathrm{Mpc \, h^{-1}}]$",fontsize=20)
if r_pow == 0:
    ax.set_ylabel(r"$\xi(r) \, [\mathrm{Mpc \, h^{-1}}]$",fontsize=20)
ax.grid()
ax.legend()
plt.show()

