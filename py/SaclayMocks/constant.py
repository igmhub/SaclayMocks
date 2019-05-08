# physical constant and parameters
import scipy as sp
from scipy import interpolate
import scipy.constants
from pkg_resources import resource_filename

c = scipy.constants.speed_of_light / 1000.  # in km/s

n_qso_exp = 100.  # expected number of qso per square degrees

lya = 1215.67 ## angstrom https://en.wikipedia.org/wiki/Lyman-alpha_line  1215.668 and 1215.674
lylimit = lya * 3 /4. # 911.75  ok with https://en.wikipedia.org/wiki/Hydrogen_spectral_series
lyb = lylimit * 9./8. # 1025.72  https://en.wikipedia.org/wiki/Hydrogen_spectral_series : 1025.7
lambda_min = 3530.  # Lambda min in A (given by Stephen Bailey - 04/05/2018)

# Not used:
R_H = 10973731.568 # m^-1  https://en.wikipedia.org/wiki/Rydberg_constant
                # pdg gives hc R_H -> 100000*13.605693009/197.3269788/2/np.pi
                # i.e. 1097.373156849444
#lylimit = 1E10/R_H # 911.27 angtroem !!!
#lya = lylimit * 4 / 3.  #  1215.0227
#lyb = lylimit * 9 / 8.  # 1025.18


# Cosmo params given by Helion in a mail (30/03/2018)
# should correspond to Planck 2015
h = 0.6731
omega_M_0 = 0.31457
omega_b_0 = 0.045  # only use in run_camb.py
omega_lambda_0 = 0.68543
omega_k_0 = 0.0
ns = 0.96  # only used in run_camb.py

QSO_bias = 3.7  # should become b(z) in data.py / QSO.py
z_QSO_bias = 2.4  # QSO_bias taken from Helion thesis
z0 = 1.70975268202
# box_pixel = 1600./1024.   # rather be stored in box.fits

H0 = 100.  # In km/s /(Mpc/h)

deg2rad = sp.pi/180.
rad2deg = 180/sp.pi

#boss_lambda_min = 3600.


### Define path to files
path_PlanckDR12 = resource_filename('SaclayMocks', '/etc/PlanckDR12.fits')


'''
class cosmo:

    def __init__(self,Om,Ok=0):
        H0 = 100. ## km/s/Mpc
        ## ignore Orad and neutrinos
        nbins=10000
        zmax=10.
        dz = zmax/nbins
        z=sp.array(range(nbins))*dz
        hubble = H0*sp.sqrt(Om*(1+z)**3+Ok*(1+z)**2+1-Ok-Om)

        chi=sp.zeros(nbins)
        for i in range(1,nbins):
            chi[i]=chi[i-1]+c*(1./hubble[i-1]+1/hubble[i])/2.*dz

        self.r_comoving = interpolate.interp1d(z,chi)

        ## da here is the comoving angular diameter distance
        da = chi
        if Ok<0:
            da = sp.sin(H0*sp.sqrt(-Ok)/c*chi)/(H0*sp.sqrt(-Ok)/c)
        if Ok>0:
            da = sp.sinh(H0*sp.sqrt(Ok)/c*chi)/(H0*sp.sqrt(Ok)/c)
        self.da = interpolate.interp1d(z,da)
        self.hubble = interpolate.interp1d(z,hubble)
        self.r_2_z = interpolate.interp1d(chi,z)
'''
