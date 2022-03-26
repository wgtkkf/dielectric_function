# Dielectric function of Au
# used as external function
# J. Chem. Phys. 127 189901 (2007)
# Last modified: June 05 2018
# Coded by Takuro TOKUNAGA

import math
import numpy as np
import cmath
import time
from scipy.integrate import trapz, simps, quad, quadrature, romberg
start = time.time()

# unit conversion
ucnm = 1.0*np.power(10.,-9)

# Constants
sc = 299792458 # small c, speed of light
#
ipu_inf = 1.54
lambda_p = 143*ucnm # [m]
omega_p = (2*np.pi*sc)/lambda_p
gamma_p = 14500*ucnm # [m]
lgamma_p = (2*np.pi*sc)/gamma_p

# i=1
la1 = 1.27 # large A1
phi1 = -np.pi/4 # [rad]
lambda1 = 470*ucnm # [m]
omega1 = (2*np.pi*sc)/lambda1
gamma1 = 1900*ucnm # [m]
lgamma1 = (2*np.pi*sc)/gamma1 # large gamma1
mu1 = -1

# i=2
la2 = 1.10 # large A2
phi2 = -np.pi/4 # [rad]
lambda2 = 325*ucnm # [m]
omega2 = (2*np.pi*sc)/lambda2
gamma2 = 1060*ucnm # [m]
lgamma2 = (2*np.pi*sc)/gamma2 # large gamma2
mu2 = -1

# parameters for drude model
omegap_drude = 13.71*np.power(10.,15) # [rad/s]
gamma_drude = 4.05*np.power(10.,13) # [s-1]

# Dielectric function
def au(omega):
    term1 = ipu_inf - np.power(omega_p,2)/(np.power(omega,2)+1j*lgamma_p*omega)
    term2 = (la1*omega)*(np.exp(1j*phi1)*np.power(omega1-omega-1j*lgamma1,mu1)+\
    np.exp(-1j*phi1)*np.power(omega1+omega+1j*lgamma1,mu1))
    term3 = (la2*omega)*(np.exp(1j*phi2)*np.power(omega2-omega-1j*lgamma2,mu2)+\
    np.exp(-1j*phi2)*np.power(omega2+omega+1j*lgamma2,mu2))

    function = term1 + term2 + term3

    return function

def au_drude(omega):
    omegap2 = np.power(omegap_drude,2)
    gammap2 = np.power(gamma_drude,2)

    real = 1-(omegap2/(np.power(omega,2)+gammap2))
    imag = omegap2*gamma_drude/(omega*(np.power(omega,2)+gammap2))

    function = real + 1j*imag

    return function

# time display
#elapsed_time = time.time()-start
#print("dc_au elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
