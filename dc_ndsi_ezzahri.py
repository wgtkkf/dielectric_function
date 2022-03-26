# Dielectric function of p-doped silicon
# used as external function
# Coded by Takuro TOKUNAGA
# Last modified: January 18 2018

import math
import numpy as np
import cmath
import time
from scipy.integrate import trapz, simps, quad, quadrature
start = time.time()

# N-type ionization constants
dpc = 1.0*np.power(10.,21) # Doping concentration
ec = 1.60217657*np.power(10.,-19) # electron charge
mzero = 9.10938215*np.power(10.,-31) # Electron Rest mass
mstar = 0.27*mzero # Effective mass of an electron
epuzero = 8.85418782*np.power(10.,-12) # Permittivity of free space

temp = 300 # temperature
tn = temp/300 # reduced temperature
epub = 11.7

# plasma frequency
omegap2 = dpc*np.power(ec,2)/(mstar*epuzero)

# N-type Mobility
numerator = 7.4*np.power(10.,8)*np.power(tn,-2.33)
denominator = 1+(0.88/1.26)*np.power(10.,-17)*dpc*np.power(tn,-2.546)
mue =  88*np.power(tn,-0.57)+(numerator/denominator)

# scattering Rate
gamma = ec/(mstar*mue)

# dielectric function
def ndsi_ezzahri(omega):
    y = epub-np.power(omegap2,2)/(omega*(omega + 1j*gamma))
    return y

# display time
#elapsed_time = time.time()-start
#print("dc_ndsi elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
