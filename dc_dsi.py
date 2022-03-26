# Dielectric function of doped silicon
# used as external function
# Phys. Rev. B 90 115433 (2014)
# Last modified: August 30 2017
# Coded by Takuro TOKUNAGA
# imperfecet

import math
import numpy as np
import cmath
import time
from scipy.integrate import trapz, simps, quad, quadrature, romberg
start = time.time()

# Parameters
epuinf = 6.7
largeN = 1.0*math.pow(10, 21)*math.pow(10, 6) # N, Doping level
se = 1.602176620898*math.pow(10,-19) # small e, electron elementary charge
mzero = 9.1093835611*math.pow(10, -31) # electron rest mass
mstar = 0.27*mzero # m*
ipuzero = 8.854187817*math.pow(10,-12) # vecuum permittivity
temperature = 300

# Dielectric function
def drude(omega): # equation (15)
    # Tn
    Tn = temperature/300

    # mue
    mue_numerator = 7.4*math.pow(10,8)*math.pow(temperature,-2.33)
    mue_denominator = 1+(0.88/1.26)*math.pow(10,-17)*largeN*math.pow(Tn,-2.546)
    mue = 88*math.pow(Tn,-0.57)+ mue_numerator/mue_denominator

    # gamma
    gamma = se/(mstar*mue)

    # omega-p
    omegap2 = largeN*np.power(se,2)/(mstar*ipuzero)

    # drude model, equation (14)
    numerator = omegap2
    denominator = np.power(omega,2)+1j*omega*gamma
    function = epuinf-(numerator/denominator)

    return function

# time display
#elapsed_time = time.time()-start
#print("dc_dsi elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
