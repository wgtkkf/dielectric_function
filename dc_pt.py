# Dielectric function of Platinum
# A. B. Djurisic et al, Appl. Opt. 28 7097 (1997)
# Frequency range
# Last modified: April 24 2018
# Coded by Takuro TOKUNAGA
# working properly

import math
import numpy as np
import cmath
import time
from scipy.integrate import trapz, simps, quad, quadrature, romberg
start = time.time()

# parameters
rh = 6.626070040*np.power(10.,-34)/(2*np.pi) # reduced Planck constant
ev = 1.602176620898*np.power(10.,-19) # 1 ev to joules

num = 7 # from 0 to 6, dont change. based on the journal.
sf = np.zeros(num,dtype=np.float64) # small f
gamma_t = np.zeros(num,dtype=np.float64) # _t: table
omega_t = np.zeros(num,dtype=np.float64) # _t: table

# Initialization, table 5 in the journal
sf[0] = 1.090
sf[1] = 1.490
sf[2] = 4.631
sf[3] = 0.026
sf[4] = 3.091
sf[5] = 6.909
sf[6] = 8.969
#
gamma_t[0] = 0.078
gamma_t[1] = 0.838
gamma_t[2] = 5.194
gamma_t[3] = 0.403
gamma_t[4] = 5.467
gamma_t[5] = 11.539
gamma_t[6] = 9.274
#
omega_t[0] = 0
omega_t[1] = 0.861
omega_t[2] = 2.464
omega_t[3] = 6.168
omega_t[4] = 9.293
omega_t[5] = 14.264
omega_t[6] = 20.012
#
omegap = 5.15
omegap2 = np.power(omegap,2)
ohmp = np.sqrt(sf[0])*omegap
ohmp2 = np.power(ohmp,2)
gamma0 = gamma_t[0]

# Dielectric function
def pt(omega):
    # initialize
    ipuf = 0
    ipub = 0
    function = 0

    # if omega = 0

    # omega conversion:
    omega = (rh/ev)*omega

    ipuf = 1-ohmp2/(omega*(omega+1j*gamma0))

    for i in range(1, num): # 1 to 6
        ipub = ipub - (sf[i]*omegap2)/(np.power(omega,2)\
        -np.power(omega_t[i],2)+1j*omega*gamma_t[i])

    function = ipuf+ipub

    return function

# time display
#elapsed_time = time.time()-start
#print("dc_pt elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
