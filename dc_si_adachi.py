# Dielectric function of Si
# used as external function
# Last modified: August 30 2017
# Coded by Takuro TOKUNAGA

import math
import numpy as np
import cmath
import time
from scipy.integrate import trapz, simps, quad, quadrature, romberg
start = time.time()

# Constants
h = 6.626070040*math.pow(10,-34)/(2*np.pi) # reduced planck constant
c = 299792458 # light speed
kb = 1.38064852*math.pow(10,-23) # boltzmann constant
ev = 1.602176620898*math.pow(10,-19)

# parameters for dielectric function of intrinsic Silicon
# Constants for eq4, J.Appl.Phys. 69 1574 (1991)
b1 = 5.22 # dimensionless TABLEIII
le1 = 3.38 # eV TABLEIII

# Constants for eq10, J.Appl.Phys. 69 1574 (1991)
lc = 3.01 # dimensionless TABLEIII
gamma = 0.127 # dimensionless TABLEIII

# Constants for eq15, J.Appl.Phys. 69 1574 (1991)
lf = 3.51 # dimensionless TABLEIII
lecl = le1 # Ecl = E1, eV TABLEIII

# Constants for eq10, eq15, J.Appl.Phys. 69 1574 (1991)
le2 = 4.27 # eV TABLEIII

# Constants for eq4, eq15, J.Appl.Phys. 69 1574 (1991)
lgamma2 = 3*math.pow(10,-2)

# Constants for eq2122, Phs.Rev.B 38 12966 (1988)
lca = 0.21 # dimensionless TABLEIII
lea = 5.32 # eV TABLEIII
gammaa = 0.089*2 # eq2122

# Additional term for real part, Phs.Rev.B 38 12966 (1988)
ipu1inf = 1.8 # dimensionless TABLEIII

# Integral part of prop
def si_adachi(omega):

    # Dielectric function of intrinsic Silicon
    # eq4, J.Appl.Phys. 69 1574 (1991)
    xi1d2 = np.power(((h/ev)*omega+1j*lgamma2)/le1,2)
    eq4 = (-b1/xi1d2)*np.log(1-xi1d2) # ignore second term

    # eq10, J.Appl.Phys. 69 1574 (1991)
    xi2 = ((h/ev)*omega)/le2
    xi22 = np.power(xi2,2)
    eq10 = lc/((1-xi22)-1j*xi2*gamma)

    # eq15, J.Appl.Phys. 69 1574 (1991)
    xicl2 = np.power(((h/ev)*omega+1j*lgamma2)/lecl,2)
    xi2g2 = np.power(((h/ev)*omega+1j*lgamma2)/le2,2)
    eq15 = -(lf/xi2g2)*np.log((1-xicl2)/(1-xi2g2))

    # eq 21 & 22, Phs.Rev.B 38 12966 (1988)
    xia = ((h/ev)*omega)/lea
    xia2 = np.power(xia,2)
    eq2122 = lca/((1-xia2)-1j*xia*gammaa)

    # ipu1inf
    ipu1inf = 1.8

    # dielectric function of film 1 & 3
    function = eq4 + eq10 + eq15 + eq2122 + ipu1inf

    return function

# time display
#elapsed_time = time.time()-start
#print("dc_si_adachi elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
