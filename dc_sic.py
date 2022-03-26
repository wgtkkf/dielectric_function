# Dielectric function of SiC
# used as external function
# J. Phys. D: Appl. Phys. 43 075501 (2010)
# Last modified: August 29 2017
# Coded by Takuro TOKUNAGA

import math
import numpy as np
import cmath
import time
start = time.time()

# Constants
ipuinf = 6.7
omegalo = 1.825*np.power(10.,14)
omegato = 1.494*np.power(10.,14)
gamma = 8.966*np.power(10.,11)

# Dielectric function
def sic(omega):
    numerator = np.power(omega,2)-np.power(omegalo,2)+1j*gamma*omega
    denominator = np.power(omega,2)-np.power(omegato,2)+1j*gamma*omega
    function = ipuinf*(numerator/denominator)

    return function

# time display
#elapsed_time = time.time()-start
#print("dc_sic elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
