# Dielectric function
# used as external function
# Last modified: Dec 11 2017
# Coded by Takuro TOKUNAGA

import math
import numpy as np
import cmath
import time
from scipy.integrate import trapz, simps, quad, quadrature, romberg
start = time.time()

# Constants
ipureal = 20
ipuimag = 0.0001

# Dielectric function
def problem2(omega):
    function = ipureal + ipuimag*1j

    return function

# time display
#elapsed_time = time.time()-start
#print("dc_problem2 elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
