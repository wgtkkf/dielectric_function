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
ipureal = 1
ipuimag = 0.0001

# Dielectric function
def problem1(omega):
    function = ipureal + ipuimag*1j

    return function

# time display
#elapsed_time = time.time()-start
#print("dc_problem1 elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
