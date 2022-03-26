# Appl. Opt. 51 6790 (2012)
# Compared with Palik's
# Coded by Takuro TOKUNAGA
# Last modified: March 18 2019

import math
import numpy as np
import cmath
import time
start = time.time()

# define parameters
number = 5 # number of j change here
sigma=np.zeros((number), dtype='float64')
nyu_zero=np.zeros((number), dtype='float64')
nyu_p=np.zeros((number), dtype='float64')
nyu_tau=np.zeros((number), dtype='float64')

# physical constant
c = 299792458 # [m/s]

# Initialization of parameters
# SO40-200, Table 2
ipu_inf = 2.09
factor = (2*c*np.pi)/np.power(10,-2.0) # cm-1 to rad/s

sigma[0] = 63*factor  # [cm-1] to [rad/s]
sigma[1] = 122*factor # [cm-1] to [rad/s]
sigma[2] = 67*factor  # [cm-1] to [rad/s]
sigma[3] = 316*factor # [cm-1] to [rad/s]
sigma[4] = 85*factor # [cm-1] to [rad/s]

nyu_zero[0] = 1046*factor # [cm-1] to [rad/s]
nyu_zero[1] = 1167*factor # [cm-1] to [rad/s]
nyu_zero[2] = 1058*factor # [cm-1] to [rad/s]
nyu_zero[3] = 434*factor  # [cm-1] to [rad/s]
nyu_zero[4] = 799*factor # [cm-1] to [rad/s]

nyu_p[0] = 544*factor # [cm-1] to [rad/s]
nyu_p[1] = 309*factor # [cm-1] to [rad/s]
nyu_p[2] = 466*factor # [cm-1] to [rad/s]
nyu_p[3] = 427*factor # [cm-1] to [rad/s]
nyu_p[4] = 233*factor # [cm-1] to [rad/s]

nyu_tau[0] = 15.53*factor  # [cm-1] to [rad/s]
nyu_tau[1] = 4.43*factor  # [cm-1] to [rad/s]
nyu_tau[2] = 0.42*factor # [cm-1] to [rad/s]
nyu_tau[3] = 54.14*factor # [cm-1] to [rad/s]
nyu_tau[4] = 12.94*factor # [cm-1] to [rad/s]

# function: begin
def begin():
    print ("begin")

# function: end
def end():
    print ("end")

# integrand of equation (3)
def eq1_sum(arg_i, nyu): # [rad/s] (parameters are converted to rad/s)
    numerator = np.power(nyu_p[arg_i],2.0)
    denominator = np.power(nyu_zero[arg_i],2.0)-np.power(nyu,2.0)-1j*(nyu_tau[arg_i])*nyu

    xk = numerator/denominator

    return xk

# equation (1) & equation (3)
def sio2_drude_lorentz(omega): # [rad/s]

    # 1st term
    function = ipu_inf

    # loop for dielectric funciton
    for i in range(0,number):
        function = function + eq1_sum(i, omega)

    return function

# time display
elapsed_time = time.time()-start
#print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
