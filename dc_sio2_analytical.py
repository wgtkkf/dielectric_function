# Coded by Takuro TOKUNAGA
# Last modified: March 18 2019

import math
import numpy as np
import cmath
import time
start = time.time()

# define parameters
number = 4 # number of j change here
sigma=np.zeros((number), dtype='float64')
nyu_zero=np.zeros((number), dtype='float64')
nyu_p=np.zeros((number), dtype='float64')
nyu_tau=np.zeros((number), dtype='float64')


# Initialization of parameters
# SO14-200, Table 2
ipu_inf = 2.09
lx = ipu_inf # large x
factor = 3.33564*np.power(10., -11)/(2*np.pi)

sigma[0] = 76
sigma[1] = 125
sigma[2] = 63
sigma[3] = 170
#sigma[4] = 0

nyu_zero[0] = 1046
nyu_zero[1] = 1167
nyu_zero[2] = 1058
nyu_zero[3] = 798
#nyu_zero[4] = 0

nyu_p[0] = 575
nyu_p[1] = 288
nyu_p[2] = 459
nyu_p[3] = 415
#nyu_p[4] = 0

nyu_tau[0] = 1.55
nyu_tau[1] = 0.51
nyu_tau[2] = 10.57
nyu_tau[3] = 55.30
#nyu_tau[4] = 0

# function: begin
def begin():
    print ("begin")

# function: end
def end():
    print ("end")

# integrand of equation (3)
def eq1_sum(num_i, nyu):
    sa = np.sqrt(np.power(nyu,2)-1j*nyu_tau[num_i]*nyu)

    term1 = 1j*np.sqrt(np.pi)*np.power(nyu_p[num_i],2)/(np.sqrt(2)*2*sa*sigma[num_i])
    term2 = np.exp(-np.power(sa-nyu_zero[num_i],2)/(2*np.power(sigma[num_i],2)))
    term3 = special.erfc((sa-nyu_zero[num_i])/(1j*sigma[num_i]))
    term4 = np.exp(-np.power(sa+nyu_zero[num_i],2)/(2*np.power(sigma[num_i],2)))
    term5 = special.erfc((sa+nyu_zero[num_i])/(1j*sigma[num_i]))

    xk = term1*(term2*term3+term4*term5)

    return xk

# equation (1) & equation (3)
def sio2_analytical(omega):
    # rad/s into cm-1
    omega = factor*omega

    # 1st term
    function = ipu_inf

    # loop for dielectric funciton
    for i in range(0,number):
        function = function + eq1_sum(i, omega)

    return function

# time display
elapsed_time = time.time()-start
#print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
