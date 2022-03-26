# 3rd order spline interporation
# Coded by Takuro TOKUNAGA
# Last modified: June 11 2018

import math
import numpy as np
import cmath
import time
import pandas as pd
import scipy.linalg as linalg
import matplotlib.pyplot as plt
from scipy.integrate import trapz, simps, quad, quadrature, romberg
start = time.time()

# define parameters
number = 4 # number of j change here
sigma=np.zeros((number), dtype='float64')
nyu_zero=np.zeros((number), dtype='float64')
nyu_p=np.zeros((number), dtype='float64')
nyu_tau=np.zeros((number), dtype='float64')

# DE method parameters
nmax = np.power(10.,3)
sn = 1 # n
sh = np.power(10.,-3) # h
conv = np.power(10.,-4) # convergence criteria
dif = 1

# Initialization of parameters
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
nyu_tau[3] = 55.3
#nyu_tau[4] = 0

# function: begin
def begin():
    print ("begin")

# function: end
def end():
    print ("end")

# integrand of equation (3)
def integrand(u, num_i, nyu): # x to u
    term1 = np.power(u-(nyu_zero[num_i]),2)/(2*np.power(sigma[num_i],2))
    denominator = np.power(u,2)-np.power(nyu,2)+(1j*nyu_tau[num_i]*nyu)
    term2 = np.power(nyu_p[num_i],2)/denominator

    integrand = np.exp(-term1)*term2

    return integrand

# DE method, -inf to +inf type
def phi(sn):
    y = np.sinh((0.5*np.pi)*np.sinh(sn*sh))
    return y

def dphi(sn):
    y = (0.5*np.pi)*np.cosh(sn*sh)*np.cosh((0.5*np.pi)*np.sinh(sn*sh))
    return y

# equation (1) & equation (3)
def sio2(omega):
    # rad/s into cm-1
    omega = factor*omega
    #print(str(omega))

    # Initialization of old, sn & dif
    old = 0
    sn = 1
    dif = 1

    # 1st term
    function = ipu_inf

    # loop for dielectric funciton
    for i in range(0,number):
        # Initialization of new
        new = integrand(phi(0), i, omega)*dphi(0)

        # integral calculation
        while dif>conv:
            new = new + integrand(phi(-sn), i, omega)*dphi(-sn)\
            + integrand(phi(sn), i, omega)*dphi(sn)

            # conv check
            dif = abs(new-old)
            #print(str(dif))

            sn = sn+1
            old = new

            if dif < conv:
                break
                
        # integral result
        integral = (sh*new)

        # update dielectric function
        function = function + (1/(np.sqrt(2*np.pi)*sigma[i]))*integral

        # Reset parameters
        old = 0
        sn = 1
        dif = 1

    return function

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
