# body code for activation of dielectric constant code
# Last modified: August 30 2017
# Coded by Takuro TOKUNAGA

import math
import numpy as np
import cmath
import time
import sys
from scipy.integrate import trapz, simps, quad, quadrature, romberg

#sys.path.append('../dc/')
from dc_dsi import drude # doped silicon
from dc_ndsi import ndsi # n-doped silicon
from dc_pdsi import pdsi # p-doped silicon
from dc_si import si # intrinsic silicon
from dc_si_adachi import si_adachi # intrinsic silicon by adachi model
from dc_sic import ipu # sic

start = time.time()

# loop number, adjust those values
nmax = np.power(10,4)

# Frequency
omegamin = 1.0*math.pow(10,12)#1.5*math.pow(10,14)
omegamax = 1.0*math.pow(10,16)#1.9*math.pow(10,14)
domega = (omegamax-omegamin)/nmax
omega = omegamin # Initialization of Frequency to omegamin

# functions for begin & finish
def begin():
    print ("begin")

def finish():
    print ("finish")

# main
begin()

# file open
f1 = open('data/dsi.txt', 'w')
f2 = open('data/ndsi.txt', 'w')
f3 = open('data/pdsi.txt', 'w')
f4 = open('data/si.txt', 'w')
f5 = open('data/si_adachi.txt', 'w')
f6 = open('data/sic.txt', 'w')

# first line
f1.write('omega dsi\n')
f2.write('omega ndsi\n')
f3.write('omega pdsi\n')
f4.write('omega si\n')
f5.write('omega si_adachi\n')
f6.write('omega sic\n')

# Frequency loop
for i in range(0,int(nmax)+1):

    df_dsi = drude(omega) # doped silicon
    df_ndsi = ndsi(omega) # n-doped silicon
    df_pdsi = pdsi(omega) # p-doped silicon
    df_si = si(omega) # intrinsic silicon
    df_si_adachi = si_adachi(omega) # intrinsic silicon by adachi model
    df_sic = ipu(omega) # sic

    # doped silicon
    f1.write(str(omega))
    f1.write(str(' '))
    f1.write(str(df_dsi.real))
    f1.write(str(' '))
    f1.write(str(df_dsi.imag))
    f1.write('\n')

    # n-doped silicon
    f2.write(str(omega))
    f2.write(str(' '))
    f2.write(str(df_ndsi.real))
    f2.write(str(' '))
    f2.write(str(df_ndsi.imag))
    f2.write('\n')

    # p-doped silicon
    f3.write(str(omega))
    f3.write(str(' '))
    f3.write(str(df_pdsi.real))
    f3.write(str(' '))
    f3.write(str(df_pdsi.imag))
    f3.write('\n')

    # silicon
    f4.write(str(omega))
    f4.write(str(' '))
    f4.write(str(df_si.real))
    f4.write(str(' '))
    f4.write(str(df_si.imag))
    f4.write('\n')

    # silicon by adachi model
    f5.write(str(omega))
    f5.write(str(' '))
    f5.write(str(df_si_adachi.real))
    f5.write(str(' '))
    f5.write(str(df_si_adachi.imag))
    f5.write('\n')

    # sic
    f6.write(str(omega))
    f6.write(str(' '))
    f6.write(str(df_sic.real))
    f6.write(str(' '))
    f6.write(str(df_sic.imag))
    f6.write('\n')

    # Frequency update
    omega = omega+domega

    # current step number
    #print("step:{:.0f}".format(i))

# file close
f1.close()
f2.close()
f3.close()
f4.close()
f5.close()
f6.close()

finish()

# time display
elapsed_time = time.time()-start
print("elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
