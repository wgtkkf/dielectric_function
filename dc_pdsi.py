# Dielectric function of P-doped silicon
# used as external function
# Coded by Takuro TOKUNAGA
# Last modified: August 30 2017

import math
import numpy as np
import cmath
import time
from scipy.integrate import trapz, simps, quad, quadrature
start = time.time()

# Doping concentration
dpconcentration = 9.0*math.pow(10,9) # N1

# parmeters for drude model
mzero = 9.10938215*math.pow(10,-31) # Electron Rest mass
mstar = 0.37*mzero # Effective mass of an electron
se = 1.60217657*math.pow(10,-19) # electron charge
epuzero = 8.85418782*math.pow(10,-12) # Permittivity of free space
epuinf = 11.7 # Limiting value of Dielectric Function

# N-type constants for Mobility
mue1=44.9
mue2=29.0
muemax=470.5
cr=2.23*math.pow(10,17)
cs=6.10*math.pow(10,20)
alpha=0.719
beta=2.0
pc = 9.23*math.pow(10,16)

# N-type ionization constants
temperature = 300 # temperature
redtemp = temperature/300 # reduced temperature
A1 = 0.2364*np.power(redtemp,-1.474)
N1 = 1.577*math.pow(10,18)*np.power(redtemp,0.46) # N01

if dpconcentration<N1:
    B1 = 0.433*np.power(redtemp,0.2213)
elif dpconcentration>N1:
    B1 = 1.268-0.338*redtemp

# ion concentration ratio
icr = 1-(A1*np.exp(-np.power(B1*np.log(dpconcentration/N1),2))) # I1

# hole carrier concentration
hcc = icr*dpconcentration # nh1

# P-type Mobility
ptm =  (mue1*np.exp(-pc/hcc))+(muemax/(1+np.power((hcc/cr),alpha)))-(mue2/(1+np.power((cs/hcc),beta))) # mue3001

# temperature consideration
tc = ptm*np.power((redtemp),1.5) # mue11

# Scattering Rate or damping
srod = se/(mstar*(tc/np.power(10,4))) # R1

# plasma frequency
pf = (hcc*np.power(10,6)*np.power(se,2))/np.sqrt((mstar*epuzero))

# dielectric function
def pdsi(omega):
    y = epuinf-((np.power(pf,2))/(omega*(omega + 1j*srod)))
    return y

# display time
#elapsed_time = time.time()-start
#print("dc_pdsi elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
