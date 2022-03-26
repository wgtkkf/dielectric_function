# Dielectric function of n-doped silicon
# S. Basu et al, J. Heat Transfer, 132 023301 (2010)
# Coded by Takuro TOKUNAGA
# Last modified: September 24 2019

import numpy as np
import time
start = time.time()

# Doping concentration
dpc = 1.0*np.power(10.,18) # cm^-3

# parmeters for drude model
mzero = 9.10938215*np.power(10.,-31)           # Electron Rest mass, [kg]
mstar = 0.27*mzero                  # Effective mass of an electron, [kg]
ec = 1.60217657*np.power(10.,-19)                 # electron charge, [C]
epuzero = 8.85418782*np.power(10.,-12) # Permittivity of free space, [F/m]
epuinf = 11.7               # Limiting value of Dielectric Function, [-]

# parameters for mobility, equation (2)
mumax=1414 # [cm^2/Vs]
mu1=68.5   # [cm^2/Vs]
mu2=56.1   # [cm^2/Vs]
cr=9.2*np.power(10.,16)  # [cm^-3]
cs=3.41*np.power(10.,20) # [cm^-3]
alpha=0.711              # [-]
beta=1.98                # [cm^-3]

# N-type ionization parameters & equation (4)
temperature = 300                      # temperature, [K]
theta = temperature/300 #        reduced temperature, [-]
la = 0.0824*np.power(theta,-1.622)      # large a =A, [-]
Nzero = 1.6*np.power(10.,18)*np.power(theta,0.7267) # [-]

if dpc<Nzero:
    lb = 0.4722*np.power(theta,0.0652) # large b = B, [-]
elif dpc>Nzero:
    lb = 1.23-0.3162*theta             # large b = B, [-]

xi = 1-la*np.exp(-np.power(lb*np.log(dpc/Nzero),2)) # [-], Eq. 4

# electron carrier concentration
ne = xi*dpc # cm^-3: [-]*[cm^-3]

# N-type Mobility
ntm =  mu1 + (mumax-mu1)/(1+np.power((ne/cr),alpha)) - mu2/(1+np.power((cs/ne),beta)) # Eq. 2, [cm^2/Vs]

# N-type Mobility, temperature consideration
ntm_tc = ntm*np.power(theta,1.5) # [cm^2/Vs]

# scattering Rate
gamma = ec/(mstar*(ntm_tc*np.power(10.,-4))) # R1, below Eq.1, V: [eV/(kg*(m^2/Vs))]

# plasma frequency
omegap2 = (ne*np.power(10.,6)*np.power(ec,2))/(mstar*epuzero) # below Eq.1, [m] , [m^-3*eV^2/(kg*(F*m^-1))]

# dielectric function
def ndsi(omega):
    y = epuinf-omegap2/(omega*(omega + 1j*gamma)) # Eq.1
    return y

# display time
#elapsed_time = time.time()-start
#print("dc_ndsi elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
