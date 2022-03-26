# Dielectric function of Silicon
# used as external function
# Coded by Takuro TOKUNAGA
# Last modified: August 30 2017

import math
import numpy as np
import cmath
import time
from scipy.integrate import trapz, simps, quad, quadrature

start = time.time()

# parameters
counter = 0

# variables
listnumber1 = 615
listnumber2 = 2522
realfreq = [0]*listnumber1 # 0 to 614
realref = [0]*listnumber1 # 0 to 614
imagfreq = [0]*listnumber2 # 0 to 2521
imagref = [0]*listnumber2 # 0 to 2521

# open files
f1 = open('../dc/Si real frequency.txt','r') # do not move the txt file to other folder
f2 = open('../dc/Si real ref index.txt','r') # do not move the txt file to other folder
f3 = open('../dc/Si imag frequency.txt','r') # do not move the txt file to other folder
f4 = open('../dc/Si imag ref index.txt','r') # do not move the txt file to other folder

# read file
# file open read mode, real frequency
for line in f1:
    realfreq[counter] = float(line)
    #print(str(realfreq[counter]))
    counter = counter+1

# counter reset
counter = 0

# file open read mode, real refractive index
for line in f2:
    realref[counter] = float(line)
    #print(str(realref[counter]))
    counter = counter+1

# counter reset
counter = 0

# file open read mode, imag frequency
for line in f3:
    imagfreq[counter] = float(line)
    counter = counter+1

# counter reset
counter = 0

# file open read mode, imag refractive index
for line in f4:
    imagref[counter] = float(line)
    counter = counter+1

# end read file

# Dielectric function
def si(omega):
    # Refractive index of Si
    Si_n=3.4155
    Si_k=0.0014

    for j in range(listnumber1-1):
        if (omega>=realfreq[j] and omega<=realfreq[j+1]):
            Si_n=(omega-realfreq[j])/(realfreq[j+1]-realfreq[j])*(realref[j+1]-realref[j])+realref[j]
            #print(str(Si_n))

    for j in range(listnumber2-1):
        if (omega>=imagfreq[j] and omega<=imagfreq[j+1]):
            Si_k=(omega-imagfreq[j])/(imagfreq[j+1]-imagfreq[j])*(imagref[j+1]-imagref[j])+imagref[j]
            #print(str(Si_k))

    # Refractive index of Si
    ref_Si=complex(Si_n,Si_k)

    # Dielectric function of silicon
    Diel_Si=np.power(ref_Si,2)

    return Diel_Si

# file close
f1.close()
f2.close()
f3.close()
f4.close()

# display time
#elapsed_time = time.time()-start
#print("dc_si elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
