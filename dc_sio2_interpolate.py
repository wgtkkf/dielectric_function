# Dielectric function of Silica
# used as external function
# Coded by Takuro TOKUNAGA
# Last modified: iune 12 2018

import math
import numpy as np
import cmath
import time
from scipy.integrate import trapz, simps, quad, quadrature

start = time.time()

# parameters
counter = 0

# variables
listnumber1 = 106 # real part
listnumber2 = 105 # imaginary part
realfreq = [0]*listnumber1
dc_real = [0]*listnumber1
imagfreq = [0]*listnumber2
dc_img = [0]*listnumber2

# open files
f1 = open('../dc/Sio2 real frequency.txt','r') # do not move the txt file to other folder
f2 = open('../dc/Sio2 real dc.txt','r') # do not move the txt file to other folder
f3 = open('../dc/Sio2 imag frequency.txt','r') # do not move the txt file to other folder
f4 = open('../dc/Sio2 imag dc.txt','r') # do not move the txt file to other folder

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
    dc_real[counter] = float(line)
    #print(str(dc_real[counter]))
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
    dc_img[counter] = float(line)
    counter = counter+1

# end read file

# Dielectric function
def sio2_interpolate(omega):

    for i in range(listnumber1-1):
        if (omega>=realfreq[i] and omega<=realfreq[i+1]):
            Sio2_real=(omega-realfreq[i])/(realfreq[i+1]-realfreq[i])*(dc_real[i+1]-dc_real[i])+dc_real[i]

    for i in range(listnumber2-1):
        if (omega>=imagfreq[i] and omega<=imagfreq[i+1]):
            Sio2_img=(omega-imagfreq[i])/(imagfreq[i+1]-imagfreq[i])*(dc_img[i+1]-dc_img[i])+dc_img[i]

    Diel_Sio2=Sio2_real+1j*Sio2_img

    return Diel_Sio2

# file close
f1.close()
f2.close()
f3.close()
f4.close()

# display time
#elapsed_time = time.time()-start
#print("dc_si elapsed_time:{:.2f}".format(elapsed_time) + "[sec]")
