##########################################
# File name: 
# Author: Harrison Agrusa, hagrusa@astro.umd.edu
# Date created: 
# Date last modified:
# Description: 
# INPUTS:
#
#
#
# OUTPUTS:
#
#
#
##########################################


#from __future__ import print_function
#import psutil
#from joblib import Parallel, delayed
#import multiprocessing
import numpy as np
import matplotlib.pyplot as plt
#import os
#import sys
import cgs as cgs #this is my own script that contains cgs constants
#import dart_utilities as util
#import re
#import glob

N = np.array([5e4, 1e5, 3e5, 1e6, 3e5, 3e5, 3e5, 3e5, 3e5, 3e5])
Mimp = 1000*1e22 * np.array([1.68, 1.68, 1.68, 1.68, 1.68, 1.68, 1.68, 3.5, 6, 10]) #grams
alpha = np.array([1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 2, 1.4, 1.4, 1.9])
theta = np.array([45, 45, 45, 45, 30, 60, 45, 45, 45, 45])
Eimp = 1e7 * 1e29 * np.array([3, 3, 3, 3, 3, 3, 6.3, 6, 10, 30]) #ergs
gamma  = np.array([0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.058, 0.1, 0.167])
Mdisk_ratio = np.array([0.026, 0.029, 0.033, 0.033, 0.024, 0.017, 0.024, 0.033, 0.023, 0.013])
Mdisk = Mdisk_ratio*Mimp

rho_imp = 3 #basalt, g/cc
rho_core = 7.874# iron, g/cc
rhoe_mant = 3.5 #olivine, g/cc

Mt = 6e23 * 1000 #mass of target
Mcore = 0.3*Mt
Mmant = 0.8*Mt

Rcore = (Mcore/(4*np.pi*rho_core/3))**(1./3)
Rmant = (Mmant/(4*np.pi*rhoe_mant/3) + Rcore**3)**(1./3)
Rt = Rmant

Mimp = 1000*1e22 * 1.68
Rimp = (Mimp/(4*np.pi*rho_imp/3))**(1./3)


#theory
angles = np.linspace(30, 60, 50)
b = np.sin(np.radians(angles))
Y = (Rimp + Rt)*b

d = Y
R = Rt
r = Rimp

#V_interact = np.pi*(R + r - d)**2 *(d**2 + 2*d*r - 3*r**2 + 2*d*R + 6*r*R - 3*R**2) /(12*d)
B = (R+r)*np.sin(np.radians(angles))
l = R+r-B
M_interact = Mimp*(3*r*l**2 - (np.pi/3)*l**3)/(4*r**3)

A = 0.28
a = -2.7
b = 3.7
c = -2.7

def f(y):
	value = y*(1-y)*(y**(1./3) + (1-y)**(1./3))**0.5
	return value

M_pass = Mimp - M_interact
#print(M_pass/Mimp)

gamma = 0.028
C_y = 2.8*(0.1/gamma)**1.25
Jcol = np.sqrt(2)*f(gamma)*np.sin(np.radians(angles))*1.4


Morb1 = Mt*C_y*(M_pass/Mt)**2
Morb2 = Mt*A*gamma**a * (M_pass/Mt)**b * 1.4**c
Morb3 = Mt * 2.15*Jcol**2.79
plt.plot(angles, Morb1/Mt, label = 'Theory')
plt.plot(angles, Morb2/Mt)
plt.plot(angles, Morb3/Mt)
plt.plot(theta[:7], Mdisk[:7]/Mt, 'o', label = 'Simulations')
plt.legend(loc=0)
plt.show()































