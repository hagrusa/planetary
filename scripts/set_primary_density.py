##########################################
# File name: carve.py
# Author: Harrison Agrusa, hagrusa@astro.umd.edu
# Date created: 
# Date last modified:
# Description: carve rubble pile to sphere of desired radius
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
#import matplotlib.pyplot as plt
#import os
import sys
#import cgs as cgs #this is my own script that contains cgs constants
import utilities as util
#import re
#import glob


rp = sys.argv[1]

N = int(sys.argv[2])#desired number of particles
rp, t = util.get_data(rp, units='raw')

r = np.sqrt(rp["x"]**2 + rp["y"]**2 + rp["z"]**2)
R = np.max(r)

for i in range(0, len(r)):
	if r <= 0.3*R: #core
		

def calc_N(rp, R):
	r = np.sqrt(rp["x"]**2 + rp["y"]**2 + rp["z"]**2)
	indices = np.where(r<R)
	N = len(rp['index1'][indices])
	return N

def carve(rp, R):
	r = np.sqrt(rp["x"]**2 + rp["y"]**2 + rp["z"]**2)
	indices = np.where(r<R)
	util.write_data(rp, 'carved', indices)
	return 1

print(rp['x'])


r = np.sqrt(rp["x"]**2 + rp["y"]**2 + rp["z"]**2)

R = np.max(r)
while calc_N(rp, R) > N:
	R *= 0.99 #decrease by 1%

carve(rp, R)

























