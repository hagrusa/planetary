##########################################
# File name: 
# Author: Harrison Agrusa, hagrusa@astro.umd.edu
# Date created: 
# Date last modified:
# Description: make text file of disk mass & escaped mass vs time
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
import sys
import cgs as cgs #this is my own script that contains cgs constants
import utilities as util
#import re
import glob
#import glob



init_ss = sys.argv[1]
data_file = sys.argv[2]

#ss_files = sorted(glob.glob(directory + "ss.[0-9]*[0-9]")) #get list of strings with full paths to each ss file

target, impactor, t = util.get_sorted_data(init_ss, units='cgs')

#get basic data
M_T = np.sum(target['mass'])
R_T = np.max(np.sqrt(target['x']**2 + target['y']**2 + target['z']**2))
r_T = np.mean(target['radius'])

M_imp = np.sum(impactor['mass'])
R_imp = np.max(np.sqrt(impactor['x']**2 + impactor['y']**2 + impactor['z']**2))
r_T = np.mean(impactor['radius'])

data = np.genfromtxt(data_file)
t = data[:,0]
M_disk = data[:,1]
M_esc = data[:,2]


M_tot = M_imp + M_T
plt.plot(t, M_disk/M_tot, label = 'M_disk/M_imp')
plt.plot(t, M_esc/M_tot, label = 'M_esc/M_imp')
plt.legend(loc=0)
plt.show()



















