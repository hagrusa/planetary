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
import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']#import os
import sys
import cgs as cgs #this is my own script that contains cgs constants
import utilities as util
#import re
import glob
#import glob

small = 10
medium = 14
large = 18


#numbers from paper
thetas = np.array([45, 45, 45, 45, 30, 60])
Mdisk_ratio = np.array([0.026, 0.029, 0.033, 0.033, 0.024, 0.017])




##############
#Plot one curve


#init_ss = sys.argv[1]#

#data_file = sys.argv[2]

#ss_files = sorted(glob.glob(directory + "ss.[0-9]*[0-9]")) #get list of strings with full paths to each ss file

#target, impactor, t = util.get_sorted_data(init_ss, units='cgs')#

##get basic data
#M_T = np.sum(target['mass'])
#R_T = np.max(np.sqrt(target['x']**2 + target['y']**2 + target['z']**2))
#r_T = np.mean(target['radius'])#

#M_imp = np.sum(impactor['mass'])
#R_imp = np.max(np.sqrt(impactor['x']**2 + impactor['y']**2 + impactor['z']**2))
#r_T = np.mean(impactor['radius'])#

#data = np.genfromtxt(data_file)
#t = data[:,0]
#M_disk = data[:,1]
#M_esc = data[:,2]#
#

#M_tot = M_imp + M_T
#plt.plot(t, M_disk/M_tot, label = 'M_disk/M_imp')
#plt.plot(t, M_esc/M_tot, label = 'M_esc/M_imp')
#plt.legend(loc=0)
##plt.show()




#######################
#Plot multiple curves
init_ss = sys.argv[1]
directory = sys.argv[2]


target, impactor, t = util.get_sorted_data(init_ss, units='cgs')

#get basic data
M_T = np.sum(target['mass'])
R_T = np.max(np.sqrt(target['x']**2 + target['y']**2 + target['z']**2))
r_T = np.mean(target['radius'])

M_imp = np.sum(impactor['mass'])
R_imp = np.max(np.sqrt(impactor['x']**2 + impactor['y']**2 + impactor['z']**2))
r_T = np.mean(impactor['radius'])

M_tot = M_imp + M_T
data_files = sorted(glob.glob(directory + "mass_time_b[0-9][0-9].txt")) #get list of strings with full paths to each ss file
theta = np.arange(0,86, 5)
M_esc_final = np.empty(len(theta))
M_disk_final = np.empty(len(theta))

color_idx = np.linspace(0, 1, len(theta))
plt.style.use('dark_background')



###########################################
########## Disk Mass vs Time ##############
###########################################
for i in range(0, len(data_files)):
	data = np.genfromtxt(str(data_files[i]))
	t = data[:,0]
	M_disk = data[:,1]
	M_esc = data[:,2]
	M_esc_final[i] = M_esc[-1]
	M_disk_final[i] = M_disk[-1]
	plt.plot(t/60, M_disk/M_imp, color = plt.cm.hsv(color_idx[i]), label = r'{0}'.format(theta[i]))
plt.legend(loc=0, fontsize = 'x-small', title = r'$\text{\bf Impact Angle [Degrees]}$', ncol=3)
plt.title(r'$\text{\bf Disk Mass vs. Time}$', fontsize=large)
plt.ylabel(r'$\mathbf{\frac{M_{disk}}{M_{imp}}}$', rotation=0, fontsize = large, labelpad=20)
plt.xlabel(r"$\text{\bf Time [minutes]}$", fontsize=medium)
#plt.show()
plt.tight_layout()

plt.savefig('disk_mass_vs_time.pdf')
plt.close()

###########################################
##########  Esc Mass vs Time ##############
###########################################
for i in range(0, len(data_files)):
	data = np.genfromtxt(data_files[i])
	t = data[:,0]
	M_disk = data[:,1]
	M_esc = data[:,2]

	plt.plot(t/60, M_esc/M_imp, color = plt.cm.hsv(color_idx[i]), label = r'{0}'.format(theta[i]))
plt.legend(loc=0, fontsize = 'x-small', title = r'$\text{\bf Impact Angle [Degrees]}$', ncol=3)
plt.title(r'$\text{\bf Escaping Mass vs. Time}$', fontsize=large)
plt.ylabel(r'$\mathbf{\frac{M_{disk}}{M_{imp}}}$', rotation=0, fontsize = large, labelpad=20)
plt.xlabel(r"$\text{\bf Time [minutes]}$", fontsize=medium)
#plt.show()
plt.tight_layout()

plt.savefig('esc_mass_vs_time.pdf')
plt.close()


###########################################
########## Final Mass vs Angle ##############
###########################################
fig, ax1 = plt.subplots()

#Simulation Results
ax1.plot(theta, M_disk_final/M_imp, 'o', color = 'b',  label = 'Final Disk Mass')
#ax1.set_xlabel('time (s)')
# Make the y-axis label, ticks and tick labels match the line color.
#ax1.set_ylabel('exp', color='b')
ax1.set_ylabel(r'$\mathbf{\frac{M_{disk}}{M_{imp}}}$', rotation=0, fontsize = large, labelpad=20)
#ax1.tick_params('y', colors='b')

#Paper Restuls
ax2 = ax1.twinx()
ax2.plot(thetas, Mdisk_ratio, 'o', color='g', label = 'Final Disk Mass [Citron 2015]')
ax2.set_ylabel(r'$\mathbf{\frac{M_{disk}}{M_{imp}}}$', rotation=0, fontsize = large, labelpad=20)

#ax2.tick_params('y', colors='r')

#fig.tight_layout()

#plt.plot(theta, M_esc_final/M_imp, 'o', label = 'Final Escaped Mass')
#plt.plot(theta, M_disk_final/M_imp, 'o', label = 'Final Disk Mass')
plt.title(r'$\text{\bf Final Disk Mass vs. Impact Angle}$', fontsize=large)
#plt.ylabel(r'$\mathbf{\frac{M}{M_{imp}}}$', rotation=0, fontsize = large, labelpad=20)
plt.xlabel(r"$\text{\bf Impact Angle [degrees]}$", fontsize=medium)
plt.legend(loc=0, fontsize = 'small')
plt.grid(linestyle='-', linewidth='0.5')
#plt.show()
plt.tight_layout() 
plt.savefig('mass_vs_theta.pdf')
plt.close()



