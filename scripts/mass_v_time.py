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


from __future__ import print_function
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



directory = sys.argv[1]
output = sys.argv[2]

ss_files = sorted(glob.glob(directory + "ss.[0-9]*[0-9]")) #get list of strings with full paths to each ss file

target, impactor, t = util.get_sorted_data(ss_files[0], units='cgs')


M_T = np.sum(target['mass'])
R_T = np.max(np.sqrt(target['x']**2 + target['y']**2 + target['z']**2))
r_T = np.mean(target['radius'])



#calculate approx kinetic and potential energy of each particle in each frame
f = open(output, 'w')

for frame in ss_files:
	sys.stdout.write('Current Frame: {0} of {1}\r'.format(frame, ss_files[-1]))
	sys.stdout.flush()
    #print("Current Frame: {0}".format(frame))
	data, t = util.get_data(frame, units='cgs')
	v2 = data['xdot']**2 + data['ydot']**2 + data['zdot']**2
	r = np.sqrt(data['x']**2 + data['y']**2 + data['z']**2)
	E_k = 0.5*data['mass']*v2
	E_pot = -cgs.G * data['mass']*M_T/r
	E = E_k + E_pot
	bound_ind = np.where((E <= 0.0) & (r > R_T + r_T)) #indices where particles is on bound orbit, but not touching surface (approx)
	unbound_ind = np.where(E > 0.0)
	N_esc = len(unbound_ind)
	M_esc = np.sum(data['mass'][unbound_ind])
	N_disk = len(bound_ind)
	M_disk = np.sum(data['mass'][bound_ind])

	in_disk = []

	if len(bound_ind) > 1:
		for particle in bound_ind:
			m = data['mass'][particle]
			r2 = np.array([data['x'][particle], data['y'][particle], data['z'][particle]])
			r1 = np.array([0,0,0]) #assumes target still at origin
			v1 = np.array([0,0,0])
			v2 = np.array([data['xdot'][particle], data['ydot'][particle], data['zdot'][particle]])
			OrbElem = util.orbital_elements(M_T, m, r1, r2, v1, v2)
			a = OrbElem[0]
			e = OrbElem[1]	

			pericenter = a*(1-e)
			if pericenter > R_T + r_T:
				in_disk += [particle]	
	


	N_disk_approx = len(in_disk)
	M_disk_approx = np.sum(data['mass'][in_disk])

	f.write("{0} {1} {2} {3} {4} {5} {6}\n".format(t, M_disk, M_esc, M_disk_approx, N_disk, N_esc, N_disk_approx))
f.close()














