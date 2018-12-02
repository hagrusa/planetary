##########################################
# File name: init_norm_primphase.py
# Author: Harrison Agrusa
# Date created: 4/10/2018
# Date last modified: 10/11/2018
# INPUTS:
#     initstate_file = initial conditions file in special format for benchmarking
#     systemdata_file = file containing onther important info - mass of each body, constants etc
#     main_ss = ss_file of primary (this should be already carved out with carvepoly)
#     moon_ss = ss_file of secondary, already carved
# OUTPUT:
#     final ss file to use as initial condition for pkdgrav run
# Description/procedure:
#		1. Read in initial conditions/ss files
#       2. Build rotation matrices from initial conditions
#       3. Set orientation of primary & secondary particle-by-particle
#       4. Save new ss files
#       5. write bash script that will call rpx
#       6. execute rpx bash script, this will set the other condtions 
#       7. rpx script saves final output as init.ss
# NOTES/OTHER:
#   There are 3 frames in this problem: 
#		-Inertial (V frame)
#		-primary body frame (V' or Vp frame)
#		-secondary body frame (V'' or Vs frame)
#   In general:
# 				variable_A refers to primary
# 				variable_B refers to secondary
#	units: All units get normalized - 
#			they are normalized based on the mass/time/length unit specified in systemdata file
##########################################

from __future__ import print_function
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import dart_utilities as util
#from mpl_toolkits.mplot3d import Axes3D
import sys
import subprocess as subp
import time
np.seterr(divide='ignore') #ignore divide by zero errors



target = sys.argv[1] 
impactor = sys.argv[2]
theta = np.deg2rad(float(sys.argv[3])) 

R_imp = 75 #km
R_t = 488 #km
x = -(R_imp + R_t + 5) #5km for buffer space
b = np.sin(theta)
y = b *(R_imp + R_t)
vx = 710 #m/s


##############################################
######## WRITE BASH SCRIPT TO RUN RPX ########
##############################################

bash_script = "temp_rpx.sh"
f = open(bash_script, 'w')
f.write('#!/bin/sh \n')
f.write('touch rpx.log init.ss \n')
#f.write('\n')
#call rpx
f.write('rpx << EOF\n')
#Overwrite log file [Y/n]? 
f.write('\n') #default yes
# Use simulation units (AU, M_sun, etc.) [y/N]? 
f.write('\n') #default no

#######################################
############### primary ###############
#######################################

#File 1 [or RETURN to quit]:
f.write(str(target) + '\n') #enter main.ss file (primary)
#Number of particles = 17000 (time 0.00222413)
#1. Total mass = 4.66444e-18 M_sun
#2. Bulk radius = 8.10743e-09 AU
#   [Bulk semi-axes: 7.73259e-09 7.56302e-09 7.55147e-09 AU]
#3. Bulk density = 2.5215e+06 M_sun/AU^3
#4. Centre-of-mass position = -7.28308e-11 9.23819e-12 2.21702e-11 AU
#5. Centre-of-mass velocity = 5.19502e-11 -9.89897e-11 3.74422e-11 x 30 km/s
#6. Orientation: a1 = 0.491202 -0.348738 0.798187
#                a2 = 0.505301 -0.632322 -0.587231
#                a3 = 0.709501 0.691774 -0.134380
#7. Spin = 0.000820355 -0.00610639 -0.00133202 x 2pi rad/yr
#   [Ang mom = 7.6975e-38 -6.01312e-37 -1.32847e-37 (sys units)]
#   [Effective spin = 0.00630345 x 2pi rad/yr]
#   [Rotation index = 0.19 (SAM)]
#8. Color = 3 (green)
#9. Aggregate ID = N/A
#10. Particle ID


#######################################
############## secondary ##############
#######################################

#Enter number to change (or 0 to continue): 
f.write('0\n')
#File 2 [or RETURN to quit]: 
f.write(str(impactor) + '\n')

### set COM position ###
f.write('4\n')
#Specify absolute position [Y/n]? 
f.write('\n') #default yes
#Enter new position [x y z in km]: 
f.write("{0} {1} {2}\n".format(x, y, 0)) 

### set COM velocity ###
f.write('5\n')
#Specify absolute velocity [Y/n]?
f.write('\n')
#Enter new velocity [vx vy vz in m/s]:
f.write("{0} {1} {2}\n".format(vx, 0, 0)) #velocity = momentum/mass

### done ###
#Enter number to change (or 0 to continue):
f.write('0\n')
#File 3 [or RETURN to quit]: 
f.write('\n')
#Recenter barycentric position and velocity [Y/n]? 
f.write('\n') #default yes
#Output file [default moon.ss]:
f.write('init.ss\n')
#Output file already exists...overwrite [y/N]?
f.write('y\n')
f.write('EOF')
f.close()

############################################
############# CALL BASH SCRIPT #############
############################################

subp.call('chmod +x ' + str(bash_script), shell=True)
subp.call('./'+str(bash_script), shell=True, stdout=subp.DEVNULL)

########################################################################################
#								END OF SCRIPT
########################################################################################
