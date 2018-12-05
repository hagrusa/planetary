##########################################
# File name: restartDEM.py
# Author: Harrison Agrusa
# Date created: 8/23/2018
# Date last modified:
# Description: rewrite ss.par file for DEM restart
# How to Run: python restartDEM.py 
##########################################

# THE LATEST OUTPUT FILE MUST ALSO MATCH THE LATEST DEM FILE

from __future__ import print_function
import numpy as np
import os
import sys
import glob
import subprocess
import argparse
import utilities as util


parser = argparse.ArgumentParser(description='Rewrite ss.par for DEM restart. (current ss.par file will be renamed to old_ss.par)')

parser.add_argument('--sspar', help = 'path to current ss.par file ')
parser.add_argument('--nSteps',  help = 'Number of additional steps for run, default is to just double original nSteps')


args=parser.parse_args()


if args.sspar:
	sspar = args.sspar
else:
	sspar = 'ss.par' 

dem_files = sorted(glob.glob("ss.[0-9]*[0-9].dem")) #get list of strings with full paths to each dem file
ss_files = sorted(glob.glob("ss.[0-9]*[0-9]"))

#find latest dem file
last_dem = dem_files[-1]
#make copy of last dem file
#subprocess.call("cp {0} {1}".format(last_dem, "backup_" + last_dem), shell=True) #make backup copy

last_dem_id = last_dem[:-4] #cut off .dem
last_ss = last_dem_id #has form ss.00001234 this is the file to restart from
last_dem_id = last_dem_id[3:] #cut off ss.
last_dem_id = int(last_dem_id) #integer value

last_ss = ss_files[-1]

if args.nSteps:
	nSteps = args.nSteps
else:
	nSteps = 2*last_dem_id #double number of steps


subprocess.call("cp {0} {1}".format(sspar, "old_ss.par"), shell=True) #make backup copy

input_file = last_ss #restart file
read_DEM = 1 #variable to set for ss.par to read DEM file
start_step = last_dem_id

#rewrite ss.par to be set up for DEM restart
#dummy = util.edit_sspar(sspar, iStartStep = start_step, nSteps=nSteps, bReadDEMData=read_DEM, achInFile=input_file)

dummy = util.edit_sspar(sspar, iStartStep = start_step, bReadDEMData=read_DEM, achInFile=input_file)












