##########################################
# File name: utilities.py
# Author: Harrison Agrusa, hagrusa@astro.umd.edu
# Date created: 2/15/2018
# Date last modified: 10/17/2018
# Description: Useful functions for doing analysis
##########################################


import numpy as np
import subprocess
import cgs as cgs
import ssio as ss
import scipy

##########################################
# To use thes functions, the 'data dictionaries' must have the following keys:
# ["index1", "index2", "mass", "radius", "x", "y", "z", \
# "xdot", "ydot", "zdot", "spin_x", "spin_y", "spin_z", "color"])
# dictionary["x"] should give you a 1D array of the x positions of particles

# WARNING:
# These functions assume that everything is in cgs units!!!!!

# WARNING #2:
# cgs.py is imported above, this is my own script that just defines
# various constants in cgs. Either make your own, or ask me for mine!
##########################################


length_to_cm = cgs.AU  # convert system units to cm
time_to_sec = (1. / (2 * np.pi)) * 365.25 * 24 * 3600  # yr/2pi to seconds
velocity_to_cms = length_to_cm / time_to_sec  # convert system units to cm/s
mass_to_g = cgs.Msun  # convert to grams

###############################################################################
#
# get_sorted_data
# author       : Harrison Agrusa
# date         : 2/15/2018
# Description  : function to create two dictionaries from a single ss frame
# Invocation   : get_sorted_data ( full_path_to_ssfile (string) )
# Input        : full_path_to_ssfile (required)
# Returns      : rp1 (dictionary), rp2 (dictionary), timestamp
# Notes        : use ONLY for 2 rubble piles that must be delineated by color
#              : rp1 will be the more massive rubble pile
###############################################################################

def get_sorted_data(ss_file, units='None'):
    # given path to ss_file,
    # construct 2 data dictionaries and return them along with the timestamp
    # (timestamp is scalar)
    ######### IMPORTANT: ###########
    # This function assumes that the two rubble piles are distinguished by their color
    keys = np.array(["index1", "index2", "mass", "radius", "x", "y", "z",
                     "xdot", "ydot", "zdot", "spin_x", "spin_y", "spin_z", "color"])
    #for flipping x and y:
    #keys = np.array(["index1", "index2", "mass", "radius", "y", "x", "z",
    #                 "ydot", "xdot", "zdot", "spin_y", "spin_x", "spin_z", "color"])
    
    num_keys = len(keys)

    rp1 = {}
    rp2 = {}

    header, data_array = ss.read_SS(ss_file, headout="yes, thank you very much")  # always be polite to your code
    timestamp = header[0]

    colors = np.unique(data_array[-1, :])  # color is the last element
    if len(colors) > 2:
        print("more than two colors/rubble-piles")
        quit()
    elif len(colors) == 1:
        print("only 1 rubble-pile, don't use this function")
        quit()
    elif len(colors) != 2:
        print("something wrong with rubble piles, read through the get_sorted_data function")
        quit()

    array1 = np.empty((0, num_keys))
    array2 = np.empty((0, num_keys))

    for i in range(len(data_array[0, :])):  # iterate through particles
        # sort based on color
        if data_array[-1, i] == colors[0]:
            array1 = np.append(array1, [data_array[:, i]], axis=0)
        elif data_array[-1, i] == colors[1]:
            array2 = np.append(array2, [data_array[:, i]], axis=0)
        else:
            print("something is wrong")
            quit()
    del data_array  # free memory
    for i in range(0, num_keys):
        # create data dictionaries
        rp1[keys[i]] = array1[:, i]
        rp2[keys[i]] = array2[:, i]
    if np.sum(rp1["mass"]) < np.sum(rp2["mass"]):
        # make sure that rp1 is the main body (most massive)
        rp1, rp2 = rp2, rp1
    del array1
    del array2
    # convert to the BEST (cgs) units
    if units == 'cgs':
        rp1, timestamp = convert_to_cgs(rp1, timestamp)
        rp2 = convert_to_cgs(rp2)
    elif units == 'SI' or units == 'mks':
        rp1, timestamp = convert_to_mks(rp1, timestamp)
        rp2 = convert_to_mks(rp2)
    elif units == 'raw':
        pass
    else:
        print("Please choose desired units: cgs, mks or raw")
        print("Usage: get__sorted_data(path/to/ss/file, units = 'cgs')")
        return None
    return rp1, rp2, timestamp


################################################################################################
##
## get_data
## author       : Harrison Agrusa
## date         : 2/15/2018
## Description  : get data dictionary of all particles
## Invocation   : get_data(path_to_ssfile)
## Input        : path to ss file (string)
## Returns      : data dictionary (see keys below), and timestamp. Qauntites returned in desired units
## Notes        : options for units are 'raw', 'mks', and 'cgs', must choose 1
################################################################################################
def get_data(ss_file, units=None):
    # given path to ss_file, construct data dictionary
    # return dictionary are timestamp (timestamp is scalar)
    # DOES NOT SEPARATE RUBBLE PILES
    keys = np.array(["index1", "index2", "mass", "radius", "x", "y", "z",
                     "xdot", "ydot", "zdot", "spin_x", "spin_y", "spin_z", "color"])
    num_keys = len(keys)
    dictionary = {}
    header, data_array = ss.read_SS(ss_file, headout="yes please")
    timestamp = header[0]
    for j in range(0, num_keys):
        dictionary[keys[j]] = data_array[j, :]
    if units == 'cgs':
        dictionary, timestamp = convert_to_cgs(dictionary, timestamp)
    elif units == 'SI' or units == 'mks':
        dictionary, timestamp = convert_to_mks(dictionary, timestamp)
    elif units == 'raw':
        pass
    else:
        print("Please choose desired units: cgs, mks or raw")
        print("Usage: get_data(path/to/ss/file, units = 'cgs')")
        return None
    return dictionary, timestamp

################################################################################################
## write_data
## 
## author       : Harrison Agrusa
## date         : 2/15/2018
## Description  : convert data dictionary back to 14xN numpy array and write as ss file
## Invocation   : write_data(data_dictionary, "my_file.ss", t = current_time)
## Input        : data dictionary [python dictionary], desired output file name [str], t (optional) [float]
## Returns      : just returns integer 1
## Notes        : 
################################################################################################

def write_data(data, outname, indices, t = 0.0):
    N = len(data["index1"][indices])
    data_array = np.empty((14, N))
    keys = np.array(["index1", "index2", "mass", "radius", "x", "y", "z",
                     "xdot", "ydot", "zdot", "spin_x", "spin_y", "spin_z", "color"])
    for i in range(0, 14):
        data_array[i,:] = data[keys[i]][indices]
    ss.write_SS(data_array, outname, t)
    return 1

################################################################################################
## normalize_units
## 
## author       : Harrison Agrusa
## date         : 2/15/2018
## Description  : take data dictionary (in mks units), and normalize to didymos system
## Invocation   : 
## Input        : 
## Returns      : 
## Notes        : Default values for DU MU and TU are the benchmarking activities
## WARNINGS     : This function assumes data dictionary in mks units, then converts to normalized units
################################################################################################


def normalize_units(data, timestamp=None , DU=1180.0, MU=5.27784750960427490234375e11, TU=6829.65691795942):
    VU = DU/TU
    data["x"] /= DU #convert AU to meters, then normalize
    data["y"] /= DU
    data["z"] /= DU
    data["xdot"] /= VU  
    data["ydot"] /= VU  
    data["zdot"] /= VU
    data["spin_x"] /=  (1./TU) 
    data["spin_y"] /=  (1./TU) 
    data["spin_z"] /=  (1./TU) 
    data["mass"] /= MU
    data["radius"] /= DU
    if timestamp != None:
        timestamp /= TU
        return data, timestamp
    else:
        return data


################################################################################################
## unnormalize_units
## 
## author       : Harrison Agrusa
## date         : 2/15/2018
## Description  : take normalized data dictionary (didymos units) and convert to mks
## Invocation   : 
## Input        : 
## Returns      : 
## Notes        : 
################################################################################################


def un_normalize_units(data, timestamp=None , DU=1180.0, MU=5.27784750960427490234375e11 , TU=6829.65691795942 ):
    VU = DU/TU
    data["x"] *= DU #convert AU to meters, then normalize
    data["y"] *= DU
    data["z"] *= DU
    data["xdot"] *= VU  #this gets changed later anyways, shouldnt matter
    data["ydot"] *= VU  
    data["zdot"] *= VU
    data["spin_x"] *= (1./TU) #this gets changed later anyways, shouldnt matter
    data["spin_y"] *= (1./TU) 
    data["spin_z"] *= (1./TU) 
    data["mass"] *= MU #this gets changed later anyways, shouldnt matter
    data["radius"] *= DU
    if timestamp != None:
        #timestamp *= time_to_sec / 3600
        timestamp *= TU
        return data, timestamp
    else:
        return data

################################################################################################
##
## convert_to_cgs
## author       : Harrison Agrusa
## date         : 2/15/2018
## Description  : takes data dictionary that has pkdgrav units and converts to cgs
## Invocation   : convert_to_cgs(dictionary, timestamp)
## Input        : data dictionary w/ pkdgrav units, and timestamp (optional)  units
## Returns      : returns same quantities in cgs units (the superior unit system)
## Notes        : 
################################################################################################

def convert_to_cgs(data, timestamp=None, length=cgs.AU, mass=cgs.Msun, time=365.25*24*3600*0.5/np.pi):
    velocity = length/time
    data["x"] *= length  # convert to cm
    data["y"] *= length  # convert to cm
    data["z"] *= length  # convert to cm
    data["xdot"] *= velocity  # convert to cm/s
    data["ydot"] *= velocity  # convert to cm/s
    data["zdot"] *= velocity  # convert to cm/s
    data["spin_x"] *= (1. / time)  # convert to 1/s
    data["spin_y"] *= (1. / time)  # convert to 1/s
    data["spin_z"] *= (1. / time)  # convert to 1/s
    data["mass"] *= mass
    data["radius"] *= length
    if timestamp != None:
        timestamp *= time
        return data, timestamp
    else:
        return data

################################################################################################
##
## convert_to_mks
## author       : Harrison Agrusa
## date         : 2/15/2018
## Description  : takes data dictionary that has pkdgrav units and converts to mks
## Invocation   : convert_to_cgs(dictionary, timestamp)
## Input        : data dictionary w/ system units, and timestamp (optional) in system units
## Returns      : returns same quantities in cgs units (the superior unit system)
## Notes        : 
################################################################################################

def convert_to_mks(dictionary, timestamp=None, length=cgs.AU/100, mass=cgs.Msun/1000, time=365.25*24*3600*0.5/np.pi):
    velocity = length/time
    dictionary["x"] *= length
    dictionary["y"] *= length
    dictionary["z"] *= length
    dictionary["xdot"] *= velocity 
    dictionary["ydot"] *= velocity 
    dictionary["zdot"] *= velocity
    dictionary["spin_x"] *= (1. / time)
    dictionary["spin_y"] *= (1. / time)  
    dictionary["spin_z"] *= (1. / time)  
    dictionary["mass"] *= mass
    dictionary["radius"] *= length
    if timestamp != None:
        timestamp *= time
        return dictionary, timestamp
    else:
        return dictionary


################################################################################################
##
## totE_sslog
## author       : Harrison Agrusa
## date         : 2/15/2018
## Description  : get 'total energy' from sslog file
## Invocation   : 
## Input        : path to sslog file (string)
## Returns      : 2 1-d arrays - the time and energy
## Notes        : This isn't a true total energy because sslog excludes rotational energy
##              : just good for a quick sanity check
################################################################################################

def totE_sslog(path_to_sslog):
    # given path to sslog file, return time and total energy array
    # sslog files do not contain rotational energy
    # so this plot is only an approximation, but useful because fast
    sslog_data = np.genfromtxt(path_to_sslog)
    times = sslog_data[:, 0]
    E_tot = sslog_data[:, 2]
    del sslog_data
    return times, E_tot

################################################################################################
##
## calc_I
## author       : Harrison Agrusa
## date         : 2/15/2018
## Description  : calculates inerta tensor of a rubble pile
## Invocation   : 
## Input        : data dictionary (dict), center of mass (3 element array)
## Returns      : inertia tensor
## Notes        : 
################################################################################################

def calc_I(data, com = np.array([0, 0, 0])):
    x = data["x"]
    y = data["y"]
    z = data["z"]
    m = data["mass"]
    r = data["radius"]
    # calculate inertia tensor given positions, mass, and radii of spheres
    # I_sphere = 2/5 * m r^2
    I = np.array([[0., 0., 0.],
                  [0., 0., 0.],
                  [0., 0., 0.]])
    for i in range(0, len(x)):
        a = np.array([x[i] - com[0], y[i] - com[1], z[i] - com[2]])
        a_mag = np.linalg.norm(a)
        I_sphere = (2./5)*m[i]*r[i]**2
        I_temp = np.array([[I_sphere, 0, 0],
                           [0, I_sphere, 0],
                           [0, 0, I_sphere]])
        for j in range(0, 3):
            for k in range(0,3):
                if j==k:
                    delta = 1.0
                else:
                    delta = 0.0
                I_temp[j, k] += m[i]*(delta*a_mag**2 - a[j]*a[k])
        I += I_temp 
    return I


################################################################################################
##
## 
## author       : Harrison Agrusa
## date         : 2/15/2018
## Description  : 
## Invocation   : 
## Input        : 
## Returns      : 
## Notes        : 
################################################################################################

# this is bad function
def calc_lib_angle(rp1, rp2, prev_vector=None):
    # calculate the libration angle of the moon
    # rp2 is nominally the moon, so be careful
    if prev_vector == None:
        print("You must specify the previous principle axis.")
        print("If this is the first angle you're calculating")
        print("Specify prev_principle_axis=False")
        quit()
    elif prev_vector is not np.ndarray:
        eig_val, eig_vec = np.linalg.eig(calc_I(rp2, calc_COM(rp2)))
        prev_vector = eig_vec[np.argmin(eig_val)]

    r_A = calc_COM(rp1)
    r_B = calc_COM(rp2)
    I = calc_I(rp2, r_A)
    eig_val, eig_vec = np.linalg.eig(I)
    principal_vector = eig_vec[np.argmin(eig_val)]# principal rotation axis has smallest eigenvalue
    if np.dot(principal_vector, prev_vector) < 0:
        principal_vector = -principal_vector
    
    #com vector is 'x' direction
    CoM_vec = r_B - r_A

    #get orbital angular momentum vector to define 'z' direction
    M_A = np.sum(rp1["mass"])
    M_B = np.sum(rp2["mass"])  
    v_A = calc_COM_vel(rp1)
    v_B = calc_COM_vel(rp2)
    L_orb = M_A * np.cross(r_A, v_A) + M_B * np.cross(r_B, v_B)

    #define coordinates
    z = L_orb/np.linalg.norm(L_orb)
    x = CoM_vec/np.linalg.norm(CoM_vec)
    y = np.cross(z, x)

    lat_lib = np.dot(principal_vector, z)
    long_lib = np.dot(principal_vector, y)
    #cosangle = np.dot(principal_vector, CoM_vector)
    return lat_lib, long_lib, principal_vector


################################################################################################
##
## calc_COM
## author       : Harrison Agrusa
## date         : 2/15/2018
## Description  : calculates center of mass of rubble pile
## Invocation   : 
## Input        : data dictionary
## Returns      : 3 element center of mass array
## Notes        : 

def calc_COM(data):
    # given data dictionary
    # return 3 element center of mass
    x = data["x"]
    y = data["y"]
    z = data["z"]
    m = data["mass"]
    x_avg = np.sum(x * m) / np.sum(m)
    y_avg = np.sum(y * m) / np.sum(m)
    z_avg = np.sum(z * m) / np.sum(m)
    CoM = np.array([x_avg, y_avg, z_avg])
    del x, y, z, m, x_avg, y_avg, z_avg
    return CoM
################################################################################################



################################################################################################
##
## calc_COM_vel
## author       : Harrison Agrusa
## date         : 2/15/2018
## Description  : calculates center of mass velocity vector
## Invocation   : 
## Input        : data dictionary
## Returns      : 3 element velocity vector of COM
## Notes        : 
################################################################################################
def calc_COM_vel(data):
    #calculate center of mass velocity
    #return 3 element vector
    x = data["x"]
    y = data["y"]
    z = data["z"]
    vx = data["xdot"]
    vy = data["ydot"]
    vz = data["zdot"]
    m = data["mass"]
    vx_avg = np.sum(vx*m)/np.sum(m)
    vy_avg = np.sum(vy*m)/np.sum(m)
    vz_avg = np.sum(vz*m)/np.sum(m)
    del x, y, z, m, vx, vy, vz
    return np.array([vx_avg, vy_avg, vz_avg])



################################################################################################
##
## 
## author       : Harrison Agrusa
## date         : 2/15/2018
## Description  : 
## Invocation   : 
## Input        : 
## Returns      : array containing kinetic energy of all particles
## Notes        : 
################################################################################################
def calc_Ekin(data):
    # given data dictionary
    xdot = data["xdot"]
    ydot = data["ydot"]
    zdot = data["zdot"]
    spin_x = data["spin_x"]
    spin_y = data["spin_y"]
    spin_z = data["spin_z"]
    mass = data["mass"]
    radius = data["radius"]

    v2 = np.square(xdot) + np.square(ydot) + np.square(zdot)
    w2 = np.square(spin_x) + np.square(spin_y) + np.square(spin_z)
    # Do kinetic energy
    E_k = 0.5*mass*v2
    # do rotational energy
    #E_rot = 0.5*(2./5)*mass * np.square(radius)*w2
    #E_tot = E_k + E_rot
    del xdot, ydot, zdot, spin_x, spin_y, spin_z, mass, radius
    del v2, w2, E_k, E_rot
    return E_k
    # return E_k + E_rot

################################################################################################
##
## 
## author       : Harrison Agrusa
## date         : 2/15/2018
## Description  : 
## Invocation   : 
## Input        : 
## Returns      : 
## Notes        : 
################################################################################################


def calc_Ebinding(data, ssID, G=1.0):
    x = data["x"]
    y = data["y"]
    z = data["z"]
    mass = data["mass"]
    # do gravitational potential
    # temp = np.array([mass, x, y, z])
    rand_ID = np.random.randint(1e5)
    fname = "xyzm_" + str(rand_ID) + str(ssID) + ".txt"
    outname = "E_bind_" + str(rand_ID) + str(ssID) + ".txt"
    f = open(fname, 'w')
    for i in range(0, len(x)):
        f.write(str(mass[i]) + ' ' +
                str(x[i]) + ' ' +
                str(y[i]) + ' ' +
                str(z[i]) + '\n')
    f.close()
    # np.savetxt(fname, temp)
    subprocess.call("./CalcBindingE {0} {1}".format(fname, outname), shell=True)
    E_grav = G*np.genfromtxt(outname)  # scalar
    subprocess.call("rm {0}".format(fname), shell=True)
    subprocess.call("rm {0}".format(outname), shell=True)
    return E_grav

################################################################################################
##
## 
## author       : Harrison Agrusa
## date         : 2/15/2018
## Description  : 
## Invocation   : 
## Input        : 
## Returns      : 
## Notes        : 
################################################################################################


def calc_MutualPotential(data, ssID, G=1.0):
    x = data["x"]
    y = data["y"]
    z = data["z"]
    mass = data["mass"]
    color = data["color"]
    rand_ID = np.random.randint(1e5) 
    # do gravitational potential
    # temp = np.array([mass, x, y, z])
    fname = "xyzmc_" + str(rand_ID) + str(ssID) + ".txt"
    outname = "E_pot_" + str(rand_ID) + str(ssID) + ".txt"
    f = open(fname, 'w')
    for i in range(0, len(x)):
        f.write(str(mass[i]) + ' ' +
                str(x[i]) + ' ' +
                str(y[i]) + ' ' +
                str(z[i]) + ' ' +
                str(int(color[i])) + '\n')
    f.close()
    # np.savetxt(fname, temp)
    subprocess.call("./CalcPotential {0} {1}".format(fname, outname), shell=True)
    E_grav = G*np.genfromtxt(outname)/2  # scalar, divide by two b/c of double counting - should fix in the C code
    subprocess.call("rm {0}".format(fname), shell=True)
    subprocess.call("rm {0}".format(outname), shell=True)
    return E_grav



################################################################################################
##
## 
## author       : Harrison Agrusa
## date         : 2/15/2018
## Description  : 
## Invocation   : 
## Input        : 
## Returns      : 
## Notes        : 
################################################################################################


def calc_Potential(data, G=1.0):
    x = data["x"]
    y = data["y"]
    z = data["z"]
    mass = data["mass"]
    color = data["color"]
    rand_ID = np.random.randint(1e5) 
    # do gravitational potential
    # temp = np.array([mass, x, y, z])
    fname = "xyzmc_" + str(rand_ID)  + ".txt"
    outname = "E_pot_" + str(rand_ID) + ".txt"
    f = open(fname, 'w')
    for i in range(0, len(x)):
        f.write(str(mass[i]) + ' ' +
                str(x[i]) + ' ' +
                str(y[i]) + ' ' +
                str(z[i]) + ' ' +
                str(int(color[i])) + '\n')
    f.close()
    # np.savetxt(fname, temp)
    subprocess.call("./CalcP {0} {1}".format(fname, outname), shell=True)
    E_grav = G*np.genfromtxt(outname)  # scalar, divide by two b/c of double counting - should fix in the C code
    subprocess.call("rm {0}".format(fname), shell=True)
    subprocess.call("rm {0}".format(outname), shell=True)
    return E_grav



################################################################################################
##
## 
## author       : Harrison Agrusa
## date         : 2/15/2018
## Description  : 
## Invocation   : 
## Input        : 
## Returns      : 
## Notes        : 
################################################################################################



def calc_Lorb(data, COM):

    ## Under progress - this technically only works for 2 spheres!!

    #calculat orbital angular momentum only
    #COM = np.array([0.0, 0.0, 0.0])
    # setting COM to zero for now
    # for this project, we want angular momentum W.R.T. the origin
    # we are ensuring that the COM starts at zero, and shouldnt have any drift
    x = data["x"]
    y = data["y"]
    z = data["z"]
    xdot = data["xdot"]
    ydot = data["ydot"]
    zdot = data["zdot"]
    spin_x = data["spin_x"]
    spin_y = data["spin_y"]
    spin_z = data["spin_z"]
    mass = data["mass"]
    radius = data["radius"]
    #w = np.array([spin_x, spin_y, spin_z])

    Lrot = np.array([0.0, 0.0, 0.0])

    # do rotational angular momentum
    for i in range(0, len(x)):
        r = np.array([x[i] - COM[0],
                      y[i] - COM[1],
                      z[i] - COM[2]])
        v = np.array([xdot[i], ydot[i], zdot[i]])
        current_Lrot = np.cross(r, mass[i] * v)
        Lrot += current_Lrot
    del x, y, z, xdot, ydot, zdot, spin_x, spin_y, spin_z, mass, radius
    return Lrot

################################################################################################
##
## 
## author       : Harrison Agrusa
## date         : 2/15/2018
## Description  : 
## Invocation   : 
## Input        : 
## Returns      : 
## Notes        : 
################################################################################################


def calc_Lvec(data, com):
    #COM = np.array([0.0, 0.0, 0.0])
    # setting COM to zero for now
    # for this project, we want angular momentum W.R.T. the origin
    # we are ensuring that the COM starts at zero, and shouldnt have any drift
    x = data["x"]
    y = data["y"]
    z = data["z"]
    xdot = data["xdot"]
    ydot = data["ydot"]
    zdot = data["zdot"]
    spin_x = data["spin_x"]
    spin_y = data["spin_y"]
    spin_z = data["spin_z"]
    m = data["mass"]
    R = data["radius"]
    #w = np.array([spin_x, spin_y, spin_z])
    com_v = calc_COM_vel(data)
    #Lrot = np.array([0.0, 0.0, 0.0])
    #Lspin = np.array([0.0, 0.0, 0.0])
    r = np.array([x - com[0], y - com[1], z - com[2]])
    #print('calculating L vector')
    #print(r)
    w = np.array([spin_x, spin_y, spin_z]) * 2*np.pi

    w = np.array([spin_x, spin_y, spin_z]) 

    
    #print(w)
    v = np.array([xdot - com_v[0], ydot - com_v[1], zdot - com_v[2]])
    Lvec = np.sum(m*np.transpose(np.cross(r, v, axisa = 0, axisb = 0)) + (2./5)*m*np.square(R) * w, axis=1)
    #print(Lvec)
    #print((2./5)*m*np.square(R) * w)
    # do rotational angular momentum
    """
    for i in range(0, len(x)):
        r = np.array([x[i] - COM[0],
                      y[i] - COM[1],
                      z[i] - COM[2]])
        #v = np.array([xdot[i], ydot[i], zdot[i]]) 
        v = np.array([xdot[i]-COM_v[0], ydot[i]-COM_v[1], zdot[i]-COM_v[2]])
        current_Lrot = np.cross(r, m[i] * v)
        Lrot += current_Lrot
    # do spin angular momentum
        
    for i in range(0, len(x)):
        a = np.array([x[i] - com[0], y[i] - com[1], z[i] - com[2]])
        a_mag = np.linalg.norm(a)
        I_sphere = (2./5)*m[i]*r[i]**2
        I = np.array([[I_sphere, 0, 0],
                           [0, I_sphere, 0],
                           [0, 0, I_sphere]])
        for j in range(0, 3):
            for k in range(0,3):
                if j==k:
                    delta = 1.0
                else:
                    dela = 0.0
                I[j, k] += m[i]*(delta*a_mag**2 - a[j]*a[k])
        w = np.array([spin_x[i], spin_y[i], spin_z[i]]) * 2*np.pi
        Lspin += np.dot(I, w)
    Lvec = Lrot + Lspin
    del I, I_xx, I_yy, I_zz, I_xy, I_yx, I_yz, I_zy, I_xz, I_zx
    del x, y, z, xdot, ydot, zdot, spin_x, spin_y, spin_z, mass, radius
    del w, Lrot, Lspin, current_Lrot
    """
    return Lvec

################################################################################################
##   
## orbital_elements
## author       : Harrison Agrusa
## date         : 2/15/2018
## Description  : calculate the orbital elements given the 2 body system in xyz coords
## Invocation   : 
## Input        : m1, m2, r1, r2, v1, v2 (where r & v are 3 element position/velocity vectors of the center of masses)
## Returns      : all returned angles are in degrees, a is in natural units, e is unitless
##              : a - semi-major axis
##              : e - eccentricity
##              : i - inclination
##              : OMEGA - longitude of the ascending node
##              : w - argument of pericenter
## Notes        : Assumes inputs are in a unit system such that G=1
################################################################################################

def orbital_elements(m1, m2, r1, r2, v1, v2):
    #r = r1 - r2
    r = r2-r1
    rmag = np.linalg.norm(r)
    rhat = r/rmag

    v = v2 - v1
    vmag = np.linalg.norm(v)
    mu = (m1+m2)

    h = np.cross(r, v)

    OMEGA = np.arctan(- h[0]/h[1])

    #i = np.arctan(np.sqrt(h[0]**2 + h[1]**2)/h[2])
    i = np.arccos(h[2]/np.linalg.norm(h))
    
    e_vec = np.cross(v, h)/mu - rhat
    e = np.linalg.norm(e_vec)

    w = np.arctan( (e_vec[2]/np.sin(i)) / ((e_vec[0]*np.cos(OMEGA) + e_vec[1]*np.sin(OMEGA))))

    a_inv = 2/rmag - (1./mu)*vmag**2
    a = a_inv**-1

    return np.array([a, e, np.degrees(i), np.degrees(OMEGA), np.degrees(w)])




################################################################################################
##
## edit_sspar
## author       : Harrison Agrusa
## date         : 2/15/2018
## Description  : change parameters for ss.par file
## Invocation   : 
## Input        : path to ss.par and whatever parameters you want to change those not specified are not changed 
## Input options: dDelta, nSteps, dTheta, achInFile, achOutName, iOutInterval, iRedOutInterval, iLogInterval, iCheckInterval
## Returns      : returns the paramters you have changed in a dictionary, and rewrites the given ss.par file with new paramters
## Notes        : 
## WARNING #1   : THIS REWRITES THE SS.PAR FILE THAT YOU GIVE IT, BE CAREFUL!
## WARNING #2   : THIS WILL DELETE ANY COMMENTS ON THE PARAMETERS THAT YOU CHANGE
################################################################################################

def edit_sspar(sspar, dDelta=None, dTheta=None, nSteps=None,
               achInFile=None, achOutName=None, iOutInterval=None,
               iRedOutInterval=None, iLogInterval=None, iCheckInterval=None,
               bReadDEMData=None, iStartStep=None):
    params = {"dDelta": dDelta, "dTheta": dTheta, 
            "nSteps": nSteps,  "achInFile": achInFile, 
            "achOutName": achOutName, "iOutInterval":iOutInterval,
            "iRedOutInterval": iRedOutInterval, "iLogInterval": iLogInterval,
            "iCheckInterval":iCheckInterval, "bReadDEMData": bReadDEMData, 
            "iStartStep": iStartStep}
    done = 0
    file = open(sspar)
    new_name = "temp_"+sspar
    new_file = open(new_name, 'w')
    for line in file:
        for key in params.keys():
            if line.startswith(key) and params[key] != None:
                line = key + "  =  " + str(params[key]) + "          #this line changed with dart_utilities.edit_sspar() " + "\n"#caution: any comments on this line will be lost
                done += 1
        new_file.write(line)
    new_file.close()
    subprocess.call("mv {0} {1}".format(new_name, sspar), shell=True) #replace old ss file
    assert done > 0
    return params


################################################################################################
##
## get_sspar
## author       : Harrison Agrusa
## date         : 2/15/2018
## Description  : get specified quantities from ss.par file, returned in dictionary.
## Invocation   : 
## Input        : path to ss.par and whatever parameters you want 
## Input options: dDelta, nSteps, dTheta, achInFile, achOutName, iOutInterval, iRedOutInterval, iLogInterval, iCheckInterval
## Returns      : returns the paramters you want in a dictionary
## Notes        : 
## WARNING      : be careful with returned data types, I think everything will just be strings so youll have to conver afterwords
################################################################################################

def get_sspar(sspar, dDelta=None, dTheta=None, nSteps=None,
               achInFile=None, achOutName=None, iOutInterval=None,
               iRedOutInterval=None, iLogInterval=None, iCheckInterval=None,
               bReadDEMData=None, iStartStep=None):

    possible_params = {"dDelta": dDelta, "dTheta": dTheta, 
            "nSteps": nSteps,  "achInFile": achInFile, 
            "achOutName": achOutName, "iOutInterval":iOutInterval,
            "iRedOutInterval": iRedOutInterval, "iLogInterval": iLogInterval,
            "iCheckInterval":iCheckInterval, "bReadDEMData": bReadDEMData, 
            "iStartStep": iStartStep}
    print(possible_params)
    params = {} #build dictionary to return
    done = 0
    file = open(sspar)
    for line in file: #iterate through lines
        for key in possible_params.keys():
            if line.startswith(key) and possible_params[key] != None:
                print(line.split()[2])
                params[key] = line.split()[2]
                done += 1
    assert done > 0
    return params

################################################################################################
##
## fit_sin
## author       : not me, fit_sin is copy pasted code from stack exchange - cannot remember where exactly but it works well
## date         : 11/15
## Description  : 
## Invocation   : 
## Input        : 
## Returns      : 
## Notes        : 
################################################################################################

def fit_sin(tt, yy):
    '''Fit sin to the input time sequence, and return fitting parameters "amp", "omega", "phase", "offset", "freq", "period" and "fitfunc"'''
    #tt = np.array(tt)
    #yy = np.array(yy)
    ff = np.fft.fftfreq(len(tt), (tt[1]-tt[0]))   # assume uniform spacing
    Fyy = abs(np.fft.fft(yy))
    guess_freq = abs(ff[np.argmax(Fyy[1:])+1])   # excluding the zero frequency "peak", which is related to offset
    guess_amp = np.std(yy) * 2.**0.5
    guess_offset = np.mean(yy)
    guess = np.array([guess_amp, 2.*np.pi*guess_freq, 0., guess_offset])

    def sinfunc(t, A, w, p, c):  return A * np.sin(w*t + p) + c
    popt, pcov = scipy.optimize.curve_fit(sinfunc, tt, yy, p0=guess)
    A, w, p, c = popt
    f = w/(2.*np.pi)
    fitfunc = lambda t: A * np.sin(w*t + p) + c
    return {"amp": A, "omega": w, "phase": p, "offset": c, "freq": f, "period": 1./f, "fitfunc": fitfunc, "maxcov": np.max(pcov), "rawres": (guess,popt,pcov)}



################################################################################################
##
## fit_period
## author       : Harrison Agrusa
## date         : 11/15
## Description  : take file in my 'allstats' format, take x position of moon 
## Invocation   : 
## Input        : 
## Returns      : 
## Notes        : 
################################################################################################

def fit_period(path_to_allstats_file):
    keys = ["time", "M1", "M2", "E_kin", "E_pot", "E_bind",
       "Ltot", "Lx", "Ly", "Lz",
        "CoMx", "CoMy", "CoMz",
        "CoM1x", "CoM1y", "CoM1z",
        "CoM2x", "CoM2y", "CoM2z",
        "CoM1vx", "CoM1vy", "CoM1vz",
        "CoM2vx", "CoM2vy", "CoM2vz", 
        "r", "a", "e", "i", "OMEGA", "w", "obliq1", "obliq2"]
    data = {}
    data_array = np.genfromtxt(path_to_allstats_file)
    for j in range(0, len(keys)): #build dictionary
        data[keys[j]] = data_array[:, j]
    t = data['time']
    x = data['CoM2x']
    y = data['CoM2y']
    z = data['CoM2z']

    fit = fit_sin(t, x) #fit x vs t to sin wave
    x_period = fit['period']
    #fit = fit_sin(t, y) #fit y vs t to sin wave
    #y_period = fit['period']
    #fit = fit_sin(t, z) #fit z vs t to sin wave
    #z_period = fit['period']
    return x_period

################################################################################################
##
## build_allstats
## author       : Harrison Agrusa
## date         : 11/21/2018
## Description  :
## Invocation   : 
## Input        : 
## Returns      : 
## Notes        : 
################################################################################################

def build_allstats(path_to_allstats_file, tmax = False):
    keys = ["time", "M1", "M2", "E_kin", "E_pot", "E_bind",
       "Ltot", "Lx", "Ly", "Lz",
        "CoMx", "CoMy", "CoMz",
        "CoM1x", "CoM1y", "CoM1z",
        "CoM2x", "CoM2y", "CoM2z",
        "CoM1vx", "CoM1vy", "CoM1vz",
        "CoM2vx", "CoM2vy", "CoM2vz", 
        "r", "a", "e", "i", "OMEGA", "w", "obliq1", "obliq2"]
    data = {}
    data_array = np.genfromtxt(path_to_allstats_file)
    if tmax:
        t = data_array[:, 0]
        indices = np.where(t < tmax) 
        for j in range(0, len(keys)): #build dictionary
            data[keys[j]] = data_array[:, j][indices]
    else:
        for j in range(0, len(keys)): #build dictionary
            data[keys[j]] = data_array[:, j]
    return data











