## this script is a post process for lammps nanoindent *****
# the motivation for making this script is to find the maximum per atom force and per atom pe of run from dumpfiles
# nano indent output gives the total force felt by a indenter, so bigger the indenter the more force (sum of all atoms below the indenter) will be felt...
# ...however, the peratom force is not seen from these output 

# this script will take in a dump file containing of potential energy and force with multiple timesteps and get the maximum and minimum as a function of timesteps
# output1: maximum and minimum potential energy as a function of timesteps
# output2: maximum and minimum forces as a function of timesteps

import matplotlib.pyplot as plt
import numpy as np
import glob


#delete the headers every timesteps (snapshots)
#the following command deletes the line with the string "ITEM:..*" and also 8 lines following that string
'''
sed -i '/ITEM: TIMESTEP/,+8d' dump.test 

'''

natoms = 1646 #2023 21054 9126
nsnaps = 601 #no. of snaps or timesteps (*actual timestep in lammps to get the real timesteps); note: 0 is also a step


filename = 'dump.test' # lammps dump file after processing with 'sed' command
flist = glob.glob(filename)


for f in flist:
    load = np.genfromtxt(f, dtype=float, usecols=(4,7)) #z-coord, fz
    data=np.array(load)

#define a 3-D array for storing snapshot, z-coordinate and fz

#each storey is for each timesteps or snapshot (you can get actual timestep from snaps; timesteps = snap * lammps timestep)
#each row is for a certain atom
#each col1 is the z coord
#each col2 is the fz

data_split = np.zeros((nsnaps, natoms, 2)) 

data_split [:] = np.split(data,nsnaps) #equally divides the data into 'nsnaps' equal size along the default axis

# check data_splits again by printing it

# print data_split [1,:,:]
#make array for storing the force and y-position along with timesteps
fz_array = np.zeros((nsnaps-1,2))
z_array = np.zeros((nsnaps-1,2))

for i in range(nsnaps-1):
    
    fz_array [i,0] = i #timesteps
    fz_array [i,1] = np.mean(data_split [i,:,1]) #instead of just taking the max you might take average force too

for i in range(nsnaps-1):
    
    z_array [i,0] = i #timesteps
    z_array [i,1] = np.mean(data_split [i,:,0]) 


#plt.plot(fz_array[:,0],fz_array[:,1])
#plt.xlabel('Time (ps)')
#plt.ylabel('force (eV/A)')

plt.plot(fz_array[:,0],fz_array[:,1])
plt.xlabel('Time (ps)')
plt.ylabel('force_z (eV/A))')

plt.show()



#np.savetxt('H2O_GP3.0_MSD',np.c_[data_avg[:,0],data_avg[:,4]])
