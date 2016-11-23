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
'''
sed -i '/ITEM: TIMESTEP/,+8d' dump.test

'''

natoms = 2023
nsnaps = 61 #no. of snaps
nH2O = 1200

filename = 'dump.test' # lammps dump file after processing with 'sed' command
flist = glob.glob(filename)


for f in flist:
    load = np.genfromtxt(f, dtype=float, usecols=(1,7,8)) #atomtype, fz, pe
    data=np.array(load)

data_split = np.zeros((nsnaps, 2023, 3))  
#print np.shape(data)

data_split [:] = np.split(data,nsnaps) 

fz_array = np.zeros((nsnaps-1,2))
pe_array = np.zeros((nsnaps-1,2))

for i in range(nsnaps-1):
    
    fz_array [i,0] = i
    fz_array [i,1] = np.max(data_split [i,:,1]) 


plt.plot(fz_array[:,0],fz_array[:,1])
plt.xlabel('Time (ps)')
plt.ylabel('force (eV/A)')
plt.show()



#np.savetxt('H2O_GP3.0_MSD',np.c_[data_avg[:,0],data_avg[:,4]])
