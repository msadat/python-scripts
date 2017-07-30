import matplotlib.pyplot as plt
import numpy as np
import glob


#delete the headers every timesteps (snapshots)
'''
sed -i '/ITEM: TIMESTEP/,+8d' dump.test

'''

natoms = 86718   #2023 21054 9126
nsnaps = 31 #no. of snaps


filename = 'dump.test' # lammps dump file after processing with 'sed' command
flist = glob.glob(filename)


for f in flist:
    load = np.genfromtxt(f, dtype=float, usecols=(1,2,3,4)) #atomtype, x, y, z
    data=np.array(load)

data_split = np.zeros((nsnaps, natoms, 4))  

data_split [:] = np.split(data,nsnaps) 


# check data_splits again by printing it

print data_split [0,1,1] #snap, row, column


#writing the New Dumpfile with only atoms needed

natoms_new = 35971 #number of atoms after deletion

xlo = -6.684e-2 #-2.061e-1 #4.2366e-1 #1.5297e-1
xhi = 1.7819e2 #1.787e2 #1.77316e2 #1.765e2
ylo = 4.284e-1 #5.0795e-1 #5.7174e-1 #8.7306e-2
yhi = 8.919599e1 #8.904e1 #8.877e1 #8.82078e1
zlo = 4.1064e-1 #6.53536e-1 #1.3324855 #1.426e-1
zhi = 8.92e1 #8.888e1 #8.8e1 #8.7307e1

outFile = open('dump.GPCSH2.0_deleted', 'w')

for s in range(nsnaps):

 #outFile = open('dump.GPCSH1.2_deleted', 'a')
 outFile.write('ITEM: TIMESTEP \n')
 outFile.write('%i \n' %(s))
 outFile.write('ITEM: NUMBER OF ATOMS \n')
 outFile.write('%i \n' %(natoms_new))
 outFile.write('ITEM: BOX BOUNDS pp pp pp \n')
 outFile.write('%f %f \n' %(xlo-0.005*s, xhi+0.005*s))
 outFile.write('%f %f \n' %(ylo, yhi))
 outFile.write('%f %f \n' %(zlo, zhi))
 outFile.write('ITEM: ATOMS id type xs ys zs \n')


 for i in range(natoms):

    if int(data_split [s,i,0]) == 1:
        outFile.write('%i %i %f %f %f \n' %(i+1, data_split [s,i,0], data_split [s,i,1], data_split [s,i,2], data_split [s,i,3]))
    elif int(data_split [s,i,0]) == 2:
        outFile.write('%i %i %f %f %f \n' %(i+1, data_split [s,i,0], data_split [s,i,1], data_split [s,i,2], data_split [s,i,3]))   
    elif int(data_split [s,i,0]) == 10:
        outFile.write('%i %i %f %f %f \n' %(i+1, 10, data_split [s,i,1], data_split [s,i,2], data_split [s,i,3]))  
    elif int(data_split [s,i,0]) == 11:
        outFile.write('%i %i %f %f %f \n' %(i+1, 11, data_split [s,i,1], data_split [s,i,2], data_split [s,i,3]))                 

outFile.close()
print "All done!" 
