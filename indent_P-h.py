
#this script plots and writes the load (nN) vs. penetration depth (nm) data from lammps output
#also outputs the H and E of the material 

import numpy as np
import glob
import matplotlib.pyplot as plt

filename = 'indent-GP2-Al2O3-15-1-10.dat' # lammps data file
flist = glob.glob(filename)
outfile = 'GP2-Al2O3_100.dat'

radius = 15 #A
area = 2*np.pi*radius**2 #A**2
K = 10 #ev/A3
Ei = K * 1.60218E-19 *1e-9 / (1e-10)**3   #indenter modulus, GPa 
vi = 0.2 #indenter poisson's ratio
v = 0.25 #susbsrate poisson's ratio


for f in flist:
    load = np.genfromtxt(f, dtype=float, skip_header =2, usecols=(0,1,2))
    data=np.array(load)

final_data = np.zeros((len(data),2))

rate = 1 # A/fs
initial_time = 100000 #need the initial timestep
for i in range(len(data)):
   if data[i,0]-initial_time <= 30000: 
   #if data[i,1] >= data[i-1,1]:
        final_data[i,0] = rate*0.001*(data[i,0]-initial_time)
   else:
        final_data[i,0] = final_data[i-1,0]-rate*0.001*(data[i,0]-data[i-1,0])

for i in range(len(data)):
        
        final_data[i,1] = data[i,2]*1.60218E-19*1e9/1e-10    #taking only the z component of the indenter force   

print "Hardness, H = ", np.max(final_data[:,1])*1e-18/(area*(1e-10)**2), " GPa"

#find slope of unloading curve
slope = (final_data[len(final_data)/2,1] - final_data[(len(final_data)/2)+10,1]) / (final_data[(len(final_data)/2)+10,0] - final_data[len(final_data)/2,0]) #nN/A
E_eff = slope*1e-18 / (2*(2**0.5)*radius*(1e-10)**2)  #GPa
E = (1-v**2)*Ei*E_eff / (E_eff*(1-vi**2)-Ei)

print "Young's modulus, E = ", E , " GPa"
                                       
plt.plot(final_data[:,0],final_data[:,1])
plt.xlabel('Depth (A)')
plt.ylabel('Force (nN)')
plt.show()

np.savetxt(outfile, np.c_[final_data[:,0],final_data[:,1]])
