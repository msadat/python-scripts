import numpy as np
import glob

filename2 = 'data.poly-d200' #lammps 'charge' data file
outputfile = 'data.GP-composite_200nm_charge'
flist2 = glob.glob(filename2)

#the 1st col. is string so its loaded seperately

for f in flist2:
    load2 = np.genfromtxt(f, dtype=float, skip_header=10, usecols=(1,3,4,5)) #dtype=("|S10", float, float, float),
    dataovito1=np.array(load2)

for f in flist2:
    load2 = np.genfromtxt(f, dtype=str, skip_header=10, usecols=(1)) #dtype=("|S10", float, float, float),
    dataovito2=np.array(load2)


size1 = len(dataovito1) #total atoms-water
size2 = len(dataovito1) #total atoms

grain_types = 131

grain_list = []
grain_new = []

for i in range(1,grain_types+1):
    grain_list.append(i)

for i in range(1,grain_types+1):
    grain_new.append(int(np.random.uniform(1,4)))

natoms = size2

#print dataovito2

xmin = 0.0 #np.min(dataovito1[:,0]) # added -1.0 to prohibit overlapping at boundaries
xmax =  1000 #np.max(dataovito1[:,0])

ymin = 0.0 #np.min(dataovito1[:,1])
ymax = 1000 #np.max(dataovito1[:,1])

zmin = 0.0 #np.min(dataovito1[:,2])
zmax = 1000 #np.max(dataovito1[:,2])

#print xmax, ymax, zmax

density = 0.0022 #2.2 gm/cm3 or 0.0022 atto-gram/nm3
volume = 1000 #cube of sc

outFile = open(outputfile, 'w')
outFile.write('PDLAMMPS data file written by Python\n')
outFile.write('\n')
outFile.write('%i %s \n' %(natoms, 'atoms'))
outFile.write('3 atom types \n')
outFile.write('\n')
outFile.write('%f %f %s %s \n' %(xmin-1, xmax+1, 'xlo', 'xhi'))
outFile.write('%f %f %s %s \n' %(ymin-1, ymax+1, 'ylo', 'yhi'))
outFile.write('%f %f %s %s \n' %(zmin-1, zmax+1, 'zlo', 'zhi'))
outFile.write('\n')
outFile.write('Atoms # charge\n')
outFile.write('\n')

for j in range(size1):
    
    #if int(dataovito2[j]) in grain_list:
     #outFile.write('%i %i %i %f %f %f \n' %(j+1, np.random.uniform(1,4), 1, dataovito1[j,1], dataovito1[j,2], dataovito1[j,3]))
    outFile.write('%i %i %i %f %f %f \n' %(j+1, grain_new[grain_list.index(int(dataovito2[j]))], 1, dataovito1[j,1], dataovito1[j,2], dataovito1[j,3]))
    '''     
    
    elif dataovito2[j]=="2":
        outFile.write('%i %i %s %s %s %s %s %i %i %i \n' %(j+1, np.random.uniform(1,5), volume, density, dataovito1[j,0], dataovito1[j,1], dataovito1[j,2], 0, 0, 0))
    elif dataovito2[j]=="3":
        outFile.write('%i %i %s %s %s %s %s %i %i %i \n' %(j+1, np.random.uniform(1,5), volume, density, dataovito1[j,0], dataovito1[j,1], dataovito1[j,2], 0, 0, 0))
    elif dataovito2[j]=="4":
        outFile.write('%i %i %s %s %s %s %s %i %i %i \n' %(j+1, np.random.uniform(1,5), volume, density, dataovito1[j,0], dataovito1[j,1], dataovito1[j,2], 0, 0, 0))
    elif dataovito2[j]=="5":
        outFile.write('%i %i %s %s %s %s %s %i %i %i \n' %(j+1, np.random.uniform(1,5), volume, density, dataovito1[j,0], dataovito1[j,1], dataovito1[j,2], 0, 0, 0))
    elif dataovito2[j]=="6":
        outFile.write('%i %i %s %s %s %s %s %i %i %i \n' %(j+1, np.random.uniform(1,5), volume, density, dataovito1[j,0], dataovito1[j,1], dataovito1[j,2], 0, 0, 0))   
    elif dataovito2[j]=="7":
        outFile.write('%i %i %s %s %s %s %s %i %i %i \n' %(j+1, np.random.uniform(1,5), volume, density, dataovito1[j,0], dataovito1[j,1], dataovito1[j,2], 0, 0, 0))
    elif dataovito2[j]=="8":
        outFile.write('%i %i %s %s %s %s %s %i %i %i \n' %(j+1, np.random.uniform(1,5), volume, density, dataovito1[j,0], dataovito1[j,1], dataovito1[j,2], 0, 0, 0))
         
    elif dataovito2[j]=="9":
        outFile.write('%i %i %s %s %s %s %s %i %i %i \n' %(j+1, 3, volume, density, dataovito1[j,0], dataovito1[j,1], dataovito1[j,2], 0, 0, 0))
    elif dataovito2[j]=="10":
        outFile.write('%i %i %s %s %s %s %s %i %i %i \n' %(j+1, 3, volume, density, dataovito1[j,0], dataovito1[j,1], dataovito1[j,2], 0, 0, 0))
    '''      
    
outFile.close()
print "All done!"                          

