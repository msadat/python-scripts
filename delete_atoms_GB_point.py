#Author: Mohammad R. Sadat
#PhD candidate
#Dept of CEEM, U of A
#this script randomly deletes atoms from a grain boundary while maintaining stoichiometry and writes data in 'full' style format after deletion
    
def delAtoms(x,t):
    """x = array of xyz data file with the atom types,
    t = list of atoms to be deleted """
    
    '''
    idx = []
    
    for i in range(len(x)):
        if x[i,0] == t:
            idx.append(i)
        else:
            idx = idx
    
    #print "no of atoms ="+ str(len(idx))
    #print "Initial points = "+str(len(x))
    #print "deleted atoms = "+str(len(idx))'''
    return np.delete(x,t,0) 
    
import numpy as np
import glob
import random

natoms_initial = int(raw_input("Enter the original no. of atoms: ")) #int(sys.argv[1])
filename = 'data.GB-test' # name of the input data file
outputfile = 'data.GB-test-defected' # name of output data file

flist2 = glob.glob(filename)

#the 1st col. is string so its loaded seperately

for f in flist2:
    load2 = np.genfromtxt(f, skip_header=16, skip_footer=natoms_initial+1, dtype=float, usecols=(1,2,3,4)) #dtype=("|S10", float, float, float),
    dataovito1=np.array(load2)

for f in flist2:
    load2 = np.genfromtxt(f, skip_header=16, skip_footer=natoms_initial+1, dtype=int, usecols=(1)) #dtype=("|S10", float, float, float),
    dataovito2=np.array(load2)

nZn = 0
nS = 0

for j in range(len(dataovito2)):

  if dataovito2[j] == 1:
            nZn +=1
  elif dataovito2[j] ==2:
            nS +=1

total = nZn+nS
print "total of type 1: ", nZn
print "total of type 2: ", nS
#to_delete = total*1/100

#extent for deletion (around GB plane)


xmin = np.min(dataovito1[:,1]) -1.0 
xmax = np.max(dataovito1[:,1]) 

ymin = np.min(dataovito1[:,2]) -1.0 
ymax = np.max(dataovito1[:,2]) 

zmin = np.min(dataovito1[:,3]) -1.0
zmax = np.max(dataovito1[:,3])

ylo = -2 
yhi = 2 

ndel = total*1/100

ids= [] #final list of atoms to be deleted
idsin1 = [] #list of type 1 atom within chosen boundary
idsin2 = [] #list of type 2 atom within chosen boundary

for j in range(len(dataovito1)):
    #if yhi < dataovito1[j,2] or dataovito1[j,2] < ylo and dataovito1[j,0]==atomtype: 
    if ylo < dataovito1[j,2] < yhi and dataovito1[j,0]==1:
        idsin1.append(j)
                
for j in range(len(dataovito1)):
    #if yhi < dataovito1[j,2] or dataovito1[j,2] < ylo and dataovito1[j,0]==atomtype: 
    if ylo < dataovito1[j,2] < yhi and dataovito1[j,0]==2:
        idsin2.append(j)

#print idsin1
#print idsin2


for j in range(len(idsin1)):

    ids.append(random.choice(idsin1))
    if len(ids) >=ndel/2:
        break

#print ids
                            
for j in range(len(idsin2)):

    ids.append(random.choice(idsin2))
    if len(ids)/2 >=ndel/2:
        break
                    
#print ids                  
                    

dataovito1= delAtoms(dataovito1, ids)
print "total deleted: ", len(ids)
#print total
print "final no of atoms: ", len(dataovito1)
#print  dataovito1 
natoms = len(dataovito1)

outFile = open(outputfile, 'w')
outFile.write('LAMMPS data file written using Python script\n')
outFile.write('\n')
outFile.write('%i %s \n' %(natoms, 'atoms'))
outFile.write('2 atom types \n')
outFile.write('\n')
outFile.write('%f %f %s %s \n' %(xmin, xmax, 'xlo', 'xhi'))
outFile.write('%f %f %s %s \n' %(ymin, ymax, 'ylo', 'yhi'))
outFile.write('%f %f %s %s \n' %(zmin, zmax, 'zlo', 'zhi'))
outFile.write('\n')
outFile.write('%s \n' %('Masses'))
outFile.write('\n')
outFile.write('%i %f \n' %(1, 65.39))
outFile.write('%i %f \n' %(2, 32.066))
outFile.write('\n')
outFile.write('Atoms\n')
outFile.write('\n')


for j in range(len(dataovito1)):
    if dataovito1[j,0]==1:
        outFile.write('%i %i %i %i %f %f %f \n' %(j+1, 0, 1, 0, dataovito1[j,1], dataovito1[j,2], dataovito1[j,3]))
    elif dataovito1[j,0]==2:
        outFile.write('%i %i %i %i %f %f %f \n' %(j+1, 0, 2, 0, dataovito1[j,1], dataovito1[j,2], dataovito1[j,3]))

outFile.close()
print "All done!"   
