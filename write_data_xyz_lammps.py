#this script is hard coded for writing lammps data file from the xyz file 
#pay attention to the water bond sequence while changing this script
#need to follow the right sequence of Hw and Ow in xyz file

import numpy as np
import glob


filename2 = 'Hamid_tober_8103889.xyz' # xyz file was given gy Roland E Pellenq 
flist2 = glob.glob(filename2)

#the 1st col. is string so its loaded seperately

for f in flist2:
    load2 = np.genfromtxt(f, skip_header=2, dtype=float, usecols=(1,2,3)) #dtype=("|S10", float, float, float),
    dataovito1=np.array(load2)

for f in flist2:
    load2 = np.genfromtxt(f, skip_header=2, dtype=str, usecols=(0)) #dtype=("|S10", float, float, float),
    dataovito2=np.array(load2)

nSi = 0
nCa = 0
nO = 0


for j in range(len(dataovito2)):

  if dataovito2[j] == 'Si':
            nSi +=1
  elif dataovito2[j] =='Ca':
            nCa +=1
  elif dataovito2[j] =='O':
            nO +=1


#totCharge =  2.1*nSi + 1.7*nCa - 1.14*nO
totCharge =  2.0*nSi + 1.5*nCa - 1.05*nO
print "Total charge= ", totCharge


size1 = len(dataovito1) #total atoms-water
size2 = len(dataovito1) #total atoms

#print dataovito2

#xy = -16.26
#xz = 0.5433
#yz = -0.4834

#xlo_bound = xlo + MIN(0.0,xy,xz,xy+xz)
#xhi_bound = xhi + MAX(0.0,xy,xz,xy+xz)
#ylo_bound = ylo + MIN(0.0,yz)
#yhi_bound = yhi + MAX(0.0,yz)
#zlo_bound = zlo
#zhi_bound = zhi

xmin = np.min(dataovito1[:,0]) #+min(0.0,xy,xz,xy+xz)
xmax = np.max(dataovito1[:,0]) #+max(0.0,xy,xz,xy+xz)

ymin = np.min(dataovito1[:,1]) #+min(0.0,yz)
ymax = np.max(dataovito1[:,1]) #+max(0.0,yz)

zmin = np.min(dataovito1[:,2])
zmax = np.max(dataovito1[:,2])

a = 6.69  
b = 7.39  
c = 22.77

alpha = 90.0 
beta  =  90.0  
gamma = 123.49


lx=a
xy = b*np.cos(gamma*np.pi/180)
ly=np.sqrt(b**2 - xy**2)
xz = c*np.cos(beta*np.pi/180)
yz = (b*c*np.cos(alpha*np.pi/180)-xy*xz)/ly
lz = np.sqrt(c**2 - xz**2 - yz**2)

outFile = open('data.Tobermorite_Hamid-11A', 'w')
outFile.write('LAMMPS data file written by Rafat Sadat using Python\n')
outFile.write('\n')
outFile.write('%i %s \n' %(len(dataovito1), 'atoms'))
outFile.write('6 atom types \n')
#outFile.write('%i %s \n' %(nbonds, 'bonds'))
#outFile.write('1 bond types \n')
#outFile.write('%i %s \n' %(nangles, 'angles'))
#outFile.write('1 angle types \n')
outFile.write('\n')
outFile.write('%f %f %s %s \n' %(0.0, lx, 'xlo', 'xhi'))
outFile.write('%f %f %s %s \n' %(0.0, ly, 'ylo', 'yhi'))
outFile.write('%f %f %s %s \n' %(0.0, lz, 'zlo', 'zhi'))
outFile.write('\n')
outFile.write('%f %f %f %s %s %s \n' %(xy, xz, yz, 'xy', 'xz', 'yz'))
outFile.write('\n')
outFile.write('Atoms\n')
outFile.write('\n')


for j in range(size1): #writing atoms without water
    if dataovito2[j]=="Si":
        outFile.write('%i %i %i %i %f %f %f \n' %(j+1, 0, 1, 0, dataovito1[j,0], dataovito1[j,1], dataovito1[j,2]))
    elif dataovito2[j]=="O":
        outFile.write('%i %i %i %i %f %f %f \n' %(j+1, 0, 4, 0, dataovito1[j,0], dataovito1[j,1], dataovito1[j,2]))
    elif dataovito2[j]=="Ca":
        outFile.write('%i %i %i %i %f %f %f \n' %(j+1, 0, 3, 0, dataovito1[j,0], dataovito1[j,1], dataovito1[j,2]))
    elif dataovito2[j]=='Xx': #water
        outFile.write('%i %i %i %i %f %f %f \n' %(j+1, 0, 6, 0, dataovito1[j,0], dataovito1[j,1], dataovito1[j,2])) 


'''
counter = -1
molid = 0
for j in range(size1,size2): #writing the atoms for water molecules
    counter +=1
    if counter%3==0:
     molid += 1 
    else:
        molid = molid 
         
    if dataovito2[j]=='Ow':
        outFile.write('%i %i %i %i %f %f %f \n' %(j+1, molid, 5, 0, dataovito1[j,0], dataovito1[j,1], dataovito1[j,2]))
    elif dataovito2[j]=='Hw':
        outFile.write('%i %i %i %i %f %f %f \n' %(j+1, molid, 6, 0, dataovito1[j,0], dataovito1[j,1], dataovito1[j,2])) 
        
              
outFile.write('\n')
outFile.write('Bonds\n')
outFile.write('\n')

count = 0
for j in range(size1+1,size2, 3):
    count +=1
    outFile.write('%i %s %i %i \n' %(count,'1',j,j+1))
    count +=1
    outFile.write('%i %s %i %i  \n' %(count,'1',j,j+2))

outFile.write('\n')
outFile.write('Angles\n')
outFile.write('\n')   

count = 0
for j in range(size1+1,size2, 3):
    count +=1
    outFile.write('%i %s %i %i %i \n' %(count,'1', j+1, j, j+2))

'''
outFile.close()
print "All done!"                          

