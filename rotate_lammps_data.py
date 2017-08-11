import numpy as np
import glob


filename2 =  'data.GPCSH1.65r1-Eq-300K' #data.GP2.0w0r2-Compacted-Eq' # xyz file was given gy Roland E Pellenq 
flist2 = glob.glob(filename2)
natoms_GBP = 103231 

rotation = (90.00*np.pi)/180.00
c=np.cos(rotation)
s=np.sin(rotation)

#the 1st col. is string so its loaded seperately

for f in flist2:
    #load2 = np.genfromtxt(f, skip_header=12, dtype=float, usecols=(0,1, 2, 3, 4, 5, 6)) #dtype=("|S10", float, float, float),
    load2 = np.genfromtxt(f, skip_header=27, skip_footer=natoms_GBP+1, dtype=float, usecols=(0,1, 2, 3, 4, 5, 6)) #dtype=("|S10", float, float, float),

    dataovito1=np.array(load2)

for f in flist2:
    #load2 = np.genfromtxt(f, skip_header=12, dtype=str, usecols=(0)) #dtype=("|S10", float, float, float),
    load2 = np.genfromtxt(f, skip_header=27, skip_footer=natoms_GBP+1, dtype=str, usecols=(0)) #dtype=("|S10", float, float, float),

    dataovito2=np.array(load2)
    
datazero = np.zeros((len(dataovito1),5))

for j in range(len(dataovito1)):   #rotating wrt y axis
    
    datazero[j,0] = dataovito1[j,0]
    datazero[j,1] = dataovito1[j,2]
    datazero[j,2] = dataovito1[j,4]*c + dataovito1[j,6]*s
    datazero[j,3] = dataovito1[j,5] 
    datazero[j,4] = -dataovito1[j,4]*s + dataovito1[j,6]*c
    
    #dataovito1[j,4] = dataovito1[j,4]*c + dataovito1[j,6]*s
    

'''
for j in range(len(dataovito1)):   #rotating wrt x axis
    dataovito1[j,4] = -dataovito1[j,5]
    dataovito1[j,5] = dataovito1[j,4] 
    dataovito1[j,6] = dataovito1[j,6]
'''
size1 = len(dataovito1) #total atoms 
size2 = len(dataovito1) #total atoms

xmin = np.min(datazero[:,2]) #+min(0.0,xy,xz,xy+xz)
xmax = np.max(datazero[:,2]) #+max(0.0,xy,xz,xy+xz)

ymin = np.min(datazero[:,3]) #+min(0.0,yz)
ymax = np.max(datazero[:,3]) #+max(0.0,yz)

zmin = np.min(datazero[:,4])
zmax = np.max(datazero[:,4])

outFile = open('data.GPCSH1.65r1-Eq-300K-rotated', 'w')
outFile.write('LAMMPS data file written by Rafat Sadat using Python\n')
outFile.write('\n')
outFile.write('%i %s \n' %(len(dataovito1), 'atoms'))
outFile.write('12 atom types \n')
#outFile.write('%i %s \n' %(nbonds, 'bonds'))
#outFile.write('1 bond types \n')
#outFile.write('%i %s \n' %(nangles, 'angles'))
#outFile.write('1 angle types \n')
outFile.write('\n')
outFile.write('%f %f %s %s \n' %(xmin, xmax, 'xlo', 'xhi'))
outFile.write('%f %f %s %s \n' %(ymin, ymax, 'ylo', 'yhi'))
outFile.write('%f %f %s %s \n' %(zmin, zmax, 'zlo', 'zhi'))
#outFile.write('\n')
#outFile.write('%f %f %f %s %s %s \n' %(-8.1554775, 0.0, 0.0, 'xy', 'xz', 'yz'))
outFile.write('\n')
outFile.write('Masses\n')
outFile.write('\n')
outFile.write('%i %f \n' %(1, 28.065))
outFile.write('%i %f \n' %(2, 26.98))
outFile.write('%i %f \n' %(3, 22.98))
outFile.write('%i %f \n' %(4, 16.0))
outFile.write('%i %f \n' %(5, 16.0))
outFile.write('%i %f \n' %(6, 1.0))
outFile.write('%i %f \n' %(7, 28.065)) #Si
outFile.write('%i %f \n' %(8, 40.078)) #Ca
outFile.write('%i %f \n' %(9, 40.078)) #Ca
outFile.write('%i %f \n' %(10, 16.0)) #Os
outFile.write('%i %f \n' %(11, 16.0)) #Ow
outFile.write('%i %f \n' %(12, 1.0)) #Hw
outFile.write('\n')
outFile.write('Atoms\n')
outFile.write('\n')


for j in range(size1): #writing atoms without water
    if datazero[j,1]==1:
        outFile.write('%i %i %i %f %f %f %f \n' %(j+1, 1, 1, 0, datazero[j,2], datazero[j,3], datazero[j,4]))
    elif datazero[j,1]==2:
        outFile.write('%i %i %i %f %f %f %f \n' %(j+1, 1, 2, 0, datazero[j,2], datazero[j,3], datazero[j,4]))
    elif datazero[j,1]==3:
        outFile.write('%i %i %i %f %f %f %f \n' %(j+1, 1, 3, 0, datazero[j,2], datazero[j,3], datazero[j,4]))
    elif datazero[j,1]==4:
        outFile.write('%i %i %i %f %f %f %f \n' %(j+1, 1, 4, 0, datazero[j,2], datazero[j,3], datazero[j,4]))
    elif datazero[j,1]==5:
        outFile.write('%i %i %i %f %f %f %f \n' %(j+1, 1, 5, 0, datazero[j,2], datazero[j,3], datazero[j,4]))
    elif datazero[j,1]==6:
        outFile.write('%i %i %i %f %f %f %f \n' %(j+1, 1, 6, 0, datazero[j,2], datazero[j,3], datazero[j,4]))
    elif datazero[j,1]==7:
        outFile.write('%i %i %i %f %f %f %f \n' %(j+1, 1, 7, 0, datazero[j,2], datazero[j,3], datazero[j,4]))
    elif datazero[j,1]==8:
        outFile.write('%i %i %i %f %f %f %f \n' %(j+1, 1, 8, 0, datazero[j,2], datazero[j,3], datazero[j,4]))             
    elif datazero[j,1]==9:
        outFile.write('%i %i %i %f %f %f %f \n' %(j+1, 1, 9, 0, datazero[j,2], datazero[j,3], datazero[j,4]))
    elif datazero[j,1]==10:
        outFile.write('%i %i %i %f %f %f %f \n' %(j+1, 1, 10, 0, datazero[j,2], datazero[j,3], datazero[j,4]))
    elif datazero[j,1]==11:
        outFile.write('%i %i %i %f %f %f %f \n' %(j+1, 1, 11, 0, datazero[j,2], datazero[j,3], datazero[j,4]))
    elif datazero[j,1]==12:
        outFile.write('%i %i %i %f %f %f %f \n' %(j+1, 1, 12, 0, datazero[j,2], datazero[j,3], datazero[j,4])) 

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

