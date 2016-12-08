
#this script will add -OH group to the undercoordinated Si and Al atoms on the surface of indentation
#motivation for this script is to conduct MD simulation of nanoindentation of geopolymer and study the effect of -OH on surface

def Chargebalance(x):
    """x is the array of x y z coordinates with the atom type column."""
    
    #calculate the charge imbalance (Si = +2.4, Al = +1.8, Na = +0.6, O = -1.2)
    print "Calculating charge balance..."

    charges = {'1.0':2.4, '2.0':1.8, '3.0':0.6, '4.0':-1.2, '5.0':-0.856, '6.0':-0.82, '7.0':0.256, '8.0':0.41, '9.0':0.6} #change the dic according to atom type
    Totalcharge = 0.0
    for i in range(len(x)):
        Totalcharge += charges[str(x[i,0])]
    return Totalcharge

def delAtoms(x,t,n):
    """x = array of xyz data file with the atom types,
    t = type of atom to be deleted,
    n = no. of atoms to be deleted."""
    
   
    idxx = []
    
    for i in range(len(x)):
        if x[i,0] == t:
            idxx.append(i)
        else:
            idxx = idxx
    
    print "no of Na atoms ="+ str(len(idxx))
    print "deleted Na = "+str(n)
    return np.delete(x,idxx[:n],0)


import numpy as np
import glob
    
#ratio = raw_input("What is the Si/Al ratio?")
##radius = raw_input("What is the radius of void?")
#cen = tuple(int(x.strip()) for x in raw_input("Enter the center coordinates of void: ").split(','))
print "Enter the center coordinates:"
cx = float(raw_input("Enter x: "))
cy = float(raw_input("Enter y: "))
cz = float(raw_input("Enter z: "))
cen = [cx,cy,cz]
zmax = 89.0  # point above surface from which the distance will be measured to add atoms
fname = "Nlist-GBP2w0.44b.dat"
natoms_GBP = 40600

num_lines = 0
num_words = 0
num_chars = 0
idx1_3 = [] #list of atoms with less than 4 neighbors
idx2_3 = []
idx1_2 = []
idx2_2= []
idx1_1 = []
idx2_1 = []
with open(fname, 'r') as f:
    for line in f:
        num_lines += 1        
        words = line.split() # list of each lines
         # count the line numbers
        if words[1]==str(2): # condition for chosing only atom type 2      
          
          num_words = len(words)
          if num_words ==5: # if the number of neighbor is 3
            idx1_3.append(num_lines)
          elif num_words ==4: # if only 2 neighbors are present
              idx1_2.append(num_lines)
        elif words[1]==str(1): # condition for chosing only atom type 1      
          
          num_words = len(words)
          if num_words ==5: # if the number of neighbor is 3
            idx2_3.append(num_lines) 
          elif num_words ==4: #if only 2 neighbors are present
              idx2_2.append(num_lines) 

       
print "Al(3): = ",len(idx1_3) # prints the list of atom IDs with less than 4 neighbors for a given atomtype
#print idx2_3        
print "Si(3): = ",len(idx2_3)
#print idx1_2      # prints the list of atom IDs with less than 3 neighbors for a given atomtype  
print "Al(2): = ",len(idx1_2) 
#print idx2_2        
print "Si(2): = ",len(idx2_2)

idx3 = idx1_3+idx2_3
idx2 = idx1_2+idx2_2

print "Total 3-coordinated atoms = ",len(idx3)
print "Total 2-coordinated atoms = ",len(idx2)
idx3_1 = [] #list of 3-coordinated atoms near surface
idx2_1 = [] #list of 2-coordinated atoms near surface


filename = "data.GP2.0w0.44br1-Eq" #serially written lammps format data file after making voids
flist = glob.glob(filename)
 
#for f in flist: 
#    load = np.loadtxt(f, usecols=(0,1,2,3,4,5,6))
    
#data=np.array(load)                        

for f in flist: 
    load1 = np.genfromtxt(f, dtype=float, skip_header=35, skip_footer=natoms_GBP+1, usecols=(4,5,6))
radius = float(raw_input("What is the range from surface for adding OH? ")) #chose atoms 1-2 A from surface i.e. for 89 zmax and 85 max z coordinate chose 6 A; 89-85 = 4, so for 6 A the penetration is only 6-4=2 A  
data_xyz=np.array(load1) 
#center = [50.0,45.0,45.0]
for i in range(len(idx3)):
    if np.linalg.norm(zmax-data_xyz[idx3[i]-1,2]) < radius: #restrict the positions of adding -OH only near the surface
        idx3_1.append(idx3[i])
    else:
        idx3_1 = idx3_1
print "3-coordinated atoms near surface = ",len(idx3_1)
#print idx3
#print idx3_1
'''for i in range(len(idx2)):
    if np.linalg.norm(data_xyz[idx2[i]]-zmax) < radius: #restrict the positions of adding -OH only near the surface
        idx2_1.append(idx2[i])
    else:
        idx2_1 = idx2_1
print "2-coordinated atoms near surface = ",len(idx2_1)'''
print "************only adding near the surface**************"
#---find the coordinates for 3-folded atoms---------------

s = (len(idx3_1),3)
coord_array = np.zeros(s,dtype=float) #array for storing the xyz coordinates
center = np.zeros(s,dtype=float) #center of the void
d = (len(idx3_1),1)
dist_array = np.zeros(d,dtype=float)

###change the coord array to contain only z vector

for i in xrange(len(idx3_1)):
        coord_array[i,:] = data_xyz[idx3_1[i]-1,:] #####deduct 1 as numpy starts at 0

for i in xrange(len(idx3_1)):
    center[i,:] = cen
#print center
vector = center - coord_array

for i in xrange(len(idx3_1)):
    
    dist_array[i] = np.linalg.norm(vector[i])

unit_vector = vector / dist_array 

new_coord_O = coord_array + unit_vector * 1.8 # adding O(5) at a dist of 1.8 A
    
new_coord_H = coord_array + unit_vector * 2.8 # adding H(7) at a dist of 2.8 A        
    
new_coord_Na = coord_array + unit_vector * 4.0 #adding Na (3) atom for charge balance



#------*****find the coordinates for 2-folded atoms*******---------------------------------
#*******************************************************************************

nOz = len(idx2_1) #no. of Oz to add 
s2 = (nOz,3)
coord_array2 = np.zeros(s2,dtype=float) #array for storing the xyz coordinates
center2 = np.zeros(s2,dtype=float) #center of the void
d = (nOz,1)
dist_array2 = np.zeros(d,dtype=float)

for i in xrange(nOz):
     #if 0.0 < data_xyz[i,0] < 88.0 or 0.0 < data_xyz[i,0] < 88.0: 
        coord_array2[i,:] = data_xyz[idx2_1[i]-1,:] #####deduct 1 as numpy starts at 0
     #else:
        #coord_array2 = coord_array2
#print len(coord_array2)
for i in xrange(nOz):
    center2[i,:] = cen  #a different center so that atoms dont overlap (i.e. Oz with O(oh))

vector2 = center2 - coord_array2

for i in xrange(nOz):
    
    dist_array2[i] = np.linalg.norm(vector2[i])


unit_vector2 = vector2 / dist_array2

new_coord_Oz = coord_array2 + unit_vector2 * 1.8 # adding O(4) at a dist of 1.8 A
new_coord_Hz = coord_array2 + unit_vector2 * 2.8 # adding H(7) at a dist of 2.8 A 
new_coord_Naz = coord_array2 + unit_vector2 * 4.0 # adding Na(9) at a dist of 3.8 A
new_coord_Nax = coord_array2 + unit_vector2 * 10.0 # extra Na (9) for charge balance

#----------------- combine the data into a single array with atomtypes------------------
#********************************************************************************

# add the atom types before the xyz data

new_coord_O = np.insert(new_coord_O,0,5,axis=1)
new_coord_H = np.insert(new_coord_H,0,7,axis=1)
new_coord_Na = np.insert(new_coord_Na,0,9,axis=1) # changed the Na to type 9 to better study its transport apart from the bulk Na
new_coord_Oz = np.insert(new_coord_Oz,0,5,axis=1) # change to 4 or 5 depending on the need
new_coord_Hz = np.insert(new_coord_Hz,0,7,axis=1)
new_coord_Naz = np.insert(new_coord_Naz,0,9,axis=1)
new_coord_Nax = np.insert(new_coord_Nax,0,9,axis=1)



#load the original data 

for f in flist: 
    load = np.genfromtxt(f, skip_header=35, skip_footer=natoms_GBP+1, usecols=(2,4,5,6))
    load_all = np.genfromtxt(f, skip_header=35, skip_footer=natoms_GBP+1, usecols=(0,1,2,4,5,6)) 
     
data_original = np.array(load) 
data_original_all = np.array(load_all) #need this for retaining all original molid, atom serial etc.


nH2O = int(raw_input("How many water molecules? "))
size1 = len(data_original)-3*nH2O #total atoms-water
size2 = len(data_original) #total atoms
natoms = size2
nbonds = nH2O*2
nangles = nH2O

#combine the final data
#do adjustment here to ensure charge balance
dNa = int(raw_input("How many Na you want to relocate? "))
data_original1 = delAtoms(data_original,3,dNa) # after deleting bulk Na
#print len(data_original), len(data_original1)

final_data = np.vstack((data_original1,new_coord_O,new_coord_H,new_coord_Oz,new_coord_Hz,new_coord_Na))#, new_coord_Na, new_coord_Naz[24:])) 
new_data = np.vstack((new_coord_O,new_coord_H,new_coord_Oz,new_coord_Hz,new_coord_Na))

#calculate the total charge balance

charge = Chargebalance(final_data) + dNa*0.6 #(add dNa using packmol)
sNa = -1*charge/0.6 ##see comment below:

#if charge balance is positive it will deduct that many Na
# if negative it will become positive and ask to add that amount of Na finally

print "Total charge imbalance: ", charge


#----------------- write the final lammps data file--------------------------------------

natoms=len(final_data)
xmin = np.min(data_xyz[:,0])
xmax = np.max(data_xyz[:,0])
ymin = np.min(data_xyz[:,1])
ymax = np.max(data_xyz[:,1])
zmin = np.min(data_xyz[:,2])
zmax = np.max(data_xyz[:,2])+80.0

outFile = open('data.GBP2w0.44-OH', 'w')
outFile.write('LAMMPS data file written by Rafat Sadat using Python\n')
outFile.write('\n')
outFile.write('%i %s \n' %(natoms, 'atoms'))
outFile.write('9 atom types \n')
outFile.write('%i %s \n' %(nbonds, 'bonds'))
outFile.write('1 bond types \n')
outFile.write('%i %s \n' %(nangles, 'angles'))
outFile.write('1 angle types \n')
outFile.write('\n')
outFile.write('%f %f %s %s \n' %(xmin, xmax, 'xlo', 'xhi'))
outFile.write('%f %f %s %s \n' %(ymin, ymax, 'ylo', 'yhi'))
outFile.write('%f %f %s %s \n' %(zmin, zmax, 'zlo', 'zhi'))
outFile.write('\n')
outFile.write('Atoms\n')
outFile.write('\n')

  
for j in range(len(data_original_all)):
        outFile.write('%i %i %i %i %f %f %f \n' %(data_original_all[j,0], data_original_all[j,1], data_original_all[j,2], 0, data_original_all[j,3], data_original_all[j,4], data_original_all[j,5])) #atom serial is being renumbered!      

for j in range(len(data_original_all),len(final_data)): #writing the new data
        outFile.write('%i %i %i %i %f %f %f \n' %(j+1, 1, final_data[j,0], 0, final_data[j,1], final_data[j,2], final_data[j,3])) #atom serial is being renumbered!            
                        
outFile.write('\n')
outFile.write('Bonds\n')
outFile.write('\n')

count = 0
for j in range(size1+1,size2, 3):
    count +=1
    outFile.write('%i %s %i %i \n' %(count,'1',j,j+2))
    count +=1
    outFile.write('%i %s %i %i  \n' %(count,'1',j+1,j+2))

outFile.write('\n')
outFile.write('Angles\n')
outFile.write('\n')   

count = 0
for j in range(size1+1,size2, 3):
    count +=1
    outFile.write('%i %s %i %i %i \n' %(count,'1', j, j+2, j+1))
   
outFile.close()
print "All data written!"                                                          
outFile.close()

# write in xyz format
'''
outFile = open('GP-1.4-0-25A.xyz', 'w')
outFile.write('%i \n' %(len(final_data)))
outFile.write('\n')

for j in range(len(final_data)):

    if final_data[j,0]==1:
        outFile.write('%s %f %f %f \n' %('Si', final_data[j,1], final_data[j,2], final_data[j,3]))
    elif final_data[j,0]==2:
        outFile.write('%s %f %f %f \n' %('Al', final_data[j,1], final_data[j,2], final_data[j,3]))
    elif final_data[j,0]==3:
        outFile.write('%s %f %f %f \n' %('Na', final_data[j,1], final_data[j,2], final_data[j,3]))
    elif final_data[j,0]==4:
        outFile.write('%s %f %f %f \n' %('O', final_data[j,1], final_data[j,2], final_data[j,3]))
    elif final_data[j,0]==5:
        outFile.write('%s %f %f %f \n' %('Oh', final_data[j,1], final_data[j,2], final_data[j,3]))
    elif final_data[j,0]==6:
        outFile.write('%s %f %f %f \n' %('Ow', final_data[j,1], final_data[j,2], final_data[j,3]))
    elif final_data[j,0]==7:
        outFile.write('%s %f %f %f \n' %('H', final_data[j,1], final_data[j,2], final_data[j,3]))
    elif final_data[j,0]==8:
        outFile.write('%s %f %f %f \n' %('Hw', final_data[j,1], final_data[j,2], final_data[j,3]))
    elif final_data[j,0]==9:
        outFile.write('%s %f %f %f \n' %('Nap', final_data[j,1], final_data[j,2], final_data[j,3]))
outFile.close()
'''
#calculating the final Si/Al/Na ratios

nSi = 0
nAl = 0
nNa = 0
nO = 0
nOh = 0
nH = 0
nNaP = dNa #these are relocated Na from bulk
for j in range(len(final_data)):

  if final_data[j,0] == 1:
            nSi +=1
  elif final_data[j,0] ==2:
            nAl +=1
  elif final_data[j,0] ==3:
            nNa +=1
  elif final_data[j,0] ==4:
            nO +=1
  if final_data[j,0] == 5:
            nOh +=1
  elif final_data[j,0] ==7:
            nH +=1
  elif final_data[j,0] ==9:
            nNaP +=1

SiAl = nSi/float(nAl)
NaAl = (nNa+nNaP) / float(nAl)
print "No. of Si = ", nSi
print "No. of Al = ", nAl
print "No. of Na = ", nNa
print "No. of O = ", nO
print "No. of NaP = ", nNaP
print "No. of hydroxide O = ", nOh
print "No. of hydroxide H = ", nH
print "Si/Al ratio=", SiAl
print "Na/Al ratio=",NaAl
print "All done!"
print "Remarks: Don't forget to add "+ str(dNa+sNa) + " no.s of Na with H2O using packmol!"
#coordination status

