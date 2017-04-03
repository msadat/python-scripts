#this python script first creates a neighborlist of atoms using Periodic CKDTree search
#next one can get the Qn distribution of a Si-O-Al network or Si-O-Si network using this script
#Author: Mohammad Rafat Sadat, Dept. of CEEM University of Arizona
#Date: 3.31.17

from periodic_kdtree import PeriodicCKDTree
import numpy as np

import glob
import time

outputfile = 'silica'
natoms_GBP = 17496
    
start_time = time.time()    
filename = 'data.Silica-Eq'

flist = glob.glob(filename)

for f in flist:
    load = np.genfromtxt(f, dtype=float, skip_header=20, skip_footer=natoms_GBP+1, usecols=(4,5,6))
    #load = np.genfromtxt(f, dtype=float, skip_header=18, usecols=(4,5,6))
x=np.array(load)

for f in flist: #read the finally written lammps data file to make the neighborlist
    load = np.genfromtxt(f, dtype=int, skip_header=20, skip_footer=natoms_GBP+1, usecols=(0,2)) #only need the serial and atom type col
    #load = np.genfromtxt(f, skip_header=18, usecols=(0,2))
y=np.array(load)


#--------------list the IDs of Si and Al-----------------

Al_list = []
Si_list = []

for t in range(len(y)):
    if y[t,1]==1:
        Si_list.append(y[t,0])
    elif y[t,1] ==2:
        Al_list.append(y[t,0])
    else:
        Si_list = Si_list
        Al_list = Al_list

SiAl_list = Si_list + Al_list

#------------------------------------------------------------

                       
xmin = np.min(x[:,0])
xmax = np.max(x[:,0])
dx = xmax-xmin
ymin = np.min(x[:,1])
ymax = np.max(x[:,1])
dy = ymax-ymin
zmin = np.min(x[:,2])
zmax = np.max(x[:,2])
dz = zmax-zmin
s = (len(x),10)
Nlist = np.zeros(s,dtype=np.int)
  # Boundaries (0 or negative means open boundaries in that dimension)
#changing bounds manually
bounds = np.array([dx, dy, dz])   # xy periodic, open along z


  # Build kd-tree
T = PeriodicCKDTree(bounds, x)

  # Find 4 closest neighbors to a random point
  # (d[j], i[j]) = distance and index of jth closest point
  # Find neighbors within a fixed distance of a point
print "Building Neighborlist..."  

neighbors = []
for i in xrange(len(x)):
  localneigh = T.query_ball_point(x[i],r=2.1) #r = cutoff (Angstrom) for making Nlist
  #localneigh.insert(0,i)
  localneigh.remove(i)
  localneigh.insert(0,i)  
  neighbors.append(localneigh)
  
#print neighbors

print "Neighborlist built! Writing data to file...."
print "***********writing with atom types*****************"
#print neighbors
outFile = open('Nlist-types'+'-'+outputfile, 'w')

for i in xrange(s[0]):
        #Slice the atomtypes using the neighbor indices, have to subtract 1 
        #from index because you added it in your neighborlist build
	neightypes = y[neighbors[i],1]
        #print neightypes
	#Now loop of the sliced neighbor types
        if y[i,1] == 1 or y[i,1] == 2 or y[i,1] == 4:
            outFile.write("%i " % y[i,0])
            for items in neightypes:
                #outFile.write("%i " % y[i,0])
		outFile.write("%i " % items)
	    outFile.write("\n")

outFile.close()

print "***********writing with atom type O only*****************"
#print neighbors
outFile = open('Nlist-type_O'+'-'+outputfile, 'w')

for i in xrange(s[0]):
        #Slice the atomtypes using the neighbor indices, have to subtract 1 
        #from index because you added it in your neighborlist build
	neightypes = y[neighbors[i],1]
        #print neightypes
	#Now loop of the sliced neighbor types
        if y[i,1] == 4:
            outFile.write("%i " % y[i,0])
            for items in neightypes:
                #outFile.write("%i " % y[i,0])
		outFile.write("%i " % items)
	    outFile.write("\n")

outFile.close()
print "***********writing with atom ID*****************"
#print neighbors
outFile = open('Nlist-ID'+'-'+outputfile, 'w')

for i in xrange(s[0]):
        #Slice the atomtypes using the neighbor indices, have to subtract 1 
        #from index because you added it in your neighborlist build
	neightypes = y[neighbors[i],0]
        #print neightypes
	#Now loop of the sliced neighbor types
        if y[i,1] == 1 or y[i,1] == 2 or y[i,1] == 4:
            outFile.write("%i " % y[i,0])
            for items in neightypes:
                #outFile.write("%i " % y[i,0])
		outFile.write("%i " % items)
	    outFile.write("\n")

outFile.close()

print "***********writing with atom ID of O only*****************"
#print neighbors
outFile = open('Nlist-ID_O'+'-'+outputfile, 'w')

for i in xrange(s[0]):
        #Slice the atomtypes using the neighbor indices, have to subtract 1 
        #from index because you added it in your neighborlist build
	neightypes = y[neighbors[i],0]
        #print neightypes
	#Now loop of the sliced neighbor types
        if y[i,1] == 4:
            outFile.write("%i " % y[i,0])
            for items in neightypes:
                #outFile.write("%i " % y[i,0])
		outFile.write("%i " % items)
	    outFile.write("\n")

outFile.close()

print "***********writing with atom ID of Si or Al only*****************"
#print neighbors
outFile = open('Nlist-ID_Si_Al'+'-'+outputfile, 'w')

for i in xrange(s[0]):
        #Slice the atomtypes using the neighbor indices, have to subtract 1 
        #from index because you added it in your neighborlist build
	neightypes = y[neighbors[i],0]
        #print neightypes
	#Now loop of the sliced neighbor types
        if y[i,1] == 1 or y[i,1] == 2 :
            outFile.write("%i " % y[i,0])
            for items in neightypes:
                #outFile.write("%i " % y[i,0])
		outFile.write("%i " % items)
	    outFile.write("\n")

outFile.close()

##*********************************FINDING THE Qn DISTRIBUTION*********************************

#1. compare between the two array and see for each Si or Al, if a neighboring O has a different Si or Al, if TRUE then count +=1
# if count==4 then Q4 +=1, if count ==3, then Q3 += 1 and so on.... (if type ==1, Q4Si, elif type ==2, Q4Al...)
#2. loop no. 1 for all the neighboring Oxygens of that particular Si or Al. 
#3. Repeat for all Si and Al.

def Qn_distrib(SiAlID, OID, Otype):
    
    #nSiAl = nSi + nAl
    #Q1_Al = []
    #Q1_Si = []
    #Q2_Al = []
    #Q2_Si = []
    #Q3_Al = []
    #Q3_Si = []
    #Q4_Al = []
    #Q4_Si = []
    Q4 = 0
    Q3 = 0
    Q2 = 0
    Q1 = 0
    Q0 = 0
    Q5 = 0
    Q6 = 0
    
    SiAl_full = []
    O_full = []
    Otype_full = []
    Obridging = [] #list of bridging O with two Si or Al neighbors (i.e. Si-O-Al)
    count = 0
    #bridgeO = 0
    with open(SiAlID,  'r') as f:
     for line in f:
                
        words_SiAl = line.split() # list of each lines
        
        for x in range(len(words_SiAl)):
            words_SiAl[x] = int(words_SiAl[x])
        SiAl_full.append(words_SiAl)

         
    with open(Otype,  'r') as f:    
      for line in f:
                
        words_Otype = line.split() # list of each lines
      
        for x in range(len(words_Otype)):
            words_Otype[x] = int(words_Otype[x]) 
        Otype_full.append(words_Otype)
    
    ###Make a list of bridging oxygens**********************
         
    for i in range(len(Otype_full)):
        bridgeO=0
        for j in range(2,len(Otype_full[i])):
            
            if Otype_full[i][j] == 1 or Otype_full[i][j] == 2: 
              bridgeO +=1
              if bridgeO >=2:  
                Obridging.append(Otype_full[i][0])
              else:
                Obridging = Obridging  
    
    #print Otype_full           
    #print Obridging
    
    with open(OID,  'r') as f:
      for line in f:
               
        words_O = line.split() # list of each lines
        
        for x in range(len(words_O)):
            words_O[x] = int(words_O[x])
        O_full.append(words_O)
    
    #now we have Nested lists of all SiAl and O neighbors
    #each element on this list is the list of neighbors that the atom has
    #for example, SiAl_full is a list of neighborlists where the first two IDs are of the Si or Al atom itself and the rests are O neighbors
    #O_full is the list of O neighbors where the first two is the O itself followed by it's neighbors
    #print O_full
    for i in range(len(SiAl_full)):
        count = 0 
        for j in range(2,len(SiAl_full[i])):
          
          for k in range(len(Obridging)):  
             #  
             if SiAl_full[i][j] == Obridging[k]:                 
               count +=1
             else:
               count = count  
        if count == 4:
             Q4 +=1
        elif count == 3:
             Q3 +=1
        elif count == 2:
             Q2 +=1
        elif count == 1:
             Q1 +=1
        elif count  == 0:
             Q0 +=1
        elif count == 5:
             Q5 +=1 
        elif count == 6:
             Q6 +=1
          
    return Q6, Q5, Q4, Q3, Q2, Q1, Q0


#SiAllist, Olist = Qn_distrib('Nlist-ID_Si_Al-GP2C2', 'Nlist-ID_O-GP2C2')
Q6, Q5, Q4, Q3, Q2, Q1, Q0 = Qn_distrib('Nlist-ID_Si_Al-'+outputfile, 'Nlist-ID_O-'+outputfile, 'Nlist-type_O-'+outputfile)
#print SiAllist
#print Olist
print "Q6 Q5 Q4 Q3 Q2 Q1 Q0: "
print Q6, Q5, Q4, Q3, Q2, Q1, Q0

print "All done!"
print("--- %s seconds ---" % (time.time() - start_time)) 
