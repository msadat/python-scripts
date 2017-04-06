#this python script first creates a neighborlist of atoms using Periodic CKDTree search
#next one can get the Qn distribution of a Si-O-Al network or Si-O-Si network using this script
#Author: Mohammad Rafat Sadat, Dept. of CEEM University of Arizona
#Date: 3.31.17

from periodic_kdtree import PeriodicCKDTree
import numpy as np

import glob
import time

outputfile = 'GP-NaAl1.0'
global natoms_GBP 

natoms_GBP = 8600
    
start_time = time.time()    

global filename

filename = 'data.GP2-w2.17-t0.25-Eq' #original lammps datafile

flist = glob.glob(filename)

for f in flist:
    load = np.genfromtxt(f, dtype=float, skip_header=19, skip_footer=natoms_GBP+1, usecols=(4,5,6))
    #load = np.genfromtxt(f, dtype=float, skip_header=18, usecols=(4,5,6))
x=np.array(load)

for f in flist: #read the finally written lammps data file to make the neighborlist
    load = np.genfromtxt(f, dtype=int, skip_header=19, skip_footer=natoms_GBP+1, usecols=(0,2)) #only need the serial and atom type col
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


#find the coordination of Si and Al atoms

def coord(fname, datafile):
    """Nlist file name, lammps data file name"""
    num_lines = 0
    num_words = 0
    
    Al6 = []
    Al5 = []
    Al4 = []
    Al3 = []
    Al2 = []
    Si6 = []
    Si5 = []
    Si4 = []
    Si3 = []
    Si2 = []
    Si1 = []
    BOs = []
    NBOs = []

    with open(fname, 'r') as f:
      for line in f:
        num_lines += 1        
        words = line.split() # list of each lines
         # count the line numbers
        if words[1]==str(2): # condition for chosing only atom type 2      

          num_words = len(words) #-count #getting rid of the like neighbors

          if num_words ==5: # if the number of neighbor is 3
            Al3.append(num_lines)
          elif num_words ==4: # if only 2 neighbors are present
              Al2.append(num_lines)
          elif num_words ==6:
              Al4.append(num_lines)
          elif num_words ==7:
              Al5.append(num_lines)  
          elif num_words ==8:
              Al6.append(num_lines)
              
        elif words[1]==str(1): # condition for chosing only atom type 1   

          num_words = len(words) #-count

          if num_words ==5: # if the number of neighbor is 3
            Si3.append(num_lines) 
          elif num_words ==4: #if only 2 neighbors are present
              Si2.append(num_lines) 
          elif num_words ==6:
              Si4.append(num_lines)
          elif num_words ==7:
              Si5.append(num_lines)  
          elif num_words ==8:
              Si6.append(num_lines)
          elif num_words ==3:
              Si1.append(num_lines)   
              
        elif words[1]==str(4): # condition for chosing only atom type 4      
          
          num_words = len(words)
          if num_words == 4: # if the number of neighbor is 2
            BOs.append(num_lines) 
          elif num_words == 3: #if only 1 neighbors are present
              NBOs.append(num_lines)
  
    flist = glob.glob(datafile)

    for f in flist: 
      load = np.genfromtxt(f, skip_header=19, skip_footer = natoms_GBP+1, usecols=(2,4,5,6))
    
    data=np.array(load)                        


    #calculating the Si/Al/Na ratios

    nSi = 0
    nAl = 0
    nNa = 0
    nO = 0


    for j in range(len(data)):

     if data[j,0] == 1:
            nSi +=1
     elif data[j,0] ==2:
            nAl +=1
     elif data[j,0] ==3:
            nNa +=1
     elif data[j,0] ==4:
            nO +=1

    SiAl = nSi/float(nAl)
    NaAl = (nNa) / float(nAl)
    print "***********elements************" 
    print "No. of Si = ", nSi
    print "No. of Al = ", nAl
    print "No. of Na = ", nNa
    print "No. of O = ", nO
    print "Si/Al ratio=", SiAl
    print "Na/Al ratio=",NaAl

    #print "***********coordination status of Al***********" 
    #print "Al(6): = ", len(Al6)
    #print "%Al(6): = ", len(Al6)*100.0/nAl
    #print "Al(5): = ", len(Al5)
    #print "%Al(5): = ", len(Al5)*100.0/nAl
    #print "Al(4): = ",len(Al4)
    print "%Al(4): = ", len(Al4)*100.0/nAl       
    #print "Al(3): = ",len(Al3)
    #print "%Al(3): = ", len(Al3)*100.0/nAl
    #print "Al(2): = ",len(Al2)
    #print "%Al(2): = ", len(Al2)*100.0/nAl
    #print "Total Al: = ", len(Al6)+len(Al5)+len(Al4)+len(Al3)+len(Al2)
    #print "***************coordination status of Si**************************"
    #print "Si(6): = ",len(Si6)
    #print "%Si(6): = ", len(Si6)*100.0/nSi  
    #print "Si(5): = ",len(Si5)
    #print "%Si(5): = ", len(Si5)*100.0/nSi 
    #print "Si(4): = ",len(Si4)
    print "%Si(4): = ", len(Si4)*100.0/nSi        
    #print "Si(3): = ",len(Si3)
    #print "%Si(3): = ", len(Si3)*100.0/nSi        
    #print "Si(2): = ",len(Si2)
    #print "%Si(2): = ", len(Si2)*100.0/nSi
    #print "Si(1): = ",len(Si1)
    #print "%Si(1): = ", len(Si1)*100.0/nSi
    #print "Total Si: = ", len(Si6)+len(Si5)+len(Si4)+len(Si3)+len(Si2)+len(Si1)
    #print "NBOs: = ",len(NBOs)
    print "%NBOs: = ", len(NBOs)*100.0/nO

    
    return len(Al4)*100.0/nAl, len(Si4)*100.0/nSi

Q6, Q5, Q4, Q3, Q2, Q1, Q0 = Qn_distrib('Nlist-ID_Si_Al-'+outputfile, 'Nlist-ID_O-'+outputfile, 'Nlist-type_O-'+outputfile)

Al4, Si4 = coord('Nlist-types-GP-NaAl1.0', filename)

def outputwrite(outfilename):
        
        outFile = open('GB-coord-Qn-'+outfilename+'.dat', 'a')
        outFile.write('%s %s %s %s %s %s \n' %(nQ5, nQ4, nQ3, nQ2, Si4, Al4))
        outFile.write('%s %f \n' %(outname, len(Al4)*100.0/nAl))
        outFile.close()
        return 0
