import numpy as np
import random
import operator as op

def ncr(n, r):
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, xrange(n, n-r, -1))
    denom = reduce(op.mul, xrange(1, r+1))
    return numer//denom

types = 10
#ntypes = ncr(types,2) + types
#pairs = np.zeros((ntypes,2))

#x = np.linspace(1,types,types)
#y = np.linspace(1,types,types)

pairsxy = np.zeros((types**2,2))

for i in range(1,types+1):
     for j in range(1,types+1):
            k = types*(i-1) + j-1
            pairsxy [k,0] = i 
            pairsxy [k,1] = j 

remove_list = []            
for x in range(len(pairsxy)):
  if pairsxy[x,0] > pairsxy[x,1]:
       remove_list.append(x)
              
pairsxy = np.delete(pairsxy,remove_list,0)

print len(pairsxy)

      

outputfile = "in.peri-lammps"
horizon = 30.01 #nm

c = np.random.normal(56,0.9) #spring constant from distribution

s00 = np.random.normal(0.005,0.0001) #critical stretch


outFile = open(outputfile, 'w')

for j in range(len(pairsxy)):
        c = np.random.normal(56,0.9)
        s00 = np.random.normal(0.005,0.0001)
        outFile.write('%s %i %i %f %f %f %f \n' %("pair_coeff", int(pairsxy[j,0]), int(pairsxy[j,1]), c, horizon, s00, 0.25 )) 

outFile.close()
