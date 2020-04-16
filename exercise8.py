###############################################################################
# Exercise 8
###############################################################################

import matplotlib . pyplot as plt
import numpy as np
#import time
#import pandas as pd
#from mpl_toolkits.mplot3d import Axes3D  
#from matplotlib import cm

no=22;
N=500;
incR=0.01;
r=0;
x0=0.1;
x1=0.1;

###############################################################################
x = np.zeros((N), dtype=np.float32)
x[0]=x0;
x[1]=x1;


k=0;
#ratio=[0]*32;
#popul=[0]*32; 
ratio = []
popul=[]; 

# slow since we're assigning memory to list as we go forward in loop
while(r<4.0):
  # calculating convergence of x for r  
  for i in range (0,N-1):
      x[i+1]= r*x[i]*(1-x[i]);
    
  # storing the final 30 values of x for r
  for j in range(N-30,N-1):
      
    #ratio[k]=r;     #popul[k]=x[j];  
    ratio.append([r])
    popul.append([x[j]])
    k=k+1;
  
  # increment in r and restart x
  r = r + incR;
  x[0] = 0.1;

np_ratio = np.array(ratio)
np_pop = np.array(popul)
###############################################################################
  

fig = plt.figure(num=no,figsize=(5, 5)) ;
plt.clf();
#plt.subplot(311);
plt.plot(ratio,popul,'.')
plt.xlabel('r')
plt.ylabel('x$_n$')
#plt.title(r'Conergence to $\pi $')
plt.title(r'Logistic Map x$_n$$_+$$_1$=rx$_n$(1-x$_n$)')
plt.grid()

###############################################################################
###############################################################################
###############################################################################
