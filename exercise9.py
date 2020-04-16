###############################################################################
# Exercise 9 Mandelbrot set
###############################################################################

import matplotlib . pyplot as plt
import numpy as np
#import time
#import pandas as pd
#from mpl_toolkits.mplot3d import Axes3D  
#from matplotlib import cm

###############################################################################
###############################################################################
# edit to include these as flags or inputs
no=5;
startr=-2; # startin real plane
startc=-1.5; # start in complex plane
finishr=1; # finish in real plane
finishc=1.5; # finish in complex plane
inc=0.01;   # increment in both planes
# to zoom in change the numbers above

###############################################################################
###############################################################################
Nx = int((finishr-startr)/inc)
Ny = int((finishc-startc)/inc)
X = np.zeros((Nx,Ny), dtype=np.float32)
Y = np.zeros((Nx,Ny), dtype=np.float32)
Z = np.zeros((Nx,Ny), dtype=np.float32)

# slow loop -why? time it and use cython to speed up!
for i in range(0,Nx):
    for j in range(0,Ny):
        x=startr+i*inc;
        y=startc+j*inc;
        X[i,j]=x
        Y[i,j]=y
        flag=0;
        u=0;
        v=0;
        r=0.0;
 
        while(flag<=2000 and r<400):
            u1=(u*u)-(v*v)+x;
            v1=(2*u*v)+y;
            u=u1;
            v=v1;
            r=np.sqrt( (u*u)+(v*v) );
            flag=flag+1;
 
        Z[i,j]=flag;

###############################################################################
###############################################################################
# pcolor plot
fig = plt.figure(num=no,figsize=(5, 5)) ;
plt.clf();
plt.pcolor(X,Y,Z,vmin=0,vmax=50);
plt.xlabel('Real Place')
plt.ylabel('Imag Plane')
plt.title(r'Mandelbrot Set')

###############################################################################
###############################################################################
