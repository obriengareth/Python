###############################################################################
# Exercise 7
###############################################################################

import matplotlib . pyplot as plt
import numpy as np
import time
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D  
from matplotlib import cm

data = pd.read_csv('cotopaxi1024x1024_15m', delimiter = r"\s+",header = None,names=['x','y','z'])

no=1;
N=1024;


# reshape data into matrices 
X=np.reshape(data.x.as_matrix(),(N,N))
Y=np.reshape(data.y.as_matrix(),(N,N))
Z=np.reshape(data.z.as_matrix(),(N,N))

###############################################################################
# pcolor plot
fig = plt.figure(num=no,figsize=(5, 5)) ;
plt.clf();
plt.pcolor(X,Y,Z);
# to complete exercise - add in recorders 
# to complete exercise - frame image ith axis labels etc

###############################################################################
# surface pot

fig = plt.figure(num=no+1,figsize=(5, 5)) ;
ax = fig.gca(projection='3d')

surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()
# to complete exercise - add in recorders 
# to complete exercise - frame image ith axis labels etc

###############################################################################
# contour plot
fig = plt.figure(num=no+2,figsize=(5, 5)) ;
plt.contour(X,Y,Z)
# to complete exercise - add in recorders 
# to complete exercise - frame image ith axis labels etc
###############################################################################
###############################################################################