##############################################################################
# Exercise 3
#############################################################################
import matplotlib . pyplot as plt
import numpy as np
import time


no=1;
x_min=-2;
x_max=2;
x_inc = 0.5;


###############################################################################
N=int((x_max-x_min)/x_inc)
X=np.linspace(x_min, x_max, num=N)
Y=X*X -X - 2;
###############################################################################

fig = plt.figure(num=no,figsize=(5, 5)) ;
plt.clf();
plt.plot(X,Y,'k')
buf = "Function" 
plt.xlabel('X')
plt.ylabel('f(X)')
plt.title(buf)


###############################################################################
###############################################################################
dYdX=X*0;
dYdX = np.diff(Y)/x_inc; 
# 
dX = np.zeros((len(X)-1), dtype=np.float32)
for i in range(0,len(X)-1):
    dX[i] = (X[i] + X[i+1])*0.5
    
intXY = np.cumsum(Y)*x_inc;

fig = plt.figure(num=no+1,figsize=(5, 5)) ;
plt.clf();
plt.subplot(311);
plt.plot(X,Y,'k')
buf = "Function" 
plt.plot(dX,dYdX,'b')
plt.plot(X,intXY,'r')
plt.xlabel('X')
plt.ylabel('f(X)')
plt.title(buf)
plt.grid()

plt.clf();
plt.subplot(312);
buf = "Function" 
plt.plot(dX,dYdX,'b')
plt.xlabel('X')
plt.ylabel('df(X)')
plt.title(buf)

plt.clf();
plt.subplot(313);
buf = "Function" 
plt.plot(X,intXY,'r')
plt.xlabel('X')
plt.ylabel('Int f(X)')
plt.title(buf)

###############################################################################
###############################################################################