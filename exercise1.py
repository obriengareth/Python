##############################################################################
# Exercise 1 
#############################################################################
import numpy as np

# Write even numbers into first row
Nx=2;
Ny=100;

np.arange(0, Ny, 2)
X = np.zeros((Nx,Ny), dtype=np.float32)
# Write even numbers into first row
a = np.arange(0, 2*Ny, 2)
# Write odd numbers into first row
b = np.arange(1, 2*Ny+1, 2)
X[0,:]=a;
X[1,:]=b;

# Display array
X

# Mean of odd numbers
np.mean(X[1,:])

# Mean of even numbers
np.mean(X[0,:])

#% Mean of both rows
np.mean(X)


#% Sum of odd numbers
np.sum(X[1,:])

#% Sum of even numbers
np.sum(X[0,:])

#% sum of both rows
np.sum(X)

#% Transpose
np.transpose(X)

###############################################################################
# Alternative way of displaying text
s=print('mean of odd ',np.mean(X[1,:]));

###############################################################################
