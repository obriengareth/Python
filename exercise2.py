##############################################################################
# Exercise 2
#############################################################################
import numpy as np

N=5;

R = np.random.rand(N,N)

detR = np.linalg.det(R)

invR = np.linalg.inv(R)

#whos

np.matmul(R,invR)

R * invR

A=np.array([1,2,3,1,2,3,1,2,3]).reshape(3,3)

detA = np.linalg.det(A)

invA = np.linalg.inv(A)

#s=print('mean of odd ',np.mean(X[1,:]));
###############################################################################
