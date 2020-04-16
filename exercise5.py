##############################################################################
# Exercise 5
#############################################################################

import numpy as np
#import matplotlib.pyplot as plt

##############################################################################
# function definition
def quadratic_roots(a,b,c):
    
    x = np.zeros((2))
    x[0] = (-b + np.sqrt(b*b-4*a*c))/2/a
    x[1] = (-b - np.sqrt(b*b-4*a*c))/2/a
      
    return x

##############################################################################
# inputs  - build flag option/ IO option or command line 
a=1
b=-2
c=-2

x=quadratic_roots(a,b,c)
print('Solution is ',x)
#############################################################################
#############################################################################