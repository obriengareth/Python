##############################################################################
# Exercise 4
#############################################################################

import numpy as np
import matplotlib.pyplot as plt



# look at profiler options as N increases (if possible)
no=1;
N=1000;
no_bins=200;
mu, sigma = 100, 15

#############################################################################
#############################################################################
# look at seeding effects
# np.random.seed(562560) 

# scale between -1 and 1
R = -1 + 2*np.random.rand(N,1)
# different distribution
R = mu + sigma * np.random.randn(N,1)

# Mean of odd numbers
s=print('mean of the distribution is ',np.mean(R));

# the histogram of the data
n, bins, patches = plt.hist(R, no_bins, normed=1, facecolor='g', alpha=0.75)

fig = plt.figure(num=no,figsize=(5, 5)) ;

plt.xlabel('R')
plt.ylabel('Probability')
plt.title('Histogram')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
#plt.axis([40, 160, 0, 0.03])
plt.grid(True)
plt.show()

#############################################################################
#############################################################################
