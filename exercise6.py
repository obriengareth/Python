##############################################################################
# Exercise 1 
#############################################################################

##############################################################################
# Exercise 1 
#############################################################################
import matplotlib . pyplot as plt
import numpy as np
import time


no=1;
N=50;


###############################################################################
t1 = time.time()
pi_calc =0 ;
for i in range(0,N):
    if(i%2==0):
        nn=1;
    else:
        nn=-1;
    pi_calc = pi_calc+ ( (nn)/(2*i+1))
pi_calc = 4*pi_calc

t2 = time.time()
print('result =',str(pi_calc),' finished in ',str(t2-t1),' seconds \n')
 
###############################################################################
t1 = time.time()
pi_est = np.zeros((N), dtype=np.float32)
for i in range(0,N):
    if(i%2==0):
        nn=1;
    else:
        nn=-1;
    pi_est[i] = ( (nn)/(2*i+1))
pi_est=pi_est*4;
res0 = np.cumsum(pi_est);
res1 = res0[N-1]; 
print('result =',str(res1),' finished in ',str(t2-t1),' seconds \n')

###############################################################################


fig = plt.figure(num=no,figsize=(5, 5)) ;
plt.clf();
#plt.subplot(311);
plt.plot(res0,'k')
plt.xlabel('i')
plt.ylabel('pi(i)')
plt.title(r'Conergence to $\pi $')
plt.grid()
