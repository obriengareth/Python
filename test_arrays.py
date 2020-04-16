# -*- coding: utf-8 -*-
"""
 gareth 
"""


import time
import matplotlib . pyplot as plt
from matplotlib import cm
import math
from scipy import ndimage, misc
from scipy import signal
import numpy as np

#from __future__ import division
#cimport numpy as np
#cimport cython

def play_arrays():
    
    Image = np.zeros((no_cdp,ns));
    t1 = time.time()

    for s in range(0,NoS):    
        print(' running migration of shot',str(s))
        for r in range(0,NoR,1):
            for cdp in range(0,no_cdp):
                for timage in range(1,ns):
                        Image[cdp,timage]=cdp*timage;
    
    t2 = time.time()
    print('time 1',str(t2-t1))
    
 
    return


###############################################################################

print(' Running testing for array speed ')

no_cdp=300;
ns=251;
NoS=20;
NoR=80;
play_arrays()
  
    
#############################################################################
   