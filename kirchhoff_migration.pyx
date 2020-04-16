# -*- coding: utf-8 -*-
"""
LSM 1 - gareth 
"""

#import os
#clear = lambda: os.system('cls')  # On Windows System
#clear()


import time
import matplotlib . pyplot as plt
from matplotlib import cm
import math
from scipy import ndimage, misc
from scipy import signal
import numpy as np

#from __future__ import division
#cimport numpy as np
cimport cython

def play_arrays():
    
    Image = np.zeros((no_cdp,ns));
    cdef ff
    t1 = time.time()

    for s in range(0,NoS):    
        print(' running migration of shot',str(s))
        for r in range(0,NoR,1):
            for cdp in range(0,no_cdp):
                for timage in range(1,ns):
                        Image[cdp,timage]=1;
    
    t2 = time.time()
    print('time 1',str(t2-t1))
    
    t1 = time.time()

    for s in range(0,NoS):    
        print(' running migration of shot',str(s))
        for r in range(0,NoR,1):
            for x in np.nditer(Image):    
                Image=1;
    
    t2 = time.time()
    print('time 2',str(t2-t1))
    
    
    return
"""
def simple_forward_migration(trace,cdpx,recx,srcx,no_cdp,aperture,vel,dx):
    
    Image = np.zeros((no_cdp,ns));
    dipA=45;
    
                  
    for s in range(5,8):    
        print(' running migration of shot',str(s))
        for r in range(0,NoR,1):
            for cdp in range(0,no_cdp):
                for timage in range(1,ns):
                    off = np.abs(cdpx[s][r]-cdp*dx);
                    if( off < aperture ):
                        if( np.arctan( off/2/(timage*dt*1000)*3.141/180.0 )<dipA ):
                            ti = timage*timage*dt*dt;
                            ss=(srcx[s][r]-cdp*dx)*(srcx[s][r]-cdp*dx);
                            rr=(recx[s][r]-cdp*dx)*(recx[s][r]-cdp*dx);
                            Td = np.sqrt( ti/4 + ss/vel/vel ) + np.sqrt( ti/4 + rr/vel/vel );
                            Tdi = int( np.floor(Td/dt) );
                            rem = (Td/dt)-np.floor(Td/dt);
                            if( Tdi > 0 and Tdi < ns-1 ):
                                Image[cdp][timage] = Image[cdp][timage] + (1-rem)*trace[s][r][Tdi];
                                Image[cdp][timage] = Image[cdp][timage] + rem*trace[s][r][Tdi+1];    
                         
    return Image
         
"""
###############################################################################
print(' Running LSM code v1.0 for testing ')

# gather size     
trace_file='trace_data1';
dt = 0.004;
vel = 2000;
dx=12.5 ;
no_cdp=300;

Nit=2;
cg_scale = 0.015;

xv = 60;
zv = 128;

s2n=0.02;

aperture=3000;
dipA =45;
    
###############################################################################
"""
fname = "%s.npy" % (trace_file)
trace = np.load(fname);
fname = "%s_CDPX.npy" % (trace_file)
cdpx = np.load(fname);
fname = "%s_SRCX.npy" % (trace_file)
srcx = np.load(fname);
fname = "%s_RECX.npy" % (trace_file)
recx = np.load(fname);


NoS = np.size(trace,0)
NoR = np.size(trace,1)
ns = np.size(trace,2)

trace = trace /np.max(np.max(trace));
trace = trace + s2n*(-1+2*np.random.rand(NoS,NoR,ns));

count = 0;
data = np.zeros((NoS*NoR,ns));
for s in range(0, NoS):
    for r in range(0, NoR):
        data[count][:] = trace[s][r][:];
        count = count + 1 ;
"""
        
"""
fig = plt.figure(figsize=(16, 8)) ;
cmap = plt.get_cmap('PiYG')
cmap = plt.get_cmap('gray')
plt.pcolormesh(np.transpose(data),cmap=cmap)
plt.gca().invert_yaxis()
plt.colorbar();
"""
"""
fig = plt.figure(figsize=(16, 8)) ;
cmap = plt.get_cmap('gray')
plt.pcolormesh(np.transpose(trace[5][:][:]),cmap=cmap)
plt.gca().invert_yaxis()
plt.colorbar();
"""
print('running migration...');
#Image = simple_forward_migration(trace,cdpx,recx,srcx,no_cdp,aperture,vel,dx)
 
play_arrays()
      
"""
fig = plt.figure(figsize=(16, 8)) ;
cmap = plt.get_cmap('gray')
plt.pcolormesh(np.transpose(Image),cmap=cmap)
plt.gca().invert_yaxis()
plt.colorbar();
"""
  
###############################################################################
#np.save(trace_file, trace)    # 

#fname = "%s_SRCX" % (trace_file)
#np.save(fname, srcx)    # 

#fname = "%s_RECX" % (trace_file)
#np.save(fname, recx)    # 

#fname = "%s_CDPX" % (trace_file)
#np.save(fname, cdpx)    # 

#print(' Done ');

#############################################################################
#M1i= np.linalg.pinv(M1,reg_no);
    
#############################################################################
   