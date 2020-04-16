# -*- coding: utf-8 -*-
"""
radon dev 1 - gareth 12/01/18
"""
#from math import *
import matplotlib . pyplot as plt
import numpy as np
from matplotlib import cm
import math
from scipy import ndimage, misc
#import numarray

#############################################################################
print(' Running Radon development code v1.0 for REVEAL plugin')

# gather size     
nt=101;
ns=501;

#############################################################################
# input parameters for gui

Np=201;
pmax=0.008;
lin_para=1;
reg_no=0.05
dt=0.005;
dh=4;

#############################################################################
# fake cig gather - comes in from reveal
cig = np.zeros((nt,ns));

nt=np.size(cig,0)
ns=np.size(cig,1)
    
print('Traces',str(nt),' samples',str(ns));

#############################################################################
# setting up needed arrays and sampling constants
 
off = np.linspace( 0, nt-1, nt ) 
off = off*dh-((nt-1)*dh/2);
pvals = np.linspace( 0, pmax, Np ) ;
mm=np.amax(off);

vel=800;
t0=100;
for h in range(0,50): 
    it = int( ( (t0*dt + off[h]/vel)  )/dt  );
    if(it>0 and it+1<ns):
        cig[h][it] =  cig[h][it] + 5 ; 
vel=400;
t0=100;
for h in range(50,nt): 
    it = int( ( (t0*dt + off[h]/vel)  )/dt  );
    if(it>0 and it+1<ns):
        cig[h][it] =  cig[h][it] + 5 ; 

vel=300;
t0=300;
for h in range(0,nt): 
    it = int( ( (t0*dt + off[h]/vel)  )/dt  );
    if(it>0 and it+1<ns):
        cig[h][it] =  cig[h][it] + 5 ; 

cig=cig + 0.0*(-1+2*np.random.rand(nt,ns));
cig = ndimage.gaussian_laplace(cig, sigma=[0,4])
fig = plt.figure(figsize=(6, 6)) ;plt.pcolormesh(np.transpose(cig))
plt.gca().invert_yaxis();
plt.colorbar();

#############################################################################
# slant stack summation for Radon
if(lin_para==2):
    pmax=pmax*pmax;
pvals = np.linspace( -pmax, pmax, Np ) ;
off1=np.power(off,lin_para);

#############################################################################

Rcig = np.zeros((Np,ns));
#for p in range(0,Np):
#    for h in range(0,nt):
#        for t in range(0,ns):
#            t1 = (t*dt + np.power(off1[h],lin_para)*pvals[p])/dt;
#            it = int( np.floor( t1 ) );
#            rem= t1-it;
#            if(it>=0 and it+1<ns):
#                Rcig[p][t] = Rcig[p][t] + ( cig[h][it]*(1-rem) + cig[h][it+1]*(rem) )

#fig = plt.figure(figsize=(6, 6)) ;plt.pcolormesh(np.transpose(Rcig))
#plt.gca().invert_yaxis()
#plt.colorbar();

#############################################################################
fs=1/dt;
pi=np.pi;
freq=np.linspace( 0, fs, ns ) ;
fcig=np.fft.fft(cig,ns,1)
nsft=np.size(fcig,1)
fu = np.zeros((ns,Np),dtype=complex);

for f in range(0,ns):
    L = np.zeros((nt,Np),dtype=complex);
    fd = np.zeros((nt),dtype=complex);
    
    for p in range(0,Np):
        for s in range(0,nt):
            L[s][p] = np.exp( -2*pi*1j*pvals[p]*off1[s]*freq[f]);

    M1 = np.dot(np.transpose(L),L) 
    M1i= np.linalg.pinv(M1,reg_no);
   # M1i= np.linalg.pinv(M1);
    M2 = np.dot(M1i,np.transpose(L))
    
    for k in range(0,nt):
        fd[k] = fcig[k][f]

    fu[f][:] = np.dot(M2,(fd))

RIcig = np.transpose(np.real(np.fft.irfft(fu,ns,0)))
#AA=RIcig[50:151,:]

fig = plt.figure(figsize=(6, 6)) ;plt.pcolormesh(np.transpose(RIcig))
plt.gca().invert_yaxis()
plt.colorbar();

#############################################################################
DF = np.zeros((ns,nt),dtype=complex);
L = np.zeros((nt,Np),dtype=complex);
for f in range(0,ns):
  
    for p in range(0,Np):
        for s in range(0,nt):
            L[s][p] = np.exp( -2*pi*1j*pvals[p]*off1[s]*freq[f]);

    DF[f][:] = np.dot(L,(fu[f][:]))

IIcig = np.transpose(np.real(np.fft.ifft(DF,ns,0)))

#fig = plt.figure(figsize=(6, 6)) ;
#plt.pcolormesh(np.transpose(IIcig),vmin=-0.12,vmax=0.12)
#plt.gca().invert_yaxis()
#plt.colorbar();

#fig = plt.figure(figsize=(6, 6)) ;
#plt.pcolormesh(np.transpose(IIcig-cig),vmin=-0.12,vmax=0.12)
#plt.gca().invert_yaxis()
#plt.colorbar();
  
#############################################################################
 