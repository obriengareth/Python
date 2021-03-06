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

#############################################################################
print(' Running Radon development code v1.0 for REVEAL plugin')

# gather size     
nt=30;
ns=500;

#############################################################################
# input parameters for gui

dt=0.004;
dp=1;
Np=50;
pmax=0.00025;
lin_para=2;


#############################################################################
# fake cig gather - comes in from reveal
cig = np.zeros((nt,ns));
#cig=cig + 0.075*(-1+2*np.random.rand(nt,ns));
#cig[0:20,200:210]=10;
#cig[0:16,400:405]=10;

nt=np.size(cig,0)
ns=np.size(cig,1)
    
print('Traces',str(nt),' samples',str(ns));

#############################################################################
# setting up needed arrays and sampling constants
#off = np.zeros((nt));
#pvals = np.zeros((Np));
  
off = np.linspace( 0, nt-1, nt ) 
off = off*100;
pvals = np.linspace( 0, pmax, Np ) ;
mm=np.amax(off);

vel=5000;
t0=100;
for h in range(0,nt): 
    it = int( np.sqrt( (t0*dt*t0*dt + off[h]*off[h]/vel/vel)  )/dt  );
    if(it>0 and it+1<ns):
        cig[h][it] =  cig[h][it] + 5 + 5*off[h]/mm; 
        cig[h][it+1] =  cig[h][it+1] - 5 - 5*off[h]/mm; 

vel=5000;
t0=400;
for h in range(0,nt): 
    it = int( np.sqrt( (t0*dt*t0*dt - off[h]*off[h]/vel/vel)  )/dt  );
    if(it>0 and it+1<ns):
        cig[h][it] =  cig[h][it] + 5 - 10*off[h]/mm; 
        cig[h][it+1] =  cig[h][it+1] - 5 + 10*off[h]/mm; 

vel=500000;
t0=200;
for h in range(0,nt): 
    it = int( np.sqrt( (t0*dt*t0*dt - off[h]*off[h]/vel/vel)  )/dt  );
    if(it>0 and it+1<ns):
        cig[h][it] =  cig[h][it] + 5 -4*off[h]/mm; 
        cig[h][it+1] =  cig[h][it+1] - 5 +4*off[h]/mm; 


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

Rcig = np.zeros((Np,ns));
for p in range(0,Np):
    for h in range(0,nt):
        for t in range(0,ns):
            t1 = (t*dt + np.power(off[h],lin_para)*pvals[p])/dt;
            it = int( np.floor( t1 ) );
            rem= t1-it;
            if(it>=0 and it+1<ns):
                Rcig[p][t] = Rcig[p][t] + ( cig[h][it]*(1-rem) + cig[h][it+1]*(rem) )
#Rcig=Rcig/nt;

fig = plt.figure(figsize=(6, 6)) ;plt.pcolormesh(np.transpose(Rcig))
plt.gca().invert_yaxis()
plt.colorbar();

#############################################################################
# slant stack inverse for radon
Icig = np.zeros((nt,ns));
for h in range(0,nt):
    for p in range(0,Np):
        for t in range(0,ns):
            t1 = (t*dt - np.power(off[h],lin_para)*pvals[p])/dt;
            it = int( np.floor( t1 ) );
            rem= t1-it;
            if(it>=0 and it+1<ns):
                Icig[h][t] = Icig[h][t] + ( Rcig[p][it]*(1-rem) + Rcig[p][it+1]*(rem) )
#Icig=Icig/Np;

fig = plt.figure(figsize=(6, 6)) ;plt.pcolormesh(np.transpose(Icig))
plt.gca().invert_yaxis()
plt.colorbar();

print(' Done ');

#############################################################################
fs=1/dt;
pi=np.pi;
freq=np.linspace( 0, fs, ns ) ;
fcig=np.fft.fft(cig,ns,1)
nsft=np.size(fcig,1)
fu = np.zeros((ns,Np),dtype=complex);
#I=np.identity(Np)
reg_no=1e-4
#for i in range (0,Np):
#    I[i][i] = I[i][i] + reg_no#*(-1 + 2*np.random.random(1));
#I=0.25*np.eye(Np,Np,k=-1) + np.eye(Np,Np,k=0) + 0.25*np.eye(Np,Np,k=1)

f=3; 
off1=np.power(off,lin_para);

for f in range(0,ns):
    L = np.zeros((nt,Np),dtype=complex);
    fd = np.zeros((nt),dtype=complex);
    
    for p in range(0,Np):
        for s in range(0,nt):
            L[s][p] = np.exp( -2*pi*1j*pvals[p]*off1[s]*freq[f]);

    M1 = np.dot(np.transpose(L),L) #+ I*reg_no#*(np.random.random(1)); 
    #M1i= np.linalg.inv(M1);
    M1i= np.linalg.pinv(M1,reg_no);
    M2 = np.dot(M1i,np.transpose(L))
    
    for k in range(0,nt):
        fd[k] = fcig[k][f]

    fu[f][:] = np.dot(M2,(fd))

RIcig = np.transpose(np.real(np.fft.irfft(fu,ns,0)))

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

fig = plt.figure(figsize=(6, 6)) ;
plt.pcolormesh(np.transpose(IIcig),vmin=-0.12,vmax=0.12)
plt.gca().invert_yaxis()
plt.colorbar();

fig = plt.figure(figsize=(6, 6)) ;
plt.pcolormesh(np.transpose(IIcig-cig),vmin=-0.12,vmax=0.12)
plt.gca().invert_yaxis()
plt.colorbar();
  
#############################################################################
# examples for prior code in case
#Rcig = np.zeros((ns,Nangle),dtype=complex);
#fcig=np.fft.fft(cig,Lnh,1)
#fcig=fcig*(np.exp( (1j)*2*pi*freq*shift));
#fcig=np.fft.fft(fcig,Lns,0)
#Rcig = np.zeros((ns,Nangle),dtype=complex);
#############################################################################
   