# -*- coding: utf-8 -*-
"""
SL2AG - gareth 9/1/18
"""
#from math import *
import matplotlib . pyplot as plt
import numpy as np
import math

#############################################################################

ntr=30;
ns=1000;

dz=5;
max_angle=60;
angle_bin = 2.0;
dhx = 12.5;
dhi = 12.5

#############################################################################
# fake cig gather - comes in from reveal
cig = np.zeros((ns,ntr,ntr));
#cig=cig + 0.075*(-1+2*np.random.rand(ns,ntr));
#cig[800:810,10:20]=10;
#cig[400:405,14:16]=10;

ns=np.size(cig,0)
ntr=np.size(cig,1)
#ntr=nn(2);
#ns=nn(1);
print('Traces',str(ntr),' samples',str(ns));


a = np.arange(6).reshape(2,3)
for x in np.nditer(a):
    print x



#############################################################################
# setting up needed arrays and sampling constants
pi=math.pi;
Lnh=ntr*2;
Lns=ns;
shift = ntr*dh/2 ;
freq = np.zeros((Lns,Lnh));
for iz in range (0,Lns):
    for ih in range (0,Lnh):
        freq[iz][ih]=ih/(Lnh*dh);

dkh = 1/Lnh/dh ;
Nangle = int( np.ceil( max_angle/angle_bin) + 1 );
       
#############################################################################
# doing the FK transform with the phase shift

fcig=np.fft.fft(cig,Lnh,1)
fcig=fcig*(np.exp( (1j)*2*pi*freq*shift));
fcig=np.fft.fft(fcig,Lns,0)

#############################################################################
# radial stack over the fk domain
facig = np.zeros((ns,Nangle),dtype=complex);
for ang in range(0,Nangle):
    for iz in range(0,Lns):
        kz = iz/Lns/dz - 1.0/Lns/dz;
        kh =  ( kz * np.tan( ang*angle_bin*(pi/180)) ); 
        ih =  int( 1 + np.floor(kh/dkh) );
        rem = kh/dkh - np.floor(kh/dkh);
        #print('ih= ',str(ih),'iz= ',str(iz));
        if(ih+1<Lnh and ih>0):
           facig[iz][ang] = facig[iz][ang] + ( fcig[iz][ih]*(1-rem) + fcig[iz][ih+1]*(rem) )

#############################################################################

acig=np.real(np.fft.irfft(facig,Lns,0))  
print('Done with',str(Nangle),'angles');

#############################################################################

#############################################################################
# plotting if needed
#sp_cig=np.log10(np.absolute(fcig));
#ph_cig=np.angle(fcig);
fig = plt.figure(figsize=(6, 6)) ;plt.pcolormesh(cig)
#fig = plt.figure(figsize=(3, 6)) ;plt.pcolormesh(sp_cig)
#fig = plt.figure(figsize=(3, 6)) ;plt.pcolormesh(ph_cig)
#sp_facig=np.flipud((np.absolute(facig)));
#fig = plt.figure(figsize=(6, 6)) ;plt.pcolormesh(sp_facig); 
fig = plt.figure(figsize=(6, 6)) ;plt.pcolormesh(acig)
#############################################################################
     