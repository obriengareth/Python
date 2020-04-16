###############################################################################
# Exercise 11 Time series basics
###############################################################################

import matplotlib . pyplot as plt
import numpy as np
from scipy import signal


Nt=10000;
dt=60 ; # dt = 1 minute

A1=1;		freq1=1/3600;
A2=1.25;	freq2=1/3600/6;
A3=0.5;		freq3=1/3600/24;
A4=0.0	;	freq4=1/3600/24/31; # Note A4 is longer than signal if not zero 

t = np.zeros((Nt), dtype=np.float32)
s = np.zeros((Nt), dtype=np.float32)

for i in range (0,Nt):
    t[i]=i*dt;
    s[i]=A1*np.sin(2*np.pi*freq1*t[i]) + A2*np.sin(2*np.pi*freq2*t[i]) + A3*np.sin(2*np.pi*freq3*t[i]) +  A4*np.sin(2*np.pi*freq4*t[i]) ;
    
###############################################################################

fs=1/dt;
L=len(s);
freq= np.linspace(0, L, L)/(L*dt);

# Calculate spectrum using fft
fft_s=np.fft.fft(s);
psd=np.abs(fft_s);
ph=np.unwrap(np.angle(fft_s));


# filter  1
Fpass1=2e-4; Fpass2=4.0e-4;
order=4;
b, a = signal.butter(order, [Fpass1/(fs/2), Fpass2/(fs/2)],btype='band')
#b, a = signal.butter(order, [0.1,0.5],btype='band')
sf1 = signal.filtfilt(b, a, s)

# filter  2
Fpass1=3e-5; Fpass2=6.0e-5;
order=2;
b, a = signal.butter(order, [Fpass1/(fs/2), Fpass2/(fs/2)],btype='band')
#b, a = signal.butter(order, [0.1,0.5],btype='band')
sf2 = signal.filtfilt(b, a, s)

# filter  3
Fpass1=0.5e-5; Fpass2=2.0e-5;
order=2;
b, a = signal.butter(order, [Fpass1/(fs/2), Fpass2/(fs/2)],btype='band')
#b, a = signal.butter(order, [0.1,0.5],btype='band')
sf3 = signal.filtfilt(b, a, s)

###############################################################################
no=1;
fig = plt.figure(num=no,figsize=(9, 6)) ;
plt.clf()
plt.subplot(211)
plt.plot(t,s,'b');
plt.ylabel('Amplitude');
plt.xlabel('Time (s)');


plt.subplot(212);
plt.plot(freq,psd);
plt.axis([0, 0.001, 0, 6000]);
plt.ylabel('PSD');
plt.xlabel('Frequency (Hz)');

###############################################################################

fig = plt.figure(num=no+1,figsize=(9, 6)) ;
plt.clf()
plt.subplot(411)
plt.plot(t,s,'b');
plt.ylabel('Amplitude');
plt.xlabel('Time (s)');

plt.subplot(412);
plt.plot(t,sf1);
plt.ylabel('Amplitude');
plt.xlabel('Time (s)');
plt.grid()

plt.subplot(413);
plt.plot(t,sf2);
plt.ylabel('Amplitude');
plt.xlabel('Time (s)');
plt.grid()

plt.subplot(414);
plt.plot(t,sf3);
plt.ylabel('Amplitude');
plt.xlabel('Time (s)');
plt.grid()


###############################################################################
#Monthly Central England precipitation (mm). Daily automated values used after 1996.
#Wigley & Jones (J.Climatol.,1987), Gregory et al. (Int.J.Clim.,1991)
#Jones & Conway (Int.J.Climatol.,1997), Alexander & Jones (ASL,2001). Values may change after QC.
#YEAR   JAN   FEB   MAR   APR   MAY   JUN   JUL   AUG   SEP   OCT   NOV   DEC   ANN
import pandas as pd

data = pd.read_csv('England_precipitation_monthly.dat', delimiter = r"\s+",header = None,names=['YEAR',   'JAN',   'FEB',   'MAR',   'APR' ,  'MAY'  , 'JUN'  , 'JUL' ,  'AUG'  , 'SEP' ,  'OCT' ,  'NOV',   'DEC' ,  'ANN'])

Y=139;

dataM = data.drop(['YEAR','ANN'],axis=1)
# reshape data into matrices 
rain=np.reshape(dataM.as_matrix(),(Y*12))

###############################################################################

#rain=taper(rain,10);

###############################################################################

dt=1; # Month
fs=1/dt;
fs=1/dt;
L=len(rain);
freq= np.linspace(0, L, L)/(L*dt);

# Calculate spectrum using fft
fft_s=np.fft.fft(rain);
psd=np.abs(fft_s);
ph=np.unwrap(np.angle(fft_s));
###############################################################################

fig = plt.figure(num=no+2,figsize=(9, 6)) ;
plt.clf()
plt.subplot(411)
plt.plot(rain,'b');
plt.ylabel('Rain');
plt.xlabel('Time (Months)');

plt.subplot(212);
plt.plot(freq,psd/(np.max(psd)));
plt.axis([0,0.5, 0, 0.1]);
plt.ylabel('PSD');
plt.xlabel('Frequency (Hz)');

###############################################################################
