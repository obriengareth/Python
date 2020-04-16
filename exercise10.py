###############################################################################
# Exercise 10 FD soln
###############################################################################

import matplotlib . pyplot as plt
import numpy as np

###############################################################################
# Parameters
no=1;
Nx=90;		#Length X-direction
Nz=60;		# Length Z-direction
run_time = 5; 	# seconds
dt=0.05;	# time sampling
dx=1;		#% spatial sampling
dz=1;
Dxx=4;  # make these larger and see what happens
Dzz=1.0;
Dxz=0.005;
###############################################################################
# Initial Conditions
p = np.zeros((Nx,Nz), dtype=np.float32)
p[44:46,29:31]=1;  # pressure spike here

###############################################################################
#X = np.zeros((Nx,Nz), dtype=np.float32)
#Z = np.zeros((Nx,Nz), dtype=np.float32)

x = np.linspace(0, dx*Nx, Nx)
z = np.linspace(0, dz*Nz, Nz)
[Z,X] = np.meshgrid(z,x);

###############################################################################
dpdx = np.zeros((Nx,Nz), dtype=np.float32)
dp2dx2 = np.zeros((Nx,Nz), dtype=np.float32)
dp2dz2 = np.zeros((Nx,Nz), dtype=np.float32)
dp2dxdz = np.zeros((Nx,Nz), dtype=np.float32)

# using np.roll to offset arrays for spatial derivatives
# quicker that explicitly writing dp(i)/dx=(p(i+1)-p(i-1))/2/dx
Nt = int(run_time/dt)
time_step=0
for time_step in range(0,Nt):

    # 6th order second derivative
    dp2dz2 =  (  1/90*np.roll(p,3,axis=1)  - (3/20)*np.roll(p,2,axis=1)  + (3/2)*np.roll(p,1,axis=1) - (49/18)*p + 
                 1/90*np.roll(p,-3,axis=1) - (3/20)*np.roll(p,-2,axis=1) + (3/2)*np.roll(p,-1,axis=1) )/dz/dz;
    dp2dx2 =  (  1/90*np.roll(p,3,axis=0)  - (3/20)*np.roll(p,2,axis=0)  + (3/2)*np.roll(p,1,axis=0) - (49/18)*p + 
                 1/90*np.roll(p,-3,axis=0) - (3/20)*np.roll(p,-2,axis=0) + (3/2)*np.roll(p,-1,axis=0) )/dx/dx;

    # 6th order second derivative in x then y
    dpdx  = (-1/60*np.roll(p,3,axis=0)+3/20*np.roll(p,2,axis=0)-3/4*np.roll(p,1,axis=0)+ 
             3/4*np.roll(p,-1,axis=0)-3/20*np.roll(p,-2,axis=0)+1/60*np.roll(p,-1,axis=0))/dx;

    dp2dxdz  = (-1/60*np.roll(dpdx,3,axis=1)  + 3/20*np.roll(dpdx,2,axis=1)  - 3/4*np.roll(dpdx,1,axis=1) 
                 +3/4*np.roll(dpdx,-1,axis=1) - 3/20*np.roll(dpdx,-2,axis=1) + 1/60*np.roll(dpdx,-1,axis=1))/dz;

    # 2nd order time update
    p = p  + dt*( Dxx*dp2dx2 + Dzz*dp2dz2 + 2*Dxz*dp2dxdz );
    
###############################################################################
###############################################################################
    '''
    fig = plt.figure(num=no,figsize=(9, 6)) ;
    plt.clf();
    plt.pcolor(X,Z,p)#vmin=0,vmax=1);
    plt.xlabel('Length')
    plt.ylabel('Depth')
    plt.title(r'Pore Pressure Evolution')
    plt.pause(0.1);
    plt.show()
    '''
###############################################################################
###############################################################################

fig = plt.figure(num=no,figsize=(9, 6)) ;
plt.clf();
plt.pcolor(X,Z,p)#vmin=0,vmax=1);
plt.xlabel('Length')
plt.ylabel('Depth')
plt.title(r'Pore Pressure Evolution')
plt.pause(0.1);
plt.show()

###############################################################################
###############################################################################


