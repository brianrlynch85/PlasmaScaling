# Author: Brian Lynch
# Edited: 6/1/16
##############################################################################

import math
import numpy as np
import matplotlib.pyplot as plt

##############################################################################
# Physical constants

echarge  = float("1.6022e-19")                 # Coloumbs
evjoule  = float("1.6022e-19")                 # eV = ev_to_joule Joules
grav     = float("9.81")                       # m s^-2
eperm    = float("8.854e-12")                  # s^4 A^2 m^-3 kg^-1
TorrConv = 133.322                             # (N m^-2) / Torr
m_Ar     = float("6.633e-26")                  # kg
T_Ar     = 0.025 * evjoule                     # eV
kb       = float("1.381e-23")                  # m^2 kg s^-2 K^-1
v_Ar     = math.sqrt(T_Ar / (m_Ar))            # m s^-1
vt_Ar    = math.sqrt(8.0 / math.pi) * v_Ar     # m s^-1
r_Ar     = float("188.0e-12")                  # m
A_Ar     = r_Ar * r_Ar * math.pi               # m^2

##############################################################################
# Functions
def omega_p(n,m,q):
    return math.sqrt(n * q**2 / (m * eperm))
    
def omega_B(q,B,m):
    return q * B / m

def debL(kT,n,q):
    return math.sqrt(eperm * kT * evjoule / (n * q**2))

def vav(kT,m):
    return math.sqrt(8.0 * kT * evjoule / (m * math.pi))
    
def gyro_rad(m,v,q,B):
    return m * v / (q * B)    

##############################################################################
# evenly sampled B 0.1 T intervals
B         = np.arange(0.0, 4.0, 0.1)
mfp_array = np.empty(40)

#plot for a variety of ion temperatures
kT_sil_list = [0.025,0.05,0.10,0.5000]
plt.rc('lines', linewidth = 2)
plt.xticks(np.arange(0,5,0.5))
plt.yticks(np.arange(0,100,5))
plt.grid(True, which='both')

#Calculate the dust neutral collision mean free path for 200mTorr case
P_Ar     = 200.0 * TorrConv / 1000.0           # N m^-2
n_Ar     = P_Ar / (T_Ar)                       # # m^-3
mfp_Ar   = 1.0 / (A_Ar * n_Ar)                 # m
mfp_array.fill(100*mfp_Ar)
print('dust mfp [200mTorr] = ' + '{:3.2e}'.format(mfp_Ar) + ' [m]')
plt.plot(B,mfp_array,'g--',label='Neutral Collision Mean Free Path[200mTorr]')

#Calculate the dust neutral collision mean free path for 100mTorr case
P_Ar     = 100.0 * TorrConv / 1000.0           # N m^-2
n_Ar     = P_Ar / (T_Ar)                       # # m^-3
mfp_Ar   = 1.0 / (A_Ar * n_Ar)
mfp_array.fill(100*mfp_Ar)
print('dust mfp [200mTorr] = ' + '{:3.2e}'.format(mfp_Ar) + ' [m]')
plt.plot(B,mfp_array,'m--',label='Neutral Collision Mean Free Path[100mTorr]')

#Calculate the dust neutral collision mean free path for 10mTorr case
P_Ar     = 10.0 * TorrConv / 1000.0            # N m^-2
n_Ar     = P_Ar / (T_Ar)                       # # m^-3
mfp_Ar   = 1.0 / (A_Ar * n_Ar)
mfp_array.fill(100*mfp_Ar)
print('dust mfp [200mTorr] = ' + '{:3.2e}'.format(mfp_Ar) + ' [m]')
plt.plot(B,mfp_array,'k--',label='Neutral Collision Mean Free Path[10mTorr]')

#Plot for a variety of different ion temperatures
for i in range (0,4,1):
    v_itr = vav(kT_sil_list[i],m_Ar)
    plt.plot(B,100*gyro_rad(m_Ar,v_itr,echarge,B),label='T = '
             + str(round(kT_sil_list[i],3)) +
             '[eV], $v_{avg}$ = ' + str(round(v_itr,2)) + '[m/s]')
  
plt.axis([0.0,4.0,0.001,20])
plt.yscale('log')
plt.legend(loc=1,ncol=1,borderaxespad=0.,prop={'size':11})
plt.xlabel('B [T]')
plt.ylabel('gyro radius [cm]')
plt.title('Argon Ions')
          
testname = 'Argon_ion_gyro.png'
          
plt.savefig(str(testname))
plt.show()
