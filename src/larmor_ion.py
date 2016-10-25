# Author: Brian Lynch
# Edited: 6/1/16
##############################################################################

import math
import numpy as np
import matplotlib.pyplot as plt
import plasma_parameters as plasma

# Plasma parameters
T_Ar     = 0.025 * plasma.evjoule              # eV
v_Ar     = math.sqrt(T_Ar / plasma.m_Ar)            # m s^-1
vt_Ar    = math.sqrt(8.0 / math.pi) * v_Ar     # m s^-1

##############################################################################
# Functions
def omega_p(n,m,q):
    return math.sqrt(n * q**2 / (m * plasma.eperm))
    
def omega_B(q,B,m):
    return q * B / m

def debL(kT,n,q):
    return math.sqrt(plasma.eperm * kT * plasma.evjoule / (n * q**2))

def vav(kT,m):
    return math.sqrt(8.0 * kT * plasma.evjoule / (m * math.pi))
    
def gyro_rad(m,v,q,B):
    return m * v / (q * B)    

##############################################################################
# Begin the plotting

# Evenly sampled B 0.1 T intervals
B         = np.arange(0.0, 4.0, 0.1)
mfp_array = np.empty(40)

#Plot for a several different ion temperatures
kTs = np.array([0.025, 0.05, 0.10, 0.5000])
plt.rc('lines', linewidth = 2)
plt.xticks(np.arange(0,5,0.5)); plt.yticks(np.arange(0,100,5))
plt.grid(True, which='both')

#Neutral pressures and ion temperatures
P_Ar = np.array([100.0, 50.0, 5.0])                    #mTorr
n_Ar = np.divide(P_Ar * plasma.TorrConv / 1000.0, T_Ar)   # # m^-3
mfp_Ar = np.divide(1.0, plasma.A_Ar * n_Ar)
marks = ['r--', 'g--' ,'b--']

#Plot the mean free paths at several different pressures
for i in xrange(0, len(P_Ar), 1):
   mfp_array.fill(100*mfp_Ar[i])
   print('P = ' + '{:3.2e}'.format(P_Ar[i]) + ' [mTorr]' + ';' + 'dust mfp = ' + \
         '{:3.2e}'.format(mfp_Ar[i]) + ' [m]')
   plt.plot(B,mfp_array,marks[i],label='Neutral Collision Mean Free Path at ' + \
            '{:3.0f}'.format(P_Ar[i]) + '[mTorr]')

#Plot the gyro-radius for several different ion thermal speeds
for i in range (0,4,1):
    v_itr = vav(kTs[i],plasma.m_Ar)
    plt.plot(B,100*gyro_rad(plasma.m_Ar,v_itr,plasma.echarge,B),label='T = '
             + str(round(kTs[i],3)) +
             '[eV], $v_{avg}$ = ' + str(round(v_itr,2)) + '[m/s]')
 
plt.axis([B[0],B[-1],0.001,20])
plt.yscale('log')
plt.legend(loc=1, ncol=1, borderaxespad=0., prop={'size':11})
plt.xlabel('B [T]'); plt.ylabel('gyro radius [cm]'); plt.title('Argon Ions')
          
figname = 'Argon_ion_gyro.png'
          
plt.savefig(figname)
plt.show()

