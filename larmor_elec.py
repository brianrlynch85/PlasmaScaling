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
eperm    = float("8.854e-12")                  # s^4 A^2 m^-3 kg^-1
m_elec   = float("9.109e-31")                  # kg
kb       = float("1.381e-23")                  # m^2 kg s^-2 K^-1

##############################################################################
# Functions   
def vav(kT,m):
    return math.sqrt(8.0 * kT * evjoule / (m * math.pi))
    
def gyro_rad(m,v,q,B):
    return m * v / (q * B)  

def vth(kT,m):
    return math.sqrt(2.0 * kT * evjoule / m)
    
# Evenly sampled B 0.1 T intervals
B = np.arange(0.0, 4.0, 0.1)

# Plot for a variety of electron temperatures
kT_sil_list = [0.1,0.5,2.0,4.0]
plt.rc('lines', linewidth = 2)
plt.xticks(np.arange(0,5,0.5))
plt.yticks(np.arange(0,10,5))
plt.grid(True, which='both')

# Make the plot
for i in range (0,4,1):
    v_itr = vav(kT_sil_list[i],m_elec)
    plt.plot(B,100*gyro_rad(m_elec,v_itr,echarge,B),label='T = '
             + str(round(kT_sil_list[i],3)) +
             '[eV], $v_{avg}$ = ' + str(round(v_itr,2)) + '[m/s]')
  
plt.yscale('log')
plt.legend(loc=1,ncol=1,borderaxespad=0.,prop={'size':11})
plt.xlabel('B [T]')
plt.ylabel('gyro radius [cm]')
plt.title('Electrons')
          
testname = 'electron_gyro.png'
          
plt.savefig(str(testname))
plt.show()
