# Author: Brian Lynch
# Edited: 6/1/16

import math
import numpy as np
import matplotlib.pyplot as plt

##############################################################################
# Physical constants
echarge  = float("1.6022e-19") # Coloumbs
evjoule  = float("1.6022e-19") # eV = ev_to_joule Joules
eperm    = float("8.854e-12")  # s^4 A^2 m^-3 kg^-1

##############################################################################
# Functions   
def vav(kT,m):
    return math.sqrt(8.0 * kT * evjoule / (m * math.pi))
    
def gyro_rad(m,v,q,B):
    return m * v / (q * B)    

##############################################################################
# Simulation Parameters
r_d  = 0.25            # microns, Dust diameter
n_d  = 2.2             # g/cm^3,  Dust mass density
kTe  = 2.0             # eV,      Electron temperature
U_d  = 1.5 * kTe       # eV,      Dust Potential
ne   = float("1.0e15") # m^-3,    Electron density

##############################################################################
# Important Calculated Global Constants
debL       = math.sqrt(eperm * kTe * evjoule / (ne * (echarge**2))) # m
R_d        = r_d * float("1.0e-6")                                  # m
Vol        = 4.0 * math.pi * (R_d**3) / 3.0                         # m^3
den_kgperm = n_d * (10**3)                                          # kg/m^3
m_d        = den_kgperm * Vol                                       # kg
q_d        = 4.0 * math.pi * eperm * U_d * R_d * (1.0 + R_d / debL) # C

# Evenly sampled B 0.1 T intervals
B         = np.arange(0.0, 4.0, 0.1)

# Plot a variety of dust temperatures
kT_sil_list = [0.025,0.100,1.000,5.000]
plt.rc('lines', linewidth = 2)
plt.xticks(np.arange(0,5,0.5))
plt.yticks(np.arange(0,100,5))
plt.grid(True, which='both')

# Make the plot
for i in range (0,4,1):
    v_itr = vav(kT_sil_list[i],m_d)
    plt.plot(B,100*gyro_rad(m_d,v_itr,q_d,B),label='T = '
             + str(round(kT_sil_list[i],2)) +
             '[eV], $v_{avg}$ = ' + str(round(100*v_itr,2)) + '[cm/s]')
  
plt.yscale('log')
plt.legend(loc=1,ncol=1,borderaxespad=0.,prop={'size':11})
plt.xlabel('B [T]')
plt.ylabel('gyro radius [cm]')
plt.title('Silica Dust: r = ' + str(r_d) + ' [$\mu$m], ' + 'Density = '
          + str(den_kgperm/10**3) + ' [g/$cm^3$], '
          + 'Q = ' + '{:3.0f}'.format(q_d/echarge) + ' [#e]')
          
testname = 'SilicaDust_r_' + str(r_d) + 'Density_' + str(den_kgperm) + '.png'
          
plt.savefig(str(testname))
plt.show()
