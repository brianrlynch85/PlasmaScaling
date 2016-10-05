# Author: Brian Lynch
# Edited: 10/5/16
##############################################################################

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

##############################################################################
# Physical constants
echarge  = float("1.6022e-19")                 # Coloumbs
evjoule  = float("1.6022e-19")                 # eV = ev_to_joule Joules
eperm    = float("8.854e-12")                  # s^4 A^2 m^-3 kg^-1
m_elec   = float("9.109e-31")                  # kg
kb       = float("1.381e-23")                  # m^2 kg s^-2 K^-1

##############################################################################
# Functions   
def debyeL(kT,n,q):
   #print 'debyeL function'
   return math.sqrt(eperm * kT * evjoule / (n * q * q))
    
def larmorR(m,v,q,B):
   #print 'larmorR function'
   return m * v / (q * B)

def vth(kT,m):
   #print 'vth function'
   return math.sqrt(8.0 * kT * evjoule / (math.pi * m))

# Citation for this electron saturation current reduction factor
# author = D. Bohm
# title  = The Characteristics of Electrical Discharges in Magnetic Fields
# year   = 1949
def deltaCorr(rp,dl,beta):
   temp1 = np.sqrt(1.0 + dl * dl * np.divide(np.multiply(beta,beta),rp**2))
   temp2 = 1.0 / (1.0 + math.pi * rp * temp1 / (8.0 * dl))
   return temp2
    
###############################################################################
# Magnetic field parameters
B = 2.0                        # Tesla

# Electron parameters
kTe   = 2.0                    # eV,      Electron temperature
ne    = float("1.0e15")        # m^-3,    Electron density
etS   = vth(kTe,m_elec)        # m/s,     Electron thermal speed
eDL   = debyeL(kTe,ne,echarge) # m        Electron Debye Length

# Dust parameters
dustR = float("4.5e-6")        # m        Dust Radius

##############################################################################
# Start the plotting

# evenly sampled B in 0.1 T intervals
arrayB     = np.arange(0.05, 4.0, 0.05)
arrayLar   = larmorR(m_elec,etS,echarge,arrayB)
arrayBeta  = np.divide(dustR,arrayLar)
arrayBetaInv  = np.divide(1.0,arrayBeta)
arrayDelta = deltaCorr(dustR,eDL,arrayBeta)

# make the plot
plt.plot(arrayBetaInv,arrayDelta)
plt.grid(True, which='both')
plt.xlabel(r'$\beta^{-1}$')
plt.ylabel('$\delta$')
plt.title('Electron Current Reduction in Magnetic Field')
# increase the font size
mpl.rcParams.update({'font.size': 16})
          
testname = 'electron_collection_reduction.png'
          
plt.savefig(str(testname))
plt.show()
