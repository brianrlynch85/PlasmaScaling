import math

##############################################################################
# Physical constants
echarge  = float("1.6022e-19")                 # Coloumbs
evjoule  = float("1.6022e-19")                 # eV = ev_to_joule Joules
grav     = float("9.81")                       # m s^-2
eperm    = float("8.854e-12")                  # s^4 A^2 m^-3 kg^-1
kb       = float("1.381e-23")                  # m^2 kg s^-2 K^-1
TorrConv = 133.322                             # (N m^-2) / Torr
m_Ar     = float("6.633e-26")                  # kg
r_Ar     = float("188.0e-12")                  # m
m_He     = float("6.6464e-27")                 # kg
r_He     = float("140.00e-12")                 # m

##############################################################################
# Neutral gas Parameters
m_NG     = m_Ar                                        # kg
r_NG     = r_Ar                                        # m
T_NG_n   = 0.025 * evjoule                             # eV
v_NG_n   = math.sqrt(T_NG_n / (m_NG))                  # m s^-1
vt_NG_n  = math.sqrt(8.0 * T_NG_n / (math.pi * m_NG))  # m s^-1
P_NG_n   = 10.0 * TorrConv / 1000.0                    # N m^-2
n_NG_n   = P_NG_n / (T_NG_n)                           # # m^-3
A_NG     = r_NG * r_NG * math.pi                       # m^2

##############################################################################
# Dust Parameters
r_d     = 4.0            # microns, Dust radius
rho_d   = 2.0             # g/cm^3,  Dust density
B       = 0.0             # Tesla,   Magnetic field strength
kTe     = 2.0             # eV,      Electron temperature
U_d     = 2.5 * kTe       # eV,      Dust Potential
ne      = float("5.0e15") # m^-3,    Electron density
delta   = 1.41            # s^-1,    Epstein drag coefficient
n_i     = ne              # # m^-3   Ion # Density
vt_NG_i = vt_NG_n         # m/s,     Ion Thermal Speed
u_i     = 35.0            # m/s,     Drift speed
phif    = 40.0            # V,       Plasma Potential

##############################################################################
# Calculated Parameters
debL       = math.sqrt(eperm * kTe * evjoule / (ne * (echarge**2))) # m
R_d        = r_d * float("1.0e-6")                                  # m
Vol        = 4.0 * math.pi * (R_d**3) / 3.0                         # m^3
den_kgperm = rho_d * (10**3)                                        # kg/m^3
m_d        = den_kgperm * Vol                                       # kg
q_d        = 550 * echarge#4.0 * math.pi * eperm * U_d * R_d * (1.0 + R_d / debL) # C
A_d        = math.pi * R_d**2
vs_i       = math.sqrt(u_i**2 + vt_NG_i**2)
b_c        = (1.0 - 2.0 * echarge * U_d * evjoule / (m_NG * vs_i**2))**2
b_c        = math.sqrt(b_c) * R_d
sigma_coll = math.pi * b_c**2
b_pid2     = q_d * echarge / (4.0 * math.pi * eperm * m_NG * vs_i**2)
sigma_E    = 2.0 * math.pi * (b_pid2**2)
sigma_E   *= math.log((debL**2 + b_pid2**2)/(b_c**2 + b_pid2**2))
