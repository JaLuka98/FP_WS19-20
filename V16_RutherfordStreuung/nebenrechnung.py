import numpy as np
from uncertainties import ufloat

# Luft

e = 1.602176487e-19
z = 2
n = 3.5e26
m_e = 0.511 # MeV
m_p = 938 # MeV
E = 5.49 # MeV
v = np.sqrt(2*E/(4*m_p))
eps_0 = 8.854187817e-12
I = 80.5 # eV

dE_by_ds = (4*np.pi*e*e*z*z*n)/(m_e*1e6*v*v*(4*np.pi*eps_0)**2) * np.log(2*m_e*1e6*v*v/I)
# now need rho
rho = 1.2041 # in kg per m^3
dE_by_ds *= rho
print(dE_by_ds*1e-6*1e-2)

# Gold, Foliendicke

b = ufloat(25.3, 1.1)
b_Folie = ufloat(16.61, 0.35)
DeltaE = E*(1-b_Folie/b)
DeltaE *= 1e6*e # nun in J
print('Energiedifferenz', DeltaE)
I = 790 # Anregungsenergie von Gold in eV
I = 790*e # Anregungsenergie nun in Joule


# Bestimmung von n:
rho = 19302 # in kg/m^3
M = 196.97*1e-3 # Molare Masse in kg/mol
N_A = 6.0221409e23
Z = 79

N = N_A*rho/M
print('Atomdichte von Gold: ',N)
n = N*Z
print('Elektronendichte von Gold',n)
m_e *= 1e6*e # m_e nun in J


Deltax = DeltaE * (m_e*v*v*(4*np.pi*eps_0)**2)/((4*np.pi*e**4*z*z*n)*np.log(2*m_e*v*v/I))

# Deltax in m
print('Deltas in micrometern:', Deltax*1e6)
