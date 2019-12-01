import numpy as np

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
