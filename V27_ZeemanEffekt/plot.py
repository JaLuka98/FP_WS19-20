import numpy as np
import scipy.stats
import scipy.optimize
import uncertainties as unc
import uncertainties.unumpy as unp
from scipy import optimize
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy import stats
from uncertainties import correlated_values
from matrix2latex import matrix2latex

import scipy.constants as const
h = const.h
c = const.speed_of_light
mu_B = const.value('Bohr magneton')



def lin(x,a,b):
    return a*x+b


def delta_lambda(delta_s, Delta_s, delta_lambda_D):
    return delta_s/(2*Delta_s)*delta_lambda_D


def g_J(lam_0, delta_lam, B):
    return h*c/(B*mu_B) * (1/lam_0 - 1/(lam_0+delta_lam))


I, B = np.genfromtxt('data/b.txt', unpack=True)
hr = ['$I$/A', '$B$/mT']
m = np.zeros((len(I), 2))
m[:,0] = I
m[:,1] = B
t=matrix2latex(m, headerRow=hr, format='%.1f')
print(t)
B*=1e-3
# I now in A and B in T

params, covariance_matrix = optimize.curve_fit(lin, I, B)
a, b = correlated_values(params, covariance_matrix)
print('a=', a)
print('b=', b)


def B_field(I): # returns the magnetic field with uncertainty
    return a*I+b


linspace=np.linspace(0, 25, 1000)
plt.plot(linspace, lin(linspace, *params), 'b-', label='Ausgleichsrechnung', linewidth=0.5)
plt.plot(I, B, 'rx', label='Messwerte')
plt.xlabel(r'$I/$A')
plt.ylabel(r'$B$/T')
plt.xlim(-1, 20)
plt.ylim(-0.1,1.2)
plt.tight_layout()
plt.legend()
plt.grid()
plt.savefig('build/b.pdf')
plt.clf()

##### Berechnung von g_J für die rote Linie

Delta_s, delta_s = np.genfromtxt('data/temp/pixel_rot.txt', unpack=True)

hr = ['$\Delta s$/px', '$\delta s$/px', '$\delta \lambda$/pm']
m = np.zeros((len(Delta_s), 3))
m[:,0] = Delta_s
m[:,1] = delta_s

Delta_s = unp.uarray(Delta_s, 5)
delta_s = unp.uarray(delta_s, 3)
print('Rote Linie:')
lam_0 = 633.8*1e-9
delta_lambda_D = 4.8912559474967786e-11
I = 9
delta_lam = delta_lambda(delta_s, Delta_s, delta_lambda_D)

m[:,2] = 1e12*unp.nominal_values(delta_lam)
t=matrix2latex(m, headerRow=hr, format='%.2f')
print(t)

delta_lam = sum(delta_lam)/len(delta_lam)
print(delta_lam)
print(g_J(lam_0, delta_lam, B_field(I)))

##### Berechnung von g_J für die blaue Linie,

Delta_s, delta_s_sigma, delta_s_pi = np.genfromtxt('data/temp/pixel_blau.txt', unpack=True)
Delta_s = unp.uarray(Delta_s, 5)
delta_s_sigma = unp.uarray(delta_s_sigma, 5)
print('Blaue Linie Sigma (0 Grad):')
lam_0 = 480*1e-9
delta_lambda_D = 2.6952020928905114e-11
I = 5.5
delta_lam = delta_lambda(delta_s_sigma, Delta_s, delta_lambda_D)
delta_lam = sum(delta_lam)/len(delta_lam)
print(delta_lam)
print(g_J(lam_0, delta_lam, B_field(I)))

delta_s_pi = unp.uarray(delta_s_pi, 5)
print('Blaue Linie Sigma (90 Grad):')
lam_0 = 480*1e-9
delta_lambda_D = 2.6952020928905114e-11
I = 16
delta_lam = delta_lambda(delta_s_pi, Delta_s, delta_lambda_D)
delta_lam = sum(delta_lam)/len(delta_lam)
print(delta_lam)
print(g_J(lam_0, delta_lam, B_field(I)))
