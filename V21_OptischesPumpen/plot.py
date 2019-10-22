import numpy as np
import scipy.optimize
import uncertainties as unc
import uncertainties.unumpy as unp
from scipy import optimize
import matplotlib.pyplot as plt
from uncertainties import ufloat
from scipy import stats
from uncertainties import correlated_values
from matrix2latex import matrix2latex

def linfit(x, a, b):
    return a*x + b


#Daten fuer die Helmholtz-Spulenpaare:

r_sweep=0.1639
N_sweep=11
r_vert=0.11735
N_vert=20

#Konstanten
mu0=4*np.pi*1e-7
e0=1.602*1e-19
me=9.11*1e-31
h=6.626*1e-34

#Berechnung des vertikalen Erdmagnetfeldes
I_vert= 0.195 #Ampere

B_vert=mu0*8*I_vert*N_vert/(np.sqrt(125)*r_vert)
print('Vertikales Magnetfeld: ', B_vert)


#Berechnung der Lande Faktoren und des horizontalen Erdmagnetfeldes
f, I_1, I_2 = np.genfromtxt('data/frequenz.txt', unpack=True)

f*=1e3 #f in Hz
#I_1*=0.1 #I in A       #Ohne diese Faktoren sind die Ergebnisse deutlich passender!!!
#I_2*=0.1


B_hor1=mu0*8/np.sqrt(125)*I_1*N_sweep/r_sweep
B_hor2=mu0*8/np.sqrt(125)*I_2*N_sweep/r_sweep

#print('Horizontales Magnetfeld mit Isotop 1:', B_hor1)
#print('Horizontales Magnetfeld mit Isotop 2:', B_hor2)

params, covariance_matrix = optimize.curve_fit(linfit, f, B_hor1)
a, b = correlated_values(params, covariance_matrix)
print('Fit zu Isotop 1:')
print('a=', a) #a entspricht B/f
print('b=', b)

linspace=np.linspace(-50000, 1050000, 1000)
plt.plot(linspace*1e-3, linfit(linspace, *params*1e6), 'b-', label='Ausgleichsrechnung', linewidth=0.5)
plt.plot(f*1e-3, B_hor1*1e6, 'rx', mew=0.5, label='Messwerte')
plt.xlabel(r'$f/$kHz')
plt.ylabel(r'$B$/µT')
plt.xlim(-50, 1050)
plt.tight_layout()
plt.legend()
plt.grid()
plt.savefig('build/frequenz1.pdf')
plt.clf()

g1=4*np.pi*me/(e0*a)
print('Lande Faktor Isotop 1:', g1)

params, covariance_matrix = optimize.curve_fit(linfit, f, B_hor2)
a, b = correlated_values(params, covariance_matrix)
print('Fit zu Isotop 2:')
print('a=', a) #a entspricht B/f
print('b=', b)

linspace=np.linspace(-50000, 1050000, 1000)
plt.plot(linspace*1e-3, linfit(linspace, *params*1e6), 'b-', label='Ausgleichsrechnung', linewidth=0.5)
plt.plot(f*1e-3, B_hor2*1e6, 'rx', mew=0.5, label='Messwerte')
plt.xlabel(r'$f/$kHz')
plt.ylabel(r'$B$/µT')
plt.xlim(-50, 1050)
plt.tight_layout()
plt.legend()
plt.grid()
plt.savefig('build/frequenz2.pdf')
plt.clf()

g2=4*np.pi*me/(e0*a)
print('Lande Faktor Isotop 2:', g2)

#Berechnung des Kernspins



#Isoptopenverhältnis


#Quadratischer Zeeman

#Ansteigende Flanke

#Oszillationen
