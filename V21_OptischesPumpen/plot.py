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

def hyperbelfit(x,a,b,c):
    return a + b/(x-c)

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
muB=9.27*1e-24
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

hr = ['$f$/kHz', '$I_1$/A', '$I_2$/A', '$B_1$/µT' ,'$B_2$/µT']
m = np.zeros((10, 5))
m[:,0] = f*1e-3
m[:,1] = I_1
m[:,2] = I_2
m[:,3] = B_hor1*1e6
m[:,4] = B_hor2*1e6
t=matrix2latex(m, headerRow=hr, format='%.2f')
print(t)

#print('Horizontales Magnetfeld mit Isotop 1:', B_hor1)
#print('Horizontales Magnetfeld mit Isotop 2:', B_hor2)

params, covariance_matrix = optimize.curve_fit(linfit, f, B_hor1)
a, b = correlated_values(params, covariance_matrix)
print('Fit zu Isotop 1:')
print('a=', a) #a entspricht B/f
print('b=', b)
B_1=b

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
B_2=b

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
g_j=2.0023
I1=1/2*(g_j/g1-1)
I2=1/2*(g_j/g2-1)

print('Kernspin Isotop 1:', I1)
print('Kernspin Isotop 2:', I2)

#Quadratischer Zeeman
E1=g1*muB*B_1
E2=g2*muB*B_2

print('Linearer Zeeman 87:', E1)
print('Linearer Zeeman 85:', E2)

E1=g1**2*muB**2*B_1**2*(1-2)/4.53e-24
E2=g2**2*muB**2*B_2**2*(1-2)/2.01e-24

print('Quadratischer Zeeman zusätzlich 87:', E1)
print('Quadratischer Zeeman zusätzlich 85:', E2)

#Ansteigende Flanke

#Oszillationen
U, T_1, T_2 = np.genfromtxt('data/amplitude.txt', unpack=True)

T_1*=1e-3 #T in s
T_2*=1e-3

hr = ['$U$/V', '$T_1$/ms', '$T_2$/ms']
m = np.zeros((9, 3))
m[:,0] = U
m[:,1] = T_1*1e3
m[:,2] = T_2*1e3
t=matrix2latex(m, headerRow=hr, format='%.2f')
print(t)

params, covariance_matrix = optimize.curve_fit(hyperbelfit, U, T_1)
a, b, c = correlated_values(params, covariance_matrix)
print('Fit zu Isotop 1:')
print('a=', a) #a entspricht B/f
print('b=', b)
print('c=', c)

linspace=np.linspace(1, 11, 1000)
plt.plot(linspace, hyperbelfit(linspace, *params)*1e3, 'b-', label='Ausgleichsrechnung', linewidth=0.5)
plt.plot(U, T_1*1e3, 'rx', mew=0.5, label='Messwerte')
plt.xlabel(r'$U/$V')
plt.ylabel(r'$T$/ms')
plt.xlim(1, 11)
plt.tight_layout()
plt.legend()
plt.grid()
plt.savefig('build/oszi1.pdf')
plt.clf()

params, covariance_matrix = optimize.curve_fit(hyperbelfit, U, T_2)
a, b, c = correlated_values(params, covariance_matrix)
print('Fit zu Isotop 2:')
print('a=', a) #a entspricht B/f
print('b=', b)
print('c=', c)

linspace=np.linspace(1, 11, 1000)
plt.plot(linspace, hyperbelfit(linspace, *params)*1e3, 'b-', label='Ausgleichsrechnung', linewidth=0.5)
plt.plot(U, T_2*1e3, 'rx', mew=0.5, label='Messwerte')
plt.xlabel(r'$U/$V')
plt.ylabel(r'$T$/ms')
plt.xlim(1, 11)
plt.tight_layout()
plt.legend()
plt.grid()
plt.savefig('build/oszi2.pdf')
plt.clf()
