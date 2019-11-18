import numpy as np
import scipy.optimize
import uncertainties as unc
import uncertainties.unumpy as unp
from scipy import optimize
import scipy.constants as const
kboltzmann = const.k
import matplotlib.pyplot as plt
from uncertainties import ufloat
from uncertainties import correlated_values, correlation_matrix
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
from matrix2latex import matrix2latex

k_B=1.3806485e-23
e=1.602e-19

#################################
### Our statistical functions ###
#################################

def mittel(x, deltadof=1):
    return ufloat(np.mean(x),np.std(x,ddof=deltadof)/np.sqrt(len(x)))

########################################
### Our helper functions (fits, etc) ###
########################################

def linfit(x,a,b):
    return a*x+b

def hyp(x,a,b):
    return a/x+b

def lin(x,a,offset=0):
    return a*x + offset

def exp(x,a,b,c):
    return a*np.exp(x*b)+c

def lnint(T,current,heatingRate):
	array = np.array([])
	for t in T:
		array = np.append(array,[np.trapz(current[T>=t],T[T>=t])/(current[T == t]*heatingRate)])
	return array

########################
### Reading the data ###
########################

temp1_2, current1_2 = np.genfromtxt('data/1_2.txt', unpack=True)
temp2, current2 = np.genfromtxt('data/2.txt', unpack=True)
# Temp in Celsius, Current in pA (10^-12)
temps = [temp1_2, temp2] # Eine Liste der beiden Temperaturmessreihen (1_2 und 2)
currents = [current1_2, current2]

############################
### Untergrund bestimmen ###   Vielleicht sind hier die Werte für i einfach so klein dass er da nichts findet... ist mir erst zu spät eingefallen....
############################

expfitx1=np.append(temp1_2[5:19],temp1_2[55:72])
expfity1=np.append(current1_2[5:19], current1_2[55:72])

params1, covariance_matrix = optimize.curve_fit(exp, expfitx1, expfity1)
a,b,c= correlated_values(params1, covariance_matrix)
print('Fit für 1.2K/min:')
print('a=', a)
print('b=', b)
print('c=', c)

expfitx2=np.append(temp2[11:17],temp2[37:51])
expfity2=np.append(current2[11:17], current2[37:51])

params2, covariance_matrix = optimize.curve_fit(exp, expfitx2, expfity2)
a, b, c= correlated_values(params2, covariance_matrix)
print('Fit für 2K/min:')
print('a=', a)
print('b=', b)
print('c=', c)

expfitx=[expfitx1, expfitx2]
expfity=[expfity1, expfity2]
params=[params1, params2]
############################
### Plot the actual data ###
############################

labels = ['Niedigere Heizrate', 'Höhere Heizrate']
names = ['build/depolarisationskurve_1_2.pdf', 'build/depolarisationskurve_2.pdf']
x=np.linspace(-100,100,1000)
xlims = [-70,65]
ylims = [-12.5,17.5]

for i in (0,1):
    plt.plot(temps[i], currents[i], 'bx', alpha=0.75, mew=0.5, label = labels[i])
    plt.plot(expfitx[i], expfity[i], 'rx', label = 'Werte für den Fit')
    plt.plot(x, exp(x, *params[i]),'g-', label='Ausgleichsrechnung', linewidth=0.5)
    plt.plot(temps[i], currents[i]-exp(temps[i],*params[i]), 'gx', mew=0.5,  label='Messwerte mit Untergrundkorrektur')
    plt.ylim(ylims)
    plt.xlim(xlims)
    plt.grid()
    plt.xlabel(r'$T/$°C')
    plt.ylabel(r'$I/$pA')
    plt.legend(loc='best')
    plt.savefig(names[i])
    plt.clf()
