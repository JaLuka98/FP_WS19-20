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

#################################
### Our statistical functions ###
#################################

def mittel(x, deltadof=1):
    return ufloat(np.mean(x),np.std(x,ddof=deltadof)/np.sqrt(len(x)))

########################################
### Our helper functions (fits, etc) ###
########################################

def lin(x,a,offset=0):
    return a*x + offset

def exp(a,b,c,x):
    return a*np.exp(b*x)+c

def lnint(T,yt,hr):
	array = np.array([])
	for t in T:
		array = np.append(array,[np.trapz(yt[T>=t],T[T>=t])/(yt[T == t]*hr)])
	return array

########################
### Reading the data ###
########################

temp1_2, current1_2 = np.genfromtxt('data/1_2.txt', unpack=True)
temp2, current2 = np.genfromtxt('data/2.txt', unpack=True)
# Temp in Celsius, Current in pA (10^-12)

##################################
### Heating rate interpolation ###
##################################

temps = [temp1_2, temp2] # Eine Liste der beiden Temperaturmessreihen (1_2 und 2)
currents = [current1_2, current2]
times = []
for t in temps:
    time = np.arange(start=0,stop=60*np.size(t),step=60)
    times.append(time)
# In times sind nun die Zeiten für die beiden Messreihen gespeichert.
b = []

for i in (0,1):
    params, cov = optimize.curve_fit(lambda x, a: lin(x, a, offset=(temps[i])[0]), times[i], temps[i])
    # The lambda thing just means to take the fit function
    # but each time with the specific offset
    params1 = correlated_values(params, cov)
    for p in params1:
        print('Heating rate b in K/s:',p)
        print('Heating rate b in K/min:',p*60)
        b.append(p) # Store the heating rates with error in K/s in the array b

###############################
### Plot the interpolations ###
###############################

labels = ['Heizrate 1_2', 'Heizrate 2']
names = ['build/interpol_1_2.pdf', 'build/interpol_2']

for i in (0,1):
    xlin = np.linspace(times[i][0], times[i][-1], 5000)
    plt.plot(times[i]/60, temps[i], 'b.', alpha=0.5, label = labels[i])
    plt.plot(xlin/60, lin(xlin, noms(b[i]), offset=temps[i][0]), 'r-', label='Ausgleichsfunktion')
    #plt.ylim(-10,45)
    #plt.xlim(-75,55)
    plt.grid()
    plt.xlabel(r'$t/$min')
    plt.ylabel(r'$T/$°C')
    plt.legend(loc='best')
    plt.savefig(names[i])
    plt.clf()

############################
### Untergrund bestimmen ###
############################

expfitx1=np.append(temp1_2[10:19],temp1_2[54:57])
expfity1=np.append(current1_2[10:19], current1_2[54:57])

print(expfitx1)
print(expfity1)

params1, covariance_matrix = optimize.curve_fit(exp, expfitx1, expfity1)
a, b, c = correlated_values(params1, covariance_matrix)
print('Fit für 1.2K/min:')
print('a=', a)
print('b=', b)
print('c=', c)

expfitx2=np.append(temp2[20:23],temp2[45:49])
expfity2=np.append(current2[20:23], current2[45:49])

params2, covariance_matrix = optimize.curve_fit(exp, expfitx2, expfity2)
a, b, c = correlated_values(params2, covariance_matrix)
print('Fit für 2K/min:')
print('a=', a)
print('b=', b)
print('c=', c)

expfitx=[expfitx1, expfitx2]
expfity=[expfity1, expfity2]
params=[params1, prams2]
############################
### Plot the actual data ###
############################

labels = ['Heizrate 1_2', 'Heizrate 2']
names = ['build/depolarisationskurve_1_2.pdf', 'build/depolarisationskurve_2']
x=np.linspace(-70,65,1000)

for i in (0,1):
    plt.plot(temps[i], currents[i], 'bx', alpha=0.75, label = labels[i])
    plt.plot(expfitx[i], expfity[i], 'rx', label = 'Werte für den Fit')
    plt.plot(x, exp(linspace, *params[i]),'g-', label='Ausgleichsrechnung', linewidth=0.5)
    #plt.ylim(-10,45)
    #plt.xlim(-75,55)
    plt.grid()
    plt.xlabel(r'$T/$°C')
    plt.ylabel(r'$I/$pA')
    plt.legend(loc='best')
    plt.savefig(names[i])
    plt.clf()






###########################################################################
### Celsius in Kelvin, um alle Probleme bei den Rechnungen zu vermeiden ###
###########################################################################

for i in (0,1):
    temps[i] += 273.15

##############################################
### Extract tau_0 from the maximum of I(T) ###
##############################################

# Beachte, dass wir für diese Berechnung W brauchen. Es wäre am konsistentesten,
# W aus der zweiten Methode zu verwenden.
# Bitte so coden, dass das W aus der zweiten Methode verwendet wird! Danke!

# Hier als Test: W ist 1eV
W = 1.602e-19

for i in (0,1):
    ### Die untenstehende Taktik ist statt argmax nötig, weil argmax
    ### nur den ersten index gibt, wo das erste mal der maximale wert kommt
    indices_Imax = np.argwhere(currents[i] == np.amax(currents[i])) # amax means maximum of array
    T_max = temps[i][indices_Imax]
    print('T_max beträgt', mittel(T_max, deltadof=0)) # Welche ddof sollen wir nehmen?
    tau_0 = kboltzmann*T_max**2/(W*b[i]) * unp.exp(-W/(kboltzmann*T_max))
    print(tau_0)
