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
    return a*np.exp(b*x)+c

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
        heatRate = b

###############################
### Plot the interpolations ###
###############################

labels = ['Heizrate 1_2', 'Heizrate 2']
names = ['build/interpol_1_2.pdf', 'build/interpol_2']

for i in (0,1):
    xlin = np.linspace(times[i][0], times[i][-1], 5000)
    plt.plot(times[i]/60, temps[i], 'b.', alpha=0.5, label = labels[i])
    plt.plot(xlin/60, lin(xlin, noms(heatRate[i]), offset=temps[i][0]), 'r-', label='Ausgleichsfunktion')
    #plt.ylim(-10,45)
    #plt.xlim(-75,55)
    plt.grid()
    plt.xlabel(r'$t/$min')
    plt.ylabel(r'$T/$°C')
    plt.legend(loc='best')
    plt.savefig(names[i])
    plt.clf()

############################
### Untergrund bestimmen ###   Vielleicht sind hier die Werte für i einfach so klein dass er da nichts findet... ist mir erst zu spät eingefallen....
############################

expfitx1=np.append(temp1_2[10:19],temp1_2[55:58])
expfity1=np.append(current1_2[10:19], current1_2[55:58])

params1, covariance_matrix = optimize.curve_fit(exp, expfitx1, expfity1, p0=[0.01,0.01,0])
a, b, c = correlated_values(params1, covariance_matrix)
print('Fit für 1.2K/min:')
print('a=', a)
print('b=', b)
print('c=', c)

expfitx2=np.append(temp2[12:15],temp2[37:41])
expfity2=np.append(current2[12:15], current2[37:41])

params2, covariance_matrix = optimize.curve_fit(exp, expfitx2, expfity2, p0=[0.01,0.01,0])
a, b, c = correlated_values(params2, covariance_matrix)
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

labels = ['Heizrate 1_2', 'Heizrate 2']
names = ['build/depolarisationskurve_1_2.pdf', 'build/depolarisationskurve_2']
x=np.linspace(-70,65,1000)

for i in (0,1):
    plt.plot(temps[i], currents[i], 'bx', alpha=0.75, mew=0.5, label = labels[i])
    plt.plot(expfitx[i], expfity[i], 'rx', label = 'Werte für den Fit')
    plt.plot(x, exp(x, *params[i]),'g-', label='Ausgleichsrechnung', linewidth=0.5)
    plt.plot(temps[i], currents[i]-exp(temps[i],*params[i]), 'gx', mew=0.5,  label='Messwerte mit Untergrundkorrektur')
    #plt.ylim(-10,45)
    #plt.xlim(-75,55)
    plt.grid()
    plt.xlabel(r'$T/$°C')
    plt.ylabel(r'$I/$pA')
    plt.legend(loc='best')
    plt.savefig(names[i])
    plt.clf()

hr = ['$T_1$/°C', '$I_1$/pA', '$I_{1,bereinigt}$/pA', '$T_2$/°C', '$I_2$/pA', '$I_{2,bereinigt}$/pA',]
m = np.zeros((108, 6))
m[:,0] = temp1_2
m[:,1] = current1_2
m[:,2] = current1_2-exp(temp1_2, *params1)
m[:,3] = np.append(temp2, np.zeros((1, 36)))
m[:,4] = np.append(current2, np.zeros((1, 36)))
m[:,5] = np.append(current2-exp(temp2, *params2), np.zeros((1, 36)))
t=matrix2latex(m, headerRow=hr, format='%.2f')
print(t)


current1_2-=exp(temp1_2, *params1)    #Untergrund Abziehen
current2-=exp(temp2, *params2)



###########################################################################
### Celsius in Kelvin, um alle Probleme bei den Rechnungen zu vermeiden ###
###########################################################################
for i in (0,1):
    temps[i] += 273.15

# Warum doppelt? ;)
#temp1_2+=273.15
#temp2+=273.15

#current1_2*=1e12  #Umrechnung in A. Vielleicht besser wegen Einheiten?
#current2*=1e12
###############################
### Fit des Anfangsbereichs ###
###############################

print(current1_2[12:41]) # Wie man sieht sind hier negative Werte am Anfang drin.
# Das ist ein Problem. Ich nehme sie raus mit dem Argument, dass es am Anfang nur "Rauschen ist" (?)
params, covariance_matrix = optimize.curve_fit(linfit, 1/temp1_2[23:38], np.log(current1_2[23:38]))
a, b = correlated_values(params, covariance_matrix)
print('Fit für die langsame Heizrate:')
print('a=', a)
print('b=', b)
print('W=', -a*k_B)


xlin = np.linspace(0.0035,0.0045)
plt.plot(1/temp1_2[16:41], np.log(current1_2[16:41]), 'bx', mew=0.5, label = 'Messwerte')
plt.plot(1/temp1_2[23:38], np.log(current1_2[23:38]), 'rx', mew=0.5, label = 'Für die Ausgleichsrechnung verwendete Messwerte')
plt.plot(xlin, linfit(xlin,*params), label= 'Ausgleichsrechnung', linewidth=0.5)
plt.xlim(0.00387,0.0044)
plt.grid()
plt.xlabel(r'$(1/T)/$(1/K)')
plt.ylabel(r'$Log I$')
plt.legend(loc='best')
plt.savefig('build/fit1.pdf')
plt.clf()

print(current2[12:30])#Hier müssen auch noch ein paar Werte rausgenommen werden!

params, covariance_matrix = optimize.curve_fit(linfit, 1/temp2[17:28], np.log(current2[17:28]))
a, b = correlated_values(params, covariance_matrix)
print('Fit für die schnelle Heizrate:')
print('a=', a)
print('b=', b)
print('W=', -a*k_B)

plt.plot(1/temp2[15:30], np.log(current2[15:30]), 'bx', mew=0.5, label = 'Messwerte')
plt.plot(1/temp2[17:28], np.log(current2[17:28]), 'rx', mew=0.5, label = 'Für die Ausgleichsrechnung verwendete Messwerte')
xlin = np.linspace(0.0035,0.0043)
plt.plot(xlin, linfit(xlin,*params), label= 'Ausgleichsrechnung', linewidth=0.5)
plt.xlim(0.0038,0.00422)
#plt.xlim(-75,55)
plt.grid()
plt.xlabel(r'$(1/T)/$(1/K)')
plt.ylabel(r'$Log(I)$')
plt.legend(loc='best')
plt.savefig('build/fit2.pdf')
plt.clf()

#####################
### Second method ###
#####################

#lnint1_2 = lnint(current1_2[16:41], temp1_2[16:41], noms(heatRate[0]))
#print(lnint1_2)
#ln1_2 = np.log(lnint1_2)
#params, covariance_matrix = optimize.curve_fit(hyp, ln1_2, 1/temp1_2[16:41])
#a, b = correlated_values(params, covariance_matrix)
#print('Fit für die langsame Heizrate:')
#print('a=', a)
#print('b=', b)

T1 = temp1_2[16:41]
yT1 = current1_2[16:41]
hr1 = 1.1

m2y1 = lnint(T1,yT1,hr1)
print(T1)
print(yT1)
m2x1 = T1[m2y1>0]
m2y1 = m2y1[m2y1>0]
m2y1 = np.log(m2y1)
m2x1 = 1/m2x1
plt.plot(m2x1, m2y1, 'rx')
plt.savefig('build/test.pdf')
plt.clf()

#plt.plot(1/temp1_2[16:41], ln1_2,'bx', alpha=0.75, label = 'Messwerte')
#plt.ylim(-10,45)
#plt.xlim(-75,55)
plt.plot(1)
plt.grid()
plt.xlabel(r'$(1/T)/$(1/K)')
plt.ylabel(r'$Log von dem großen Scheißintegral$')
plt.legend(loc='best')
plt.savefig('build/testFuerIntfit.pdf')
plt.clf()


##############################################
### Extract tau_0 from the maximum of I(T) ###
##############################################

# Beachte, dass wir für diese Berechnung W brauchen. Es wäre am konsistentesten,
# W aus der zweiten Methode zu verwenden.
# Bitte so coden, dass das W aus der zweiten Methode verwendet wird! Danke!

# Hier als Test: W ist 1eV
W = 1.602e-19

#for i in (0,1):
#    ### Die untenstehende Taktik ist statt argmax nötig, weil argmax
#    ### nur den ersten index gibt, wo das erste mal der maximale wert kommt
#    indices_Imax = np.argwhere(currents[i] == np.amax(currents[i])) # amax means maximum of array
#    T_max = temps[i][indices_Imax]
#    print('T_max beträgt', mittel(T_max, deltadof=0)) # Welche ddof sollen wir nehmen?
#    tau_0 = kboltzmann*T_max**2/(W*b[i]) * unp.exp(-W/(kboltzmann*T_max))
#    print(tau_0)
#
