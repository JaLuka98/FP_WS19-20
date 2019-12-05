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


#au(2 µm):
N_au1=905
t_au1=480
#au(4µm):
N_au2=836
t_au2=480

#alu(2µm):
N_pt=901
t_pt=480

#bis(1µm):
N_bis= 985
t_bis=480

########################
###Fits für die Dicke###
########################
print('##############Foliendicke################')

p_f, U_f = np.genfromtxt('data/mit_folie.txt', unpack=True)
p, U = np.genfromtxt('data/ohne_folie.txt', unpack=True)


hr = ['$p_\text{Folie}$/mbar', '$U_{\text{Folie}}$/V', '$p$/mbar', '$U$/V']
m = np.zeros((18, 4))
m[:,0] = np.append(p_f, np.zeros(8))
m[:,1] = np.append(U_f, np.zeros(8))
m[:,2] = p
m[:,3] = U
t=matrix2latex(m, headerRow=hr, format='%.2f')
print(t)

params1, covariance_matrix1 = optimize.curve_fit(linfit, p_f, U_f)
a1, b1 = correlated_values(params1, covariance_matrix1)
print('Fit für mit Folie:')
b_folie = b1
print('a=', a1)
print('b=', b1)

params2, covariance_matrix2 = optimize.curve_fit(linfit, p[8:], U[8:])
a2, b2 = correlated_values(params2, covariance_matrix2)
print('Fit für mit Folie:')
b = b2
print('a=', a2)
print('b=', b2)

linspace=np.linspace(0, 300, 1000)
plt.plot(linspace, linfit(linspace, *params1), 'r-', label='Fit mit Folie', linewidth=0.5)
plt.plot(linspace, linfit(linspace, *params2), 'b-', label='Fit ohne Folie', linewidth=0.5)
plt.plot(p_f, U_f, 'rx', mew=0.5, label='Messwerte mit Folie')
plt.plot(p, U, 'bx', mew=0.5, label='Messwerte ohne Folie')
plt.plot(p[:8], U[:8], 'gx', mew=0.5)
plt.xlabel(r'$p/$mBar')
plt.ylabel(r'$U$/V')
plt.xlim(0, 300)
plt.tight_layout()
plt.legend()
plt.grid()
plt.savefig('build/dicke.pdf')
plt.clf()

#############################################
###Bestimung von Dicke siehe Nebenrechnung###
#############################################

E_alpha = 5.49
DeltaE = E_alpha * (1 - b_folie/b2)
print('DeltaE beträgt: ', DeltaE)

##########################
###Raumwinkel bestimmen###
##########################
print('##############Raumwinkel################')
a = 2
c = 41
d = 45
b = a*d/c #Breite
l = 10
h = l*d/c #Höhe

omega = 4 * np.arctan((b*h)/(2*d*np.sqrt(4*d**2 + b**2 + h**2)))

print('Raumwinkel für die Folie:', omega)

c = 97
d = 101
omega = 4 * np.arctan((b*h)/(2*d*np.sqrt(4*d**2 + b**2 + h**2)))

print('Raumwinkel für die Probe:', omega)

#########################
###Aktivität bestimmen###
#########################
print('##############Aktivität################')
A_1994=330 #kBq
t=25+1/6 #Jahre
t_halb=432 #Jahre

lamda=np.log(2)/t_halb

A_heute= A_1994*np.exp(-lamda*t)

print('Zerfallskonstante:', lamda)
print('Aktivität heute:', A_heute)



#####################################
###Plot für den differentiellen WQ###
#####################################
print('##############Differentieller WQ################')
I0=2359/150
sigmaI0=np.sqrt(2359)/150

A_heute=4*np.pi*I0/omega
sigmaA=4*np.pi*sigmaI0/omega
print('Experimentall: A_heute=', A_heute, '±', sigmaA)

print('I_0=:', I0, '+/-', sigmaI0)

counts, t, theta = np.genfromtxt('data/winkel.txt', unpack=True)

sigmaCounts=np.sqrt(counts)
I=counts/t #Aktivität pro Raumwinkel
sigmaI=sigmaCounts/t

N_A = 6.022*1e23
M=196.97 #molmasse
dicke = 2*1e-6       #Foliendicke
rho= 19.32

dichte =  N_A*rho/M  #Atomdichte

print('Atomdichte:', dichte, 'pro cm^-3')
print('Atomdichte:',  dichte*1e6,'pro m^-3')

dichte*=1e6 #Teilchenanzahl pro m^3


differentiellerWQ = I/(I0*dichte*dicke*omega)
sigmaWQ=np.sqrt((sigmaI/(I0*dichte*dicke*omega))**2 + ((I*sigmaI0/(I0**2*dichte*dicke*omega)))**2)
#Die Fehler hier eskalieren irgendwie ein bisschen, aber das I_0 ist ja auch fehlerbehaftet...


hr = ['$N$', '','$t$/s', '$I/\mathrm{s^-1}$','' ,'$\theta$/°', '\frac{d\sigma}{d \Omega}',]
m = np.zeros((33, 8))
m[:,0] = counts
m[:,1] = sigmaCounts
m[:,2] = t
m[:,3] = I
m[:,4] = sigmaI
m[:,5] = theta
m[:,6] = differentiellerWQ *1e21
m[:,7] = sigmaWQ*1e21
t=matrix2latex(m, headerRow=hr, format='%.2f')
print(t)

########Wert für E hier noch anpassen und Quelle direkt mit raussuchen!!!

def rutherford(theta):
   return (1/(4*np.pi*8.854*1e-12))**2 *((2*79*(1.602*1e-19)**2)/(4*5.638*1e6*1.602*1e-19))**2 *(1/(np.sin(np.deg2rad(theta/2))))**4 #hier noch E_0 anpassen???

def rutherfordverschoben(theta):
    return (1/(4*np.pi*8.854*1e-12))**2 *((2*79*(1.602*1e-19)**2)/(4*5.638*1e6*1.602*1e-19))**2 *(1/(np.sin(np.deg2rad((theta-2.5)/2))))**4

linspace=np.linspace(-7,-0.01,1000)
linspace2=np.linspace(0.01, 11, 1000)
plt.errorbar(theta, differentiellerWQ, yerr=sigmaWQ, fmt = 'x', elinewidth=0.5, markeredgewidth=0.5, label='Experimenteller Wirkungsquerschnitt')
plt.plot(linspace, rutherford(linspace), label='Theoretischer Wirkungsquerschnitt', color='red', linewidth=0.5)
plt.plot(linspace, rutherfordverschoben(linspace), label='Theoretischer Wirkungsquerschnitt (verschoben)', color='green', linewidth=0.5)
plt.plot(linspace2, rutherford(linspace2), color='red',linewidth=0.5)
plt.plot(linspace2, rutherfordverschoben(linspace2), color='green',linewidth=0.5)
plt.xlabel(r'$\theta/$°')
plt.ylabel(r'$\frac{d\sigma}{d\Omega}/\mathrm{m^2}$')
plt.xlim(-7, 11)
plt.ylim(-0.05*1e-20, 0.7*1e-20)
plt.tight_layout()
plt.legend()
plt.grid()
plt.savefig('build/wq.pdf')
plt.clf()


#######################
###Mehrfachstreuuung###
#######################
print('##############Mehrfachstreuung################')

dicke, theta, t, counts=np.genfromtxt('data/mehrfach.txt', unpack=True)
dicke*=1e-6 #Dicke in m
sigmaCounts=np.sqrt(counts)
I=counts/t #Aktivität
sigmaI=sigmaCounts/t

print('Nominalwerte für Mehrfachstreuung:', I)
print('Unsicherheiten für Mehrfachstreuung:', sigmaI)

differentiellerWQ = I/(I0*dichte*dicke*omega)
sigmaWQ=np.sqrt((sigmaI/(I0*dichte*dicke*omega))**2 + ((I*sigmaI0/(I0**2*dichte*dicke*omega)))**2)
print('WQ:', differentiellerWQ)
print('sigmaWQ:', sigmaWQ)


####################
###Z-Abhängigkeit###
####################
print('##############Z-Abhängigkeit################')
#au(2 µm):
N_au1=905
t_au1=480
#au(4µm):
N_au2=836
t_au2=480
#platin(2µm):
N_pt=901
t_pt=480
#bis(1µm):
N_bis= 985
t_bis=480

N=np.array([905, 836, 901, 985])
I=N/480
sigmaI=np.sqrt(N)/480

####Quellen für die Werte suchen und Werte ggf. anpassen!

rho = np.array([19.32, 19.32 ,21.45, 9.8]) * 1e6 #für in m^3
M=np.array([196.97, 196.97, 195.09, 208.98])   #molare Massen
dicke=np.array([2,4,2,1])*1e-6 #Dicke in µm
Z=np.array([79,79,78,83])
N_A = 6.022*1e23

N_abs=N_A *rho/M

komischegroesse=I/(N_abs*dicke)
sigmakomischegroesse=sigmaI/(N_abs*dicke)

print('Fehler komische Groesse:', sigmakomischegroesse)

hr = ['$N$','' ,'A/\mathrm{s^{-1}$', '', '$\rho/\mathrm{\frac{kg}{m^3}$', '$M/\mathrm{\frac{g}{Mol}}$', '$x$/µm', '$Z$', '$\frac{A}{Nx}/10^-23\mathrm{\frac{m^2}{s}}$']
m = np.zeros((4, 9))
m[:,0] = N
m[:,1] = np.sqrt(N)
m[:,2] = I
m[:,3] = sigmaI
m[:,4] = rho
m[:,5] = M
m[:,6] = dicke*1e6
m[:,7] = Z
m[:,8] = komischegroesse *1e23
t=matrix2latex(m, headerRow=hr, format='%.2f')
print(t)

#Ist Intensität hier Aktivität?

#Jetzt die komischegroesse gegen die kernladungszahl auftragen

plt.errorbar(Z, komischegroesse, yerr=sigmakomischegroesse, fmt = 'x', elinewidth=0.5, markeredgewidth=0.5)
plt.ylabel(r'$\frac{I}{N \Delta x}/\frac{\mathrm{m^2}}{\mathrm{s}}$')
plt.xlabel(r'$Z$')
#plt.xlim(-10, 10)
plt.tight_layout()
#plt.legend()
plt.grid()
plt.savefig('build/z.pdf')
plt.clf()

#Jetzt irgendwie die Theoriewerte dazu berechnen:
def rutherford(theta, Z):
   return (1/(4*np.pi*8.854*1e-12))**2 *((2*Z*(1.602*1e-19)**2)/(4*5.5*1e6*1.602*1e-19))**2 *(1/(np.sin(np.deg2rad(theta/2))))**4

theogroesse=rutherford(6, Z)*(I0*omega)
sigmaTheogroesse=rutherford(6, Z)*sigmaI0*omega

print('theogroesse:', theogroesse)
print('sigma theogroesse:', sigmaTheogroesse)
