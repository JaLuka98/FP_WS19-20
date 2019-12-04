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
print('a=', a1)
print('b=', b1)

params2, covariance_matrix2 = optimize.curve_fit(linfit, p[8:], U[8:])
a2, b2 = correlated_values(params2, covariance_matrix2)
print('Fit für mit Folie:')
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

##########################
###Raumwinkel bestimmen###
##########################
a = 2
c = 97
#####Hier noch mal überprüfen welche Werte jetzt wofür genau richtig sind!!!#######
d = 101
b = a*d/c #Breite
l = 10
h = l*d/c #Höhe

omega = 4 * np.arctan((b*h)/(2*d*np.sqrt(4*d**2 + b**2 + h**2)))

#########################
###Aktivität bestimmen###
#########################
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
#N0=

#################################################
I0= 10 #DAS IST EIN DUMMY WERT!!!
#################################################


counts, t, theta = np.genfromtxt('data/winkel.txt', unpack=True)

sigmaCounts=np.sqrt(counts)
I=counts/t #Aktivität pro Raumwinkel
sigmaI=sigmaCounts/t

dicke = 2*1e-6       #Foliendicke
dichte =5.907*1e28    #Atomdichte

differentiellerWQ = I/(I0*dichte*dicke*omega)
sigmaWQ=np.sqrt((sigmaI/(I0*dichte*dicke*omega))**2 + ((0.26*I/(I0**2*dichte*dicke*omega)))**2) #Woher kommt die 0.26?

hr = ['$N$', '','$t$/s', '$I/\mathrm{s^-1}$','' ,'$\theta$/°', '\frac{d\sigma}{d \Omega}',]
m = np.zeros((33, 8))
m[:,0] = counts
m[:,1] = sigmaCounts
m[:,2] = t
m[:,3] = I
m[:,4] = sigmaI
m[:,5] = theta
m[:,6] = differentiellerWQ
m[:,7] = sigmaWQ
t=matrix2latex(m, headerRow=hr, format='%.2f')
print(t)

########Wert für E hier noch anpassen und Quelle direkt mit raussuchen!!!

def rutherford(theta):
   return (1/(4*np.pi*8.854*1e-12)**2)*(((2*79*(1.602*1e-19)**2)/(4*5.638*1e6))**2)*1/(np.sin(np.deg2rad(theta/2)))**4 #hier noch E_0 anpassen???

linspace=np.linspace(-7,11,1000)
plt.errorbar(theta, differentiellerWQ, yerr=sigmaWQ, fmt = 'x', label='Aus Messwerten berechneter Differentieller Wirkungsquerschnitt')
plt.plot(linspace, rutherford(linspace), label='Theoriekurve', linewidth=0.5) #Die funktioniert irgendwie noch nicht...
plt.xlabel(r'$\theta/$°')
plt.ylabel(r'$\frac{d\sigma}{d\Omega}/\mathrm{m^2}$')
plt.xlim(-7, 11)
plt.tight_layout()
plt.legend()
plt.grid()
plt.savefig('build/wq.pdf')
plt.clf()


#######################
###Mehrfachstreuuung###
#######################

dicke, theta, t, counts=np.genfromtxt('data/mehrfach.txt', unpack=True)
dicke*=1e-6 #Dicke in m
sigmaCounts=np.sqrt(counts)
I=counts/t #Aktivität
sigmaI=sigmaCounts/t

differentiellerWQ = I/(I0*dichte*dicke*omega)
sigmaWQ=np.sqrt((sigmaI/I0*dichte*dicke*omega)**2 + ((0.26* I/I0**2*dichte*dicke*omega))**2) #Woher kommt die 0.26?

print('WQ:', differentiellerWQ)
print('sigmaWQ:', sigmaWQ)


####################
###Z-Abhängigkeit###
####################
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

rho = np.array([19.30, 19.30 ,21.45, 9.78]) * 1e6 #für in m^3
M=np.array([196.97, 196.97, 195.08, 208.98])   #molare Massen
dicke=np.array([2,4,2,1])*1e-6 #Dicke in µm
Z=np.array([79,79,78,83])
N_A = 6.022*1e23

N_abs=N_A *rho/M

komischegroesse=I/(N_abs*dicke)

hr = ['$N$','' ,'A/\mathrm{s^{-1}$', '', '$\rho/\mathrm{\frac{kg}{m^3}$', '$M/\mathrm{\frac{g}{Mol}}$', '$x$/µm', '$Z$', '$\frac{A}{Nx}/\mathrm{\frac{m^2}{s}}$']
m = np.zeros((4, 9))
m[:,0] = N
m[:,1] = np.sqrt(N)
m[:,2] = I
m[:,3] = sigmaI
m[:,4] = rho
m[:,5] = M
m[:,6] = dicke
m[:,7] = Z
m[:,8] = komischegroesse
t=matrix2latex(m, headerRow=hr, format='%.2f')
print(t)

#Ist Intensität hier Aktivität?

#Jetzt die komischegroesse gegen die kernladungszahl auftragen

plt.plot(Z, komischegroesse, 'rx', mew=0.5)
plt.ylabel(r'$I/(N*x)\frac{\mathrm{m^2}}{\mathrm{s}}$')
plt.xlabel(r'$Z$')
#plt.xlim(-10, 10)
plt.tight_layout()
plt.legend()
plt.grid()
plt.savefig('build/z.pdf')
plt.clf()
