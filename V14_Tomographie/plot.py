import numpy as np
import scipy.optimize
import uncertainties as unc
import uncertainties.unumpy as unp
from scipy import optimize
import matplotlib.pyplot as plt
from uncertainties import ufloat

t=180
t_0=300
t_1=120

I_0= unp.uarray([55224], [262])/300
print('I_0=', I_0)

####################
######Spektrum######
####################
print('______________________SPEKTRUM_____________________')
hoehe=np.genfromtxt('data/spektrum.txt', unpack=True)

x=np.arange(1,513)
plt.plot(x, hoehe)
plt.xlabel(r'channel')
plt.ylabel(r'Amplitude')
plt.tight_layout()
plt.grid()
plt.savefig('build/spektrum.pdf')
plt.clf()


#####################
######Wuerfel 1######
#####################
print('______________________WUERFEL_1_____________________')
N0, sigma0=np.genfromtxt('data/wuerfel1.txt', unpack=True)
N0/=120
sigma0/=120

I_0= unp.uarray(N0, sigma0)

print('I_0=', I_0)
print(N0[1])

#####################
######Wuerfel 2######
#####################
print('______________________WUERFEL_2_____________________')
N2, sigma2=np.genfromtxt('data/wuerfel2.txt', unpack=True)
N2/=180
sigma2/=180

I_2=unp.uarray(N2, sigma2)

l=np.array([3,3,3*np.sqrt(2), 2*np.sqrt(2)])
N0_i=np.insert(N0, 0, N0[0])
sigma0_i=np.insert(sigma0, 0, sigma0[0])

mu=1/l*np.log(N0_i/N2)
sigmamu=np.sqrt((sigma2/(l*N2))**2 + (sigma0_i/(l*N0_i)**2))

mittelwert=1/4*sum(mu)
sigmamittelwert=np.sqrt(1/12*sum(mittelwert-mu)**2)

print('mu=', mu, '+/-', sigmamu)

#####################
######Wuerfel 3######
#####################
print('______________________WUERFEL_3_____________________')
N3, sigma3=np.genfromtxt('data/wuerfel3.txt', unpack=True)
N3/=180
sigma3/=180

mu=1/l*np.log(N0_i/N3)
sigmamu=np.sqrt((sigma2/(l*N2))**2 + (sigma0_i/(l*N0_i))**2)

mittelwert=1/4*sum(mu)
sigmamittelwert=np.sqrt(1/12*sum(mittelwert-mu)**2)

print('mu=', mu, '+/-', sigmamu)

#####################
######Wuerfel 4######
#####################
print('______________________WUERFEL_4_____________________')
s2 = np.sqrt(2)
A1 = np.matrix([[1,1,1,0,0,0,0,0,0],[0,0,0,1,1,1,0,0,0],[0,0,0,0,0,0,1,1,1]])
A2 = np.matrix([[1,0,0,1,0,0,1,0,0],[0,1,0,0,1,0,0,1,0],[0,0,1,0,0,1,0,0,1]])
A3 = np.matrix([[0,s2,0,s2,0,0,0,0,0],[0,0,s2,0,s2,0,s2,0,0],[0,0,0,0,0,s2,0,s2,0]])
A4 = np.matrix([[0,0,0,s2,0,0,0,s2,0],[s2,0,0,0,s2,0,0,0,s2],[0,s2,0,0,0,s2,0,0,0]])
A = np.vstack((A1,A2,A3,A4))
print('A=',A)

B=np.linalg.inv(A.T@A)
print('C=1(d^2*', B)

C=1/9*B

N4, sigma4=np.genfromtxt('data/wuerfel4.txt', unpack=True)
N4/=180
sigma4/=180

N0_i=np.array([N0[0], N0[0], N0[0], N0[0], N0[0], N0[0], N0[1],N0[2], N0[1], N0[1], N0[2], N0[1]])
sigma0_i=np.array([sigma0[0],sigma0[0],sigma0[0],sigma0[0],sigma0[0],sigma0[0],sigma0[1],sigma0[2],sigma0[1],sigma0[1],sigma0[2],sigma0[1]])

mu= 3*C@A.T@np.log(N0_i.T/N4.T)

print('mu=', mu)
