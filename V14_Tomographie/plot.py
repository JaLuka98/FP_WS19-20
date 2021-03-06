import numpy as np
import scipy.optimize
import uncertainties as unc
import uncertainties.unumpy as unp
from scipy import optimize
import matplotlib.pyplot as plt
from uncertainties import ufloat
from matrix2latex import matrix2latex
from uncertainties import ufloat
from uncertainties import correlated_values

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
plt.xlabel(r'Channel')
plt.ylabel(r'Anzahl der Ereignisse')
plt.axis((0,150,-250,7000))
plt.text(37,1500,'1')
plt.text(55,550,'2')
plt.text(76,1200,'3')
plt.text(108,6500,'4')
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

N0_i_with_errs = unp.uarray(N0_i, sigma0_i)
N2_with_errs = unp.uarray(N2, sigma2)

mu2_with_errs = 1/l*unp.log(N0_i_with_errs/N2_with_errs)
#mu2=1/l*np.log(N0_i/N2)
#sigmamu2=np.sqrt((sigma2/(l*N2))**2 + (sigma0_i/(l*N0_i)**2))
#
#mittelwert=1/4*sum(mu2)
#sigmamittelwert=np.sqrt(1/12*sum(mittelwert-mu2)**2)
#sigmamittelwert=1/4*sum(sigma2)
#
##IRGENDWAS STIMMT HIER MIT DEM MITTELWERT NOCH NICHT!!!
#
#print('mu=', mu2, '+/-', sigmamu2)
#print('Mittelwert=', mittelwert, '+/-', sigmamittelwert)
print('mu=', mu2_with_errs)
print('Mittelwert', mu2_with_errs.mean())
mu2 = unp.nominal_values(mu2_with_errs)
sigmamu2 = unp.std_devs(mu2_with_errs)

#####################
######Wuerfel 3######
#####################
print('______________________WUERFEL_3_____________________')
N3, sigma3=np.genfromtxt('data/wuerfel3.txt', unpack=True)
N3/=180
sigma3/=180

#mu3=1/l*np.log(N0_i/N3)
#sigmamu3=np.sqrt((sigma2/(l*N2))**2 + (sigma0_i/(l*N0_i))**2)

#mittelwert=1/4*sum(mu3)
#sigmamittelwert=np.sqrt(1/12*sum(mittelwert-mu3)**2)
#sigmamittelwert=1/4*sum(sigma3)

N0_i_with_errs = unp.uarray(N0_i, sigma0_i)
N3_with_errs = unp.uarray(N3, sigma3)

mu3_with_errs = 1/l*unp.log(N0_i_with_errs/N3_with_errs)

#print('mu=', mu3, '+/-', sigmamu3)
#print('Mittelwert=', mittelwert, '+/-', sigmamittelwert)

print('mu=', mu3_with_errs)
print('Mittelwert', mu3_with_errs.mean())
mu3 = unp.nominal_values(mu3_with_errs)
sigmamu3 = unp.std_devs(mu3_with_errs)

hr = ['$N_2$', '','$I_2$/(1/s)', '', '$N_3$', '', '$I_3$/(1/s)', '']
m = np.zeros((4, 8))
m[:,0] = N2*180
m[:,1] = sigma2*180
m[:,2] = N2
m[:,3] = sigma2
m[:,4] = N3*180
m[:,5] = sigma3*180
m[:,6] = N3
m[:,7] = sigma3
t=matrix2latex(m, headerRow=hr, format='%.2f')
print(t)

hr = ['$mu_2$', '', '$mu_3$', '']
m = np.zeros((4, 4))
m[:,0] = mu2
m[:,1] = sigmamu2
m[:,2] = mu3
m[:,3] = sigmamu3
t=matrix2latex(m, headerRow=hr, format='%.3f')
print(t)

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

N4, sigma4=np.genfromtxt('data/wuerfel4.txt', unpack=True)
N4/=180
sigma4/=180

hr = ['$N_4$', '', '$I_4$', '']
m = np.zeros((12, 4))
m[:,0] = N4*180
m[:,1] = sigma4*180
m[:,2] = N4
m[:,3] = sigma4
t=matrix2latex(m, headerRow=hr, format='%.2f')
print(t)


N0_i=np.array([N0[0], N0[0], N0[0], N0[0], N0[0], N0[0], N0[1],N0[2], N0[1], N0[1], N0[2], N0[1]])
sigma0_i=np.array([sigma0[0],sigma0[0],sigma0[0],sigma0[0],sigma0[0],sigma0[0],sigma0[1],sigma0[2],sigma0[1],sigma0[1],sigma0[2],sigma0[1]])

#N0_iges=np.array([N0ges[0], N0ges[0], N0ges[0], N0ges[0], N0ges[0], N0ges[0], N0ges[1],N0ges[2], N0ges[1], N0ges[1], N0ges[2], N0ges[1]])

#print((N0_i/N4).T)
N0_i_with_errs = unp.uarray(N0_i, sigma0_i)
N4_with_errs = unp.uarray(N4, sigma4)

mu= np.linalg.inv(A.T@A)@A.T@np.log((N0_i/N4)).T
#####FEHLER VON MU NOCH UNBEDINGT AUSRECHNEN!!!
V_I = unp.std_devs(unp.log((N0_i_with_errs/N4_with_errs)))
V_I = np.diag(V_I**2) # Covariance matrix has variance on main diagonal
V_mu = np.linalg.inv(A.T@A)@A.T@V_I@A@(np.linalg.inv(A.T@A))
sigma_mu = np.sqrt(np.diag(V_mu)) # we want the std devs
print('mu', mu)
print('sigma_mu', sigma_mu)
