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
print(I_0)

####################
######Spektrum######
####################

#Muss man anscheinend gegen einen linspace auftragen weil es einfach keine x-Werte gibt


#####################
######Wuerfel 1######
#####################



#####################
######Wuerfel 2######
#####################



#####################
######Wuerfel 3######
#####################



#####################
######Wuerfel 4######
#####################
s2 = 2
A1 = np.matrix([[1,1,1,0,0,0,0,0,0],[0,0,0,1,1,1,0,0,0],[0,0,0,0,0,0,1,1,1]])
A2 = np.matrix([[1,0,0,1,0,0,1,0,0],[0,1,0,0,1,0,0,1,0],[0,0,1,0,0,1,0,0,1]])
A3 = np.matrix([[0,s2,0,s2,0,0,0,0,0],[0,0,s2,0,s2,0,s2,0,0],[0,0,0,0,0,s2,0,s2,0]])
A4 = np.matrix([[0,0,0,s2,0,0,0,s2,0],[s2,0,0,0,s2,0,0,0,s2],[0,s2,0,0,0,s2,0,0,0]])
A = np.vstack((A1,A2,A3,A4))
print(np.linalg.inv(A.T@A))
