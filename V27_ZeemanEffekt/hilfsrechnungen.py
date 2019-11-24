import numpy as np
import scipy.optimize
import uncertainties as unc
import uncertainties.unumpy as unp
from scipy import optimize
import matplotlib.pyplot as plt
from uncertainties import ufloat


def delta_lam_D(lam, n, d=4e-3):
    return lam*lam/(2*d)/np.sqrt(n*n-1)


def A(lam, n, L=120e-3):
    return L/lam*(n*n-1)


def n(lam):
    if lam==643.8e-9: return 1.4567
    elif lam==480e-9: return 1.4635
    else: return 0


def B(delta_lam, n, lam_0, h=6.62607e-34, c=299.792458e6, mu_B=9.274e-24):
    return h*c/(delta_lam*n*mu_B*lam_0)*(np.sqrt(lam_0*lam_0-delta_lam*delta_lam)-lam_0)


## Hilfsrechner für Vorbereitungsaufgaben 6 und 7

print('-----------------------------------')
print('Werte für Dispersionsgebiete und Auflösungsvermögen der Lummer-Gehrke-Platte (Vorbereitungsaufgabe 6):')
print('-----------------------------------')
print('Für lambda=643.8nm:')
lam = 643.8e-9
delta_lam_rot = delta_lam_D(lam=lam,n=n(lam))
print('delta_lam_D=', delta_lam_rot)
print('A=', A(lam, n=n(lam)))
print('-----------------------------------')
print('Für lambda=480nm:')
lam = 480e-9
delta_lam_blau = delta_lam_D(lam=lam,n=n(lam))
print('delta_lam_D=', delta_lam_blau)
print('A=', A(lam, n=n(lam)))

print('-----------------------------------')
print('Maximale B-Feldstärken (Vorbereitungsaufgabe 7):')
print('-----------------------------------')
print('Für die rote Linie (normaler Zeeman-Effekt):')
print('B=', B(delta_lam_rot, 1, 643.8e-9), ' (Keine Ahnung warum hier ein Minus steht)')
print('Für die blaue Linie, sigma Polarisation:')
print('B=', B(delta_lam_blau, 2, 480e-9))
print('Für die blaue Linie, pi Polarisation:')
print('B=', B(delta_lam_blau, 1/2, 480e-9))
