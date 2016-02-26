import numpy as np
import matplotlib.pyplot 
import sys
import scipy.optimize as  slv

"""
Referencias para solucion

http://www.astro.uu.se/~hoefner/astro/teach/ch10.pdf
http://www.phys.lsu.edu/~tohline/PHYS7412/sod.html
"""



#Constantes
gamma = 1.4
g2 = (gamma - 1)/(gamma + 1)
g3 = (gamma-1)/(2*gamma)

#Condiciones iniciales

rho_l  = 1.0
Pl = 1.0

rho_r = 0.125
Pr = 0.1

u = 1.0

#Ecuacion para encontrar las raices

def eq_implicita(P):
    parte1 = ((1- g2)**2)/(rho_r*(P + g2*Pr))
    parte2 = (2*np.sqrt(gamma))/((gamma -1)*(1 - P**g3))
    return (P - Pr)*np.sqrt(parte1) - parte2
   

Ppost = slv.fsolve(eq_implicita, 0.31)

print Ppost
