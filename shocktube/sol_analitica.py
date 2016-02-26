import numpy as np
import matplotlib.pyplot 
import sys
import scipy.optimize as  slv

"""
Referencias para solucion

http://www.astro.uu.se/~hoefner/astro/teach/ch10.pdf
http://www.phys.lsu.edu/~tohline/PHYS7412/sod.html
"""

n_points = 5000 #Numero de puntos en discretizacion

x0 = 0.5
t = 0.2

#Constantes
gamma = 1.4
g2 = (gamma - 1)/(gamma + 1)
g3 = (gamma-1)/(2*gamma)

#Condiciones iniciales

rho_l  = 1.0
Pl = 1.0

rho_r = 0.125
Pr = 0.1

u = 0.0

al   = np.sqrt(gamma*Pl/rho_l) #Velocidad del sonido en parte izquierda

#Ecuacion para encontrar las raices

def eq_implicita(Ms):
    """parte1 = ((1- g2)**2)/(rho_r*(P + g2*Pr))
    parte2 = (2*np.sqrt(gamma))/((gamma -1))*(1 - P**g3)
    return (P - Pr)*np.sqrt(parte1) - parte2"""
    parte1 = Ms-(1/Ms)
    parte2 = (Pr/Pl)*((Ms*Ms/g3)- g2)
    parte3 = (al/g2)*(1-parte2**g3)
    return parte3 - parte1
   
#Solucion para Ppost usando solver de scipy optimize
Ms = slv.fsolve(eq_implicita, 0.31)

#Variables necesarias para calculo de velocidades
U2 = (Ms-(1/Ms))/(g3*gamma)
a2 = al- U2*g3*gamma
v_2 = U2 - a2 


#Calculo de posiciones de onda
x1 = x0 - al*t
x2 = x0 + v_2*t
x3 = x0 + U2*t
x4 = x0 + Ms*t 

print x1, x2, x3 , x4

#Llenado de variables en las regiones

rho = np.zeros(n_points)
P = np.zeros(n_points)
u = np.zeros(n_points)
e = np.zeros(n_points)

for i in range(0, n_points):
    if(i < x1):
        rho[i] = rho_l
        P[i] = Pl
        
