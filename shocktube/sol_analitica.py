import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as  slv

"""
Referencia para solucion
http://www.astro.uu.se/~hoefner/astro/teach/ch10.pdf
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

U = 0.0

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



#Llenado de variables en las regiones

x = np.linspace(0, 1, n_points)
rho = np.zeros(n_points)
P = np.zeros(n_points)
u = np.zeros(n_points)
e = np.zeros(n_points)


for i in range(0, n_points):
    #Region izquierda
    if(x[i] < x1):
        rho[i] = rho_l
        P[i] = Pl
        u[i] = U
        e[i] = P[i]/(gamma - 1) + 0.5*rho[i]*u[i]**2
        
    #Region expansion
    elif(x[i] >= x1 and x[i] < x2 ):
        u[i] = (2/(gamma + 1))*(al + (x[i] - x0)/t)
        a = al - (gamma - 1)*u[i]*0.5
        P[i] = Pl*(a/al)**(1/g3)
        rho[i] = gamma*P[i]/a**2        
        e[i] = P[i]/(gamma - 1) + 0.5*rho[i]*u[i]**2
    
    #Region 2
    elif(x[i] >= x2 and x[i] < x3):       
        P[i] = Pl*(a2/al)**(1/g3)
        rho[i] = rho_l*(P[i]/Pl)**(1/gamma)
        u[i] = U2
        e[i] = P[i]/(gamma - 1) + 0.5*rho[i]*u[i]**2

    #Region 1
    elif(x[i] >= x3 and x[i] < x4):
        
        rho[i] = rho_r/((2/(gamma+1))*(1/Ms**2) + g2)
        P[i] = Pl*(a2/al)**(1/g3)
        u[i] = U2
        e[i] = P[i]/(gamma - 1) + 0.5*rho[i]*u[i]**2
    #Region derecha
    else:
        rho[i] = rho_r
        P[i] = Pr
        u[i] = U
        e[i] = P[i]/(gamma - 1) + 0.5*rho[i]*u[i]**2


#Importacion de datos en archivo
datos = np.loadtxt('UpwindGodunov_step_1.dat')
rho_num = datos[:,1]
u_num = datos[:,2]
e_num = datos[:,3]
P_num = datos[:,4]

#Creacion de graficas
plt.figure(1)
plt.plot(x, rho)
plt.plot(x, rho_num)
plt.title("Grafica de Densidad")
plt.xlabel(r'$x$')
plt.ylabel(r'$\rho$')
plt.savefig("Densidad.png")

plt.figure(2)
plt.plot(x, u)
plt.plot(x, u_num)
plt.title("Grafica de Velocidad")
plt.xlabel(r'$x$')
plt.ylabel(r'$Velocidad$')
plt.savefig("Velocidad.png")

plt.figure(3)
plt.plot(x, e)
plt.plot(x, e_num)
plt.title("Grafica de Energia")
plt.xlabel(r'$x$')
plt.ylabel(r'$Energia$')
plt.savefig("Energia.png")

plt.figure(4)
plt.plot(x, P)
plt.plot(x, P_num)
plt.title("Grafica de Presion")
plt.xlabel(r'$x$')
plt.ylabel(r'$P$')
plt.savefig("Presion.png")



        
