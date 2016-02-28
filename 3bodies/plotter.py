import numpy as np 
import matplotlib.pyplot as plt
import os as os

data_directory = './Results/'
plots_directory = './Figures/'

if not os.path.exists(plots_directory):
    os.makedirs(plots_directory)

# plots Figure 9
j_s, t_s, q3_s, p3_s, E_s = np.loadtxt(data_directory + 'simplectic_results_911.txt', unpack=True)
j_rk, t_rk, q3_rk, p3_rk, E_rk = np.loadtxt(data_directory + 'rk4_results_911.txt', unpack=True)

fig = plt.figure(figsize=(30, 30))

fig9a = fig.add_subplot(211)
plt.scatter(q3_s, p3_s, s=1)
plt.xlim(-3.2, 3.2)
plt.ylim(-2.5, 2.5)
plt.title('Fourth Order Symplectic Integration Algorithm', fontsize=30)
plt.xlabel('$q_3$', fontsize=40)
plt.ylabel('$p_3$', fontsize=40)

fig9b = fig.add_subplot(212)
plt.scatter(q3_rk, p3_rk, s=1)
plt.xlim(-3.2, 3.2)
plt.ylim(-2.5, 2.5)
plt.title('Fourth Order Runge-Kutta Algorithm', fontsize=30)
plt.xlabel('$q_3$', fontsize=40)
plt.ylabel('$p_3$', fontsize=40)

plt.savefig(plots_directory + 'Figure9.png')

# plots Figure 10
fig = plt.figure(figsize=(30, 30))

fig10a = fig.add_subplot(211)
plt.scatter(q3_s, p3_s, s=1)
plt.xlim(-2.4, -0.9)
plt.ylim(-0.65, 0.65)
plt.title('Fourth Order Symplectic Integration Algorithm', fontsize=30)
plt.xlabel('$q_3$', fontsize=40)
plt.ylabel('$p_3$', fontsize=40)

fig10b = fig.add_subplot(212)
plt.scatter(q3_rk, p3_rk, s=1)
plt.xlim(-2.4, -0.9)
plt.ylim(-0.65, 0.65)
plt.title('Fourth Order Runge-Kutta Algorithm', fontsize=30)
plt.xlabel('$q_3$', fontsize=40)
plt.ylabel('$p_3$', fontsize=40)

plt.savefig(plots_directory + 'Figure10.png')

# plots Figure 11
fig = plt.figure(figsize=(30, 30))

fig11a = fig.add_subplot(211)
plt.scatter(q3_s, p3_s, s=1)
plt.xlim(-0.25, 0.9)
plt.ylim(-0.35, 0.35)
plt.title('Fourth Order Symplectic Integration Algorithm', fontsize=30)
plt.xlabel('$q_3$', fontsize=40)
plt.ylabel('$p_3$', fontsize=40)

fig11b = fig.add_subplot(212)
plt.scatter(q3_rk, p3_rk, s=1)
plt.xlim(-0.25, 0.9)
plt.ylim(-0.35, 0.35)
plt.title('Fourth Order Runge-Kutta Algorithm', fontsize=30)
plt.xlabel('$q_3$', fontsize=40)
plt.ylabel('$p_3$', fontsize=40)

plt.savefig(plots_directory + 'Figure11.png')

# plots Figure 12
j_12, t_12, q3_12, p3_12, E_12 = np.loadtxt(data_directory + 'simplectic_results_1214.txt', unpack=True)
fig = plt.figure(figsize=(40, 20))
plt.scatter(q3_12, p3_12, s=1)
plt.xlim(-2, 2)
plt.ylim(-2.5, 2.5)
plt.title('Orbits of the Poincare map $P_a$, $a = 0.4325$', fontsize=30)
plt.xlabel('$q_3$', fontsize=40)
plt.ylabel('$p_3$', fontsize=40)
plt.savefig(plots_directory + 'Figure12.png')

# plots Figure 13
fig = plt.figure(figsize=(40, 20))
plt.scatter(q3_12, p3_12, s=1)
plt.xlim(0.2, 0.7)
plt.ylim(-0.65, 0.65)
plt.title('Zooming in on the region near $(0.45,0)$ in Figure 12', fontsize=30)
plt.xlabel('$q_3$', fontsize=40)
plt.ylabel('$p_3$', fontsize=40)
plt.savefig(plots_directory + 'Figure13.png')

# plots Figure 14
fig = plt.figure(figsize=(40, 20))
plt.scatter(q3_12, p3_12, s=1)
plt.xlim(-0.08, 0.08)
plt.ylim(-0.1, 0.1)
plt.title('$P_a$-orbits in a neighborhood of the origin in Figure 12', fontsize=30)
plt.xlabel('$q_3$', fontsize=40)
plt.ylabel('$p_3$', fontsize=40)
plt.savefig(plots_directory + 'Figure14.png')

# plots Figure 15
j_15, t_15, q3_15, p3_15, E_15 = np.loadtxt(data_directory + 'simplectic_results_15.txt', unpack=True)
fig = plt.figure(figsize=(40, 20))
plt.scatter(q3_15, p3_15, c=j_15, s=1, edgecolors='none')
plt.xlim(-3.2, 3.2)
plt.ylim(-2.5, 2.5)
plt.title('Ergodic orbits of $P_a$ ($a = 0.425$)', fontsize=30)
plt.xlabel('$q_3$', fontsize=40)
plt.ylabel('$p_3$', fontsize=40)
plt.savefig(plots_directory + 'Figure15.png')

# plots Figure 16
j_16, t_16, q3_16, p3_16, E_16 = np.loadtxt(data_directory + 'simplectic_results_16.txt', unpack=True)
fig = plt.figure(figsize=(40, 20))
plt.scatter(q3_16, p3_16, c=j_16, s=1, edgecolors='none')
plt.xlim(-3.2, 3.2)
plt.ylim(-2.5, 2.5)
plt.title('Ergodic orbits of $P_a$ ($a = 0.46$)', fontsize=30)
plt.xlabel('$q_3$', fontsize=40)
plt.ylabel('$p_3$', fontsize=40)
plt.savefig(plots_directory + 'Figure16.png')

# plots Figure 17 (using data from Figure 9)
selected_indexes_s = np.where(j_s == 0) 
selected_indexes_rk = np.where(j_rk == 0) 

fig = plt.figure(figsize=(20, 20))

fig17a = fig.add_subplot(211)
plt.plot(t_s[selected_indexes_s], E_s[selected_indexes_s])
plt.xlim(0, 2600)
plt.title('Energy Conservation in Symplectic Integration', fontsize=30)
plt.xlabel('$t$', fontsize=40)
plt.ylabel('$E$', fontsize=40)

fig17b = fig.add_subplot(212)
plt.plot(t_rk[selected_indexes_rk], E_rk[selected_indexes_rk])
plt.xlim(0, 2600)
plt.title('Energy Conservation in Runge-Kutta', fontsize=30)
plt.xlabel('$t$', fontsize=40)
plt.ylabel('$E$', fontsize=40)

plt.savefig(plots_directory + 'Figure17.png')

