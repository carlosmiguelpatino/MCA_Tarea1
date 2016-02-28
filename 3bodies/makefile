all: ./Figures/Figure9.png ./Figures/Figure10.png ./Figures/Figure11.png ./Figures/Figure12.png ./Figures/Figure13.png ./Figures/Figure14.png ./Figures/Figure15.png ./Figures/Figure16.png ./Figures/Figure17.png 

# Plots Figure 9
./Figures/Figure9.png ./Figures/Figure10.png ./Figures/Figure11.png ./Figures/Figure12.png ./Figures/Figure13.png ./Figures/Figure14.png ./Figures/Figure15.png ./Figures/Figure16.png ./Figures/Figure17.png: ./Results/simplectic_results_911.txt ./Results/simplectic_results_1214.txt ./Results/simplectic_results_15.txt ./Results/simplectic_results_16.txt ./Results/rk4_results_911.txt
	python plotter.py

# Results of Simplectic Integration for Figures 9 to 11.
./Results/simplectic_results_911.txt: 3bodiesSimpl.x ./InitialConditions/initialConditions911.txt 
	./3bodiesSimpl.x 0.006 2800 ./InitialConditions/initialConditions911.txt 127 ./Results/simplectic_results_911.txt

# Results of Simplectic Integration for Figures 12 to 14.
./Results/simplectic_results_1214.txt: 3bodiesSimpl.x ./InitialConditions/initialConditions1214.txt 
	./3bodiesSimpl.x 0.006 2800 ./InitialConditions/initialConditions1214.txt 84 ./Results/simplectic_results_1214.txt

# Results of Simplectic Integration for Figures 15.
./Results/simplectic_results_15.txt: 3bodiesSimpl.x ./InitialConditions/initialConditions15.txt 
	./3bodiesSimpl.x 0.006 1800 ./InitialConditions/initialConditions15.txt 114 ./Results/simplectic_results_15.txt

# Results of Simplectic Integration for Figures 16. 
./Results/simplectic_results_16.txt: 3bodiesSimpl.x ./InitialConditions/initialConditions16.txt 
	./3bodiesSimpl.x 0.006 2800 ./InitialConditions/initialConditions16.txt 114 ./Results/simplectic_results_16.txt

# Results of Runge-Kutta integration for Figures 9 to 11.
./Results/rk4_results_911.txt: 3bodiesRK4.x ./InitialConditions/initialConditions911.txt 
	./3bodiesRK4.x 0.006 2800 ./InitialConditions/initialConditions911.txt 127 ./Results/rk4_results_911.txt

# Compiles C code.
3bodiesSimpl.x: 3bodiesSimpl.c
	cc 3bodiesSimpl.c -o 3bodiesSimpl.x

3bodiesRK4.x: 3bodiesRK4.c
	cc 3bodiesRK4.c -o 3bodiesRK4.x

clean: 
	rm -r -f Results
	rm -r -f Figures
	rm -f 3bodiesRK4.x 3bodiesSimpl.x
	mkdir Results
