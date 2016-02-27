all: simplectic_results_911.txt rk4_results_911.txt 

# Results of simplectic integration for plots 9 to 11
simplectic_results_911.txt: 3bodiesSimpl.x initialConditions911.txt 
	./3bodiesSimpl.x 0.006 2800 initialConditions911.txt 127 simplectic_results_911.txt

3bodiesSimpl.x: 3bodiesSimpl.c
	cc 3bodiesSimpl.c -o 3bodiesSimpl.x

# Results of runge kutta integration for plots 9 to 11
rk4_results_911.txt: 3bodiesRK4.x initialConditions911.txt 
	./3bodiesRK4.x 0.006 2800 initialConditions911.txt 127 rk4_results_911.txt

3bodiesRK4.x: 3bodiesRK4.c
	cc 3bodiesRK4.c -o 3bodiesRK4.x

clean: 
	rm -f simplectic_results.txt rk4_results.txt simplectic_results_911.txt rk4_results_911.txt 3bodiesRK4.x 3bodiesSimpl.x
