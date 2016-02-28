all: a.out *.dat *.png

#Creacion solucion analitica y graficas

*.png: *.dat sol_analitica.py
	python sol_analitica.py

#Resultados de solucion numerica
*.dat: a.out
	./a.out

#Compilacion archivos en C

a.out: shocktube.c steps.c riemann.c
	gcc shocktube.c steps.c riemann.c -lm

clean:
	rm -f *.dat *.png a.out