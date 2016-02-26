#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "riemann.h"
#define MIN(x,y) x<y?x:y

void compute_r_vectors(double *r1, double *r2, double *r3, double u14, double rho14, double a14, double h14);


void Riemann(double *U1, double *U4, double *F) {

  const double gamma = 1.4;
  
  // Calculo de variable basicas (paso 1)
  double rho1 = U1[0];
  double u1 = U1[1] / rho1;
  double p1 = (U1[2] - rho1 * u1 * u1 / 2) * (gamma - 1);
  double ht1 = (U1[2] + p1) / rho1;
  double rho4 = U4[0];
  double u4 = U4[1] / rho4;
  double p4 = (U4[2] - rho4 * u4 * u4 / 2) * (gamma - 1);
  double ht4 = (U4[2] + p4) / rho4;
  
  //Variables delta
  double delta_u = u4 - u1;
  double delta_rho = rho4 - rho1;
  double delta_p = p4 - p1;
  
  // Calcular cantidades promedio de Roe (paso 2)
  double rho14 = sqrt(rho1*rho4);
  double u14 = (u4 * sqrt(rho4) + u1 * sqrt(rho1))/(sqrt(rho1) + sqrt(rho4));
  double h14 = (ht4 * sqrt(rho4) + ht1 * sqrt(rho1))/(sqrt(rho1) + sqrt(rho4));
  double a14 = sqrt((gamma - 1) * (h14 - (1.0/2.0) * u14 * u14));
  //Calcular velocidad de onda de promedio de Roe (paso 3)
  double lambda1 = u14;
  double lambda2 = u14 + a14;
  double lambda3 = u14 - a14;

  //Calcular fuerzas de onda (paso 4) 
  double delta_v1 = delta_rho - delta_p/(a14 * a14);
  double delta_v2 = delta_u + delta_p/(rho14 * a14);
  double delta_v3 = delta_u - delta_p/(rho14 * a14);

  //Calcular los vectores caracteristicos (paso 5)
  double *r1, *r2, *r3;

  r1 = malloc(3*sizeof(double));  
  r2 = malloc(3*sizeof(double));  
  r3 = malloc(3*sizeof(double));

  compute_r_vectors(r1, r2, r3, u14, rho14, a14, h14);

  //Calcular flujos (paso 6)
  double f1, f2, f3;
  f1 = rho1 * u1 + r1[0]*MIN(0, lambda1)*delta_v1 + r2[0] * MIN(0,  lambda2)*delta_v2 + r3[0]*MIN(0, lambda3)*delta_v3;
  f2 = rho1 * u1 * u1 + p1 + r1[1]*MIN(0, lambda1)*delta_v1 + r2[1] * MIN(0,  lambda2)*delta_v2 + r3[1]*MIN(0, lambda3)*delta_v3;
  f3 = rho1 * u1 * ht1 + r1[2]*MIN(0, lambda1)*delta_v1 + r2[2] * MIN(0,  lambda2)*delta_v2 + r3[2]*MIN(0, lambda3)*delta_v3;
 

  //Retorno de flujos calculados
  F[0] = f1;
  F[1] = f2;
  F[2] = f3;
}


void compute_r_vectors(double *r1, double *r2, double *r3, double u14, double rho14, double a14, double h14){
  
 
  //Inicializar vector 1
  r1[0] = 1;
  r1[1] = u14;
  r1[2] = (1.0/2.0) * u14 * u14;

  //Inicializar vector 2
  r2[0] = rho14/(2 * a14);    
  r2[1] = (rho14/(2 * a14)) * (u14 + a14);    
  r2[2] = rho14/(2 * a14) * (h14 + a14 * u14); 

  //Inicialiar vector 3
 
  r3[0] = -rho14/(2 * a14);    
  r3[1] = (-rho14/(2 * a14)) * (u14 - a14);    
  r3[2] = (-rho14/(2 * a14)) * (h14 - a14 * u14);

   
} 












