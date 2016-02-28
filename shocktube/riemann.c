#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "riemann.h"
#define MIN(x) x<0.0?x:0.0

/*

Referencias

Roe solver tomado del libro de Computational Gasdynamics de Culbert Laney

*/


void Riemann(double *U1, double *U4, double *F) {

  const double gamma = 1.4;

  //Calcular variables primitivas
  double rho1 = U1[0];
  double u1 = U1[1] / rho1;
  double p1 = (U1[2] - rho1 * u1 * u1 / 2.0) * (gamma - 1.0);
  double ht1 = (U1[2] + p1)/rho1;
  double rho4 = U4[0];
  double u4 = U4[1] / rho4;
  double p4 = (U4[2] - rho4 * u4 * u4 / 2.0) * (gamma - 1.0);
  double ht4 = (U4[2] + p4)/rho4;

  //Calcular cantidades promedio de Roe
  double rho14 = sqrt(rho1*rho4);
  double u14 = (sqrt(rho4) * u4 + sqrt(rho1) * u1)/(sqrt(rho4) + sqrt(rho1)); 
  double h14 = (sqrt(rho4) * ht4  + sqrt(rho1) * ht1)/(sqrt(rho4) + sqrt(rho1));
  double a14 = sqrt((gamma - 1)*(h14 - u14*u14/2.0)); //NaN se produce cuando es raiz de numero negativo.

  //Calcular velocidades de onda promedio
  double lambda1 = u14;
  double lambda2 = u14 + a14;
  double lambda3 = u14 - a14;

  //Variables delta
  double drho = rho4 - rho1;
  double dp = p4 - p1;
  double du = u4 - u1;

  //Calcular fuerzas de onda
  double dv1 = drho - dp/(a14*a14);
  double dv2 = du + dp/(rho14*a14);
  double dv3 = du - dp/(rho14*a14);

  //Calcular flujos
  
  double alpha = rho14/(2.0*a14);

  F[0] = rho1*u1 + (MIN(lambda1))*dv1 + alpha*((MIN(lambda2))*dv2 - (MIN(lambda3))*dv3);
  F[1] = rho1*u1*u1 + p1 + u14*(MIN(lambda1))*dv1 + alpha*(lambda2*(MIN(lambda2))*dv2 - lambda3*(MIN(lambda3))*dv3);
  F[2] = rho1*ht1*u1 + u14*u14*(MIN(lambda1))*dv1/2.0 + alpha*((h14 + a14*u14)*(MIN(lambda2))*dv2 - (h14 - a14*u14)*(MIN(lambda3))*dv3);
 

  

}
