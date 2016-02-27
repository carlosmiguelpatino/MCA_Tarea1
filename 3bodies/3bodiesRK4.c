#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define USAGE "./a.out delta_t T initialConditions.txt n_conditions results.txt\n"
#define alpha_0 -1*pow(2.0, 1.0/3.0)/(2.0 - pow(2.0, 1.0/3.0))
#define alpha_1 1.0/(2.0 - pow(2.0, 1.0/3.0))

void initialConditions(char* line, double* q1, double* p1, double* q2, double* p2, double* q3, double* p3);
char* getField(char* line, int num);
double deriv_q1(double t, double q1, double p1, double epsilon);
double deriv_p1(double t, double q1, double p1, double epsilon);
double deriv_q3(double t, double q1, double p1, double q3, double p3, double epsilon);
double deriv_p3(double t, double q1, double p1, double q3, double p3, double epsilon);
void bigMass_RK4_step(double t, double delta_t, double *q1, double *p1, double epsilon);
void littleMass_RK4_step(double t, double delta_t, double *q1, double *p1, double *q3, double *p3, double epsilon);
double calculateEnergy(double q1, double p1, double q2, double p2, double epsilon);

int main(int argc, char **argv) {

  // Variables for RK4 integration.
  double q1_rk, q2_rk, q3_rk;
  double p1_rk, p2_rk, p3_rk;

  double t, epsilon;

  double T;
  double delta_t;
  int n_steps;
  int i, j;

  double previous_p1_rk;
  double energy;

  FILE* input;
  int n_conditions;
  char line[256];

  FILE *rk4_results;

  // Without the correct number of attributes, the algorithm aborts.
  if(argc != 6)
  {
    printf("Usage: %s", USAGE);
    exit(1);
  }

  delta_t = atof(argv[1]);
  T = atof(argv[2]);
  input = fopen(argv[3], "r");
  n_conditions = atof(argv[4]);
  rk4_results = fopen(argv[5], "w");

  if(!input)
     exit(1);
  // Read header.
  fscanf(input, "%s", line); 

  n_steps = (int)(T/delta_t);
  t = 0.0;

  epsilon = 1.0;
  
  for (j = 0 ; j < n_conditions ; j++) {

    // Read initial conditions.
    fscanf(input, "%s", line); 
    //printf("%d\n", j);
    
    // Extract conditions.
    initialConditions(line, &q1_rk, &p1_rk, &q2_rk, &p2_rk, &q3_rk, &p3_rk);

    for (i = 0 ; i < n_steps ; i++) {        
      previous_p1_rk = p1_rk;

      // Solve for large masses.
      bigMass_RK4_step(t, delta_t, &q1_rk, &p1_rk, epsilon);

      q2_rk = -q1_rk;
      p2_rk = -p1_rk;

      // Solve for little mass.
      littleMass_RK4_step(t, delta_t, &q1_rk, &p1_rk, &q3_rk, &p3_rk, epsilon);

      t += delta_t;
      energy = calculateEnergy(q1_rk, p1_rk, q2_rk, p2_rk, epsilon);

      if (p1_rk == 0 || (previous_p1_rk < 0 && p1_rk > 0) ||  (previous_p1_rk > 0 && p1_rk < 0)) {
        fprintf(rk4_results, "%d %f %.15e %.15e %.15e\n", j, t, q3_rk, p3_rk, energy);
      }
    }
  }

  return 0;
}

void initialConditions(char* line, double* q1, double* p1, double* q2, double* p2, double* q3, double* p3) {
    char* tmp;

    *q1 = atof(getField(strdup(line), 1));
    *p1 = atof(getField(strdup(line), 2));

    *q2 = -(*q1);
    *p2 = -(*p1);

    *q3 = atof(getField(strdup(line), 3));
    *p3 = atof(getField(strdup(line), 4));
}

char* getField(char* line, int num)
{
    char* tok;
    for (tok = strtok(line, ";"); tok && *tok ; tok = strtok(NULL, ";\n"))
    {
      if (!--num)
        return tok;
    }
    printf("%s\n", "tokenizing failed");
    return NULL;
}

double deriv_q1(double t, double q1, double p1, double epsilon) {
  return p1;
}

double deriv_p1(double t, double q1, double p1, double epsilon) {
  double num;
  double radicand;
  double denom;

  num = -2.0*q1;
  radicand = 4.0*q1*q1 + epsilon*epsilon;
  denom = pow(radicand, 3.0/2.0);
  return num/denom;
}

double deriv_q3(double t, double q1, double p1, double q3, double p3, double epsilon) {
  return p3;
}

double deriv_p3(double t, double q1, double p1, double q3, double p3, double epsilon) {
  double num1, num2;
  double radicand1, radicand2;
  double denom1, denom2;
  double term1, term2;

  num1 = q1 - q3;
  radicand1 = (q1 - q3)*(q1 - q3) + 0.25*epsilon*epsilon;
  denom1 = pow(radicand1, 3.0/2.0);
  term1 = num1/denom1;

  num2 = q1 + q3;
  radicand2 = (q1 + q3)*(q1 + q3) + 0.25*epsilon*epsilon;
  denom2 = pow(radicand2, 3.0/2.0);
  term2 = num2/denom2;

  return term1 - term2;
}

void bigMass_RK4_step(double t, double delta_t, double *q1, double *p1, double epsilon) {
  double k1, k2, k3, k4;
  double l1, l2, l3, l4;
  double q1_in;
  double p1_in;

  q1_in = *q1;
  p1_in = *p1;

  k1 = deriv_q1(t, q1_in, p1_in, epsilon);  
  l1 = deriv_p1(t, q1_in, p1_in, epsilon);  

  k2 = deriv_q1(t + 0.5*delta_t, q1_in + 0.5*k1*delta_t, p1_in + 0.5*l1*delta_t, epsilon);
  l2 = deriv_p1(t + 0.5*delta_t, q1_in + 0.5*k1*delta_t, p1_in + 0.5*l1*delta_t, epsilon);

  k3 = deriv_q1(t + 0.5*delta_t, q1_in + 0.5*k2*delta_t, p1_in + 0.5*l2*delta_t, epsilon);
  l3 = deriv_p1(t + 0.5*delta_t, q1_in + 0.5*k2*delta_t, p1_in + 0.5*l2*delta_t, epsilon);

  k4 = deriv_q1(t + delta_t, q1_in + k3*delta_t, p1_in + l3*delta_t, epsilon);
  l4 = deriv_p1(t + delta_t, q1_in + k3*delta_t, p1_in + l3*delta_t, epsilon);

  q1_in += (k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0)*delta_t;
  p1_in += (l1/6.0 + l2/3.0 + l3/3.0 + l4/6.0)*delta_t;

  *q1 = q1_in;
  *p1 = p1_in;
}

void littleMass_RK4_step(double t, double delta_t, double *q1, double *p1, double *q3, double *p3, double epsilon) {
  double k1, k2, k3, k4;
  double l1, l2, l3, l4;
  double q3_in;
  double p3_in;

  q3_in = *q3;
  p3_in = *p3;

  k1 = deriv_q3(t, *q1, *p1, q3_in, p3_in, epsilon);  
  l1 = deriv_p3(t, *q1, *p1, q3_in, p3_in, epsilon);  

  k2 = deriv_q3(t + 0.5*delta_t, *q1, *p1, q3_in + 0.5*k1*delta_t, p3_in + 0.5*l1*delta_t, epsilon);
  l2 = deriv_p3(t + 0.5*delta_t, *q1, *p1, q3_in + 0.5*k1*delta_t, p3_in + 0.5*l1*delta_t, epsilon);

  k3 = deriv_q3(t + 0.5*delta_t, *q1, *p1, q3_in + 0.5*k2*delta_t, p3_in + 0.5*l2*delta_t, epsilon);
  l3 = deriv_p3(t + 0.5*delta_t, *q1, *p1, q3_in + 0.5*k2*delta_t, p3_in + 0.5*l2*delta_t, epsilon);

  k4 = deriv_q3(t + delta_t, *q1, *p1, q3_in + k3*delta_t, p3_in + l3*delta_t, epsilon);
  l4 = deriv_p3(t + delta_t, *q1, *p1, q3_in + k3*delta_t, p3_in + l3*delta_t, epsilon);

  q3_in += (k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0)*delta_t;
  p3_in += (l1/6.0 + l2/3.0 + l3/3.0 + l4/6.0)*delta_t;

  *q3 = q3_in;
  *p3 = p3_in;
}

// The article states that only the "big masses" system is Hamiltonian. 
// This system's energy is then the one we focus on.
double calculateEnergy(double q1, double p1, double q2, double p2, double epsilon) {
  double K = 0;
  double U = 0;
  int i,j;

  K += 0.5*p1*p1;
  K += 0.5*p2*p2;

  U += -1.0/(2.0*sqrt(4*q1*q1 + epsilon));
  U += -1.0/(2.0*sqrt(4*q2*q2 + epsilon)); 

  return K + U;
}
