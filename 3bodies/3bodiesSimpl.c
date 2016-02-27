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
void bigMass_simplectic_step(double t, double delta_t, double* q1, double* p1, double epsilon);
void littleMass_simplectic_step(double t, double delta_t, double* q1, double* p1, double* q3, double* p3, double epsilon);

int main(int argc, char **argv) {

  // Variables for simplectic integration.
  double q1_s, q2_s, q3_s;
  double p1_s, p2_s, p3_s;

  double t, epsilon;

  double T;
  double delta_t;
  int n_steps;
  int i, j;

  double previous_p1_s;

  FILE* input;
  int n_conditions;
  char line[256];

  FILE *simplectic_results;

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
  simplectic_results = fopen(argv[5], "w");

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
    initialConditions(line, &q1_s, &p1_s, &q2_s, &p2_s, &q3_s, &p3_s);

    for (i = 0 ; i < n_steps ; i++) {        
      previous_p1_s = p1_s;

      // Solve for large masses.
      bigMass_simplectic_step(t, delta_t, &q1_s, &p1_s, epsilon);

      q2_s = -q1_s;
      p2_s = -p1_s;

      // Solve for little mass.
      littleMass_simplectic_step(t, delta_t, &q1_s, &p1_s, &q3_s, &p3_s, epsilon);

      t += delta_t;

      // Print relevant results.
      if (p1_s == 0 || (previous_p1_s < 0 && p1_s > 0) ||  (previous_p1_s > 0 && p1_s < 0)) {
        fprintf(simplectic_results, "%d %f %.15e %.15e\n", j, t, q3_s, p3_s);
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

void bigMass_simplectic_step(double t, double delta_t, double* q1, double* p1, double epsilon) {
  double q1_in;
  double p1_in;

  q1_in = *q1;
  p1_in = *p1;

  p1_in += 0.5 * alpha_1 * deriv_p1(t, q1_in, p1_in, epsilon) * delta_t;
  q1_in += alpha_1 * deriv_q1(t, q1_in, p1_in, epsilon) * delta_t;
  p1_in += 0.5 * alpha_1 * deriv_p1(t, q1_in, p1_in, epsilon) * delta_t;

  p1_in += 0.5 * alpha_0 * deriv_p1(t, q1_in, p1_in, epsilon) * delta_t;
  q1_in += alpha_0 * deriv_q1(t, q1_in, p1_in, epsilon) * delta_t;
  p1_in += 0.5 * alpha_0 * deriv_p1(t, q1_in, p1_in, epsilon) * delta_t;

  p1_in += 0.5 * alpha_1 * deriv_p1(t, q1_in, p1_in, epsilon) * delta_t;
  q1_in += alpha_1 * deriv_q1(t, q1_in, p1_in, epsilon) * delta_t;
  p1_in += 0.5 * alpha_1 * deriv_p1(t, q1_in, p1_in, epsilon) * delta_t;

  *q1 = q1_in;
  *p1 = p1_in;
}

void littleMass_simplectic_step(double t, double delta_t, double* q1, double* p1, double* q3, double* p3, double epsilon) {
  double q3_in;
  double p3_in;

  q3_in = *q3;
  p3_in = *p3;

  p3_in += 0.5 * alpha_1 * deriv_p3(t, *q1, *p1, q3_in, p3_in, epsilon) * delta_t;
  q3_in += alpha_1 * deriv_q3(t, *q1, *p1, q3_in, p3_in, epsilon) * delta_t;
  p3_in += 0.5 * alpha_1 * deriv_p3(t, *q1, *p1, q3_in, p3_in, epsilon) * delta_t;

  p3_in += 0.5 * alpha_0 * deriv_p3(t, *q1, *p1, q3_in, p3_in, epsilon) * delta_t;
  q3_in += alpha_0 * deriv_q3(t, *q1, *p1, q3_in, p3_in, epsilon) * delta_t;
  p3_in += 0.5 * alpha_0 * deriv_p3(t, *q1, *p1, q3_in, p3_in, epsilon) * delta_t;

  p3_in += 0.5 * alpha_1 * deriv_p3(t, *q1, *p1, q3_in, p3_in, epsilon) * delta_t;
  q3_in += alpha_1 * deriv_q3(t, *q1, *p1, q3_in, p3_in, epsilon) * delta_t;
  p3_in += 0.5 * alpha_1 * deriv_p3(t, *q1, *p1, q3_in, p3_in, epsilon) * delta_t;

  *q3 = q3_in;
  *p3 = p3_in;
}