#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "../../../lib/ode.h"

const double a0 = 1.0;
const double T = 2000;
const double omega = 2*M_PI/T;



using namespace std;

//use this to simulate the tsunami

void fcts_step(double* previous_x, double* x, double* beta, size_t dimension_x, double t)
{
  double* new_x = new double[dimension_x];

  new_x[0] = t <= T ? a0*sin(omega*t) : 0;

  for (size_t i=1; i<dimension_x-1; i++)
  {
      new_x[i] = 2*(1-pow(beta[i],2))*x[i]-previous_x[i]+pow(beta[i],2)*(x[i+1]+x[i-1])+1.0/4.0*(pow(beta[i+1],2)-pow(beta[i-1],2))*(x[i+1]-x[i-1]);
  }
  new_x[dimension_x-1] = 0.0;

  for (size_t i=0; i<dimension_x; i++)
  {
    previous_x[i] = x[i];
    x[i] = new_x[i];
  }
}

void print(ofstream& outputFile, double* x, double delta_x, size_t dimension_x)
{
  for (size_t i=0; i<dimension_x; i++)
  {
    outputFile << "\t" << i*delta_x << "\t" << x[i] << endl;
  }
}

int main(int argc, char* argv[])
{
  cout << "Usage: path, x_length, dimension_x, t_length, delta_t, v_const" << endl;
  ofstream outputFile;
	outputFile.open(argv[1]);
	outputFile << fixed << setprecision(5);

  double x_length = atof(argv[2]);
  size_t dimension_x = atoi(argv[3]);
  double delta_x = x_length/dimension_x;
  double t_max = atof(argv[4]);
  double delta_t = atof(argv[5]);

  double a0 = atof(argv[6]);

  cout << "x_length: " << x_length << "\tdimension_x" << dimension_x << "\tt_max: " << t_max << "\tdelta_t: " << delta_t << endl;

  double* beta = new double[dimension_x];
  for (size_t i=0; i<dimension_x; i++)
  {
      beta[i] = sqrt((2500*(tanh(4*(i*delta_x)/1000-17.2)+1)+0.05*1000)*9.81)*delta_t/delta_x;
      cout << beta[i] << endl;
  }

  double* x = new double[dimension_x];
  double* previous_x = new double[dimension_x];
  for (size_t i=0; i<dimension_x; i++)
  {
      x[i] = 0.0;
      previous_x[i] = 0.0;
  }

  for (double t=0; t<t_max; t+=delta_t)
  {
    print(outputFile, x, delta_x, dimension_x);
    outputFile << flush;
    fcts_step(x, previous_x, beta, dimension_x, t);
  }
}

///usage: task02 res/out.dat 5000 100 1000 0.1 0.1
