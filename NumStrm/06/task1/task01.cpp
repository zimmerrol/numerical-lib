#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "../../../lib/ode.h"

using namespace std;

double analytic_solution(double x, double D, double v)
{
  return (exp(v/D*x)-1)/(exp(v/D)-1);
}

//returns the error after one step
double ftcs_time_step(double* values, double delta_t, double delta_x, double D, double v, size_t dimension_x)
{
  double* old_values = new double[dimension_x];
  for (size_t i=0; i<dimension_x; i++)
  {
    old_values[i] = values[i];
  }

  double first_derivation, second_derivation;

  values[0] = 0;
  values[dimension_x-1] = 1;

  for (size_t i=1; i<dimension_x-1; i++)
  {
    first_derivation = (old_values[i+1]-old_values[i-1])/(2*delta_x);
    second_derivation = (old_values[i+1]-2*old_values[i]+old_values[i-1])/pow(delta_x,2);
    values[i] += D*delta_t*second_derivation - v*delta_t*first_derivation;
  }

  double error = 0.0;

  for (size_t i=1; i<dimension_x-1; i++)
  {
    error += pow(values[i] - old_values[i],2);
  }
  delete old_values;

  return sqrt(error);
}

int main(int argc, char* argv[])
{
  cout << "Expects:" << endl;
  cout << "\t" << "outputPath\tdimension_x\tmax_time\tdelta_t\tv\tD" << endl;
  ofstream outputFile;
	outputFile.open(argv[1]);
	outputFile << fixed << setprecision(5);

  size_t dimension_x = (atoi)(argv[2]);
  double delta_x = 1.0/(dimension_x-1);

  int max_time = (int)atof(argv[3]);
  double delta_t = atof(argv[4]);

  double v = atof(argv[5]);
  double D = atof(argv[6]);

  cout << "delta_t="<<delta_t << "\tmax_time=" << max_time << "\tdelta_x=" << delta_x << "\tv=" << v << "\tD=" << D << endl;

  double* values = new double[dimension_x];
  for (size_t i=0; i<dimension_x-1; i++)
  {
    values[i] = i*delta_x;
  }
  values[dimension_x-1] = 1.0;

  for (double t=0; t<max_time; t+= delta_t)
  {
    ftcs_time_step(values, delta_t, delta_x, D, v, dimension_x);
  }

  double diff = 0.0;
  for (size_t i=0; i<dimension_x; i++)
  {
    outputFile << i*delta_x << "\t" << values[i] << "\t" << analytic_solution(i*delta_x, D, v) << endl;
    diff += pow(analytic_solution(i*delta_x, D, v) - values[i],2);
  }

  cout << "MSE = " << sqrt(diff) << endl;
  cout << flush;
}
