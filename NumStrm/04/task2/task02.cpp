#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include "../../../lib/linear_system.h"
#include "../../../lib/ode.h"
#include "../../../lib/square_matrix.h"


using namespace std;
using namespace numerical;

double inhom_q(double x)
{
  return 0;
}

//returns the error after one step
void btcs_time_step(double* values, square_matrix& coeff_matrix, double t, double delta_t, double delta_x, size_t dimension_x)
{
  double* x = new double[coeff_matrix.get_n()];
  double* b = new double[coeff_matrix.get_n()];

  double s = delta_t/pow(delta_x,2);

  b[0] = -values[1]-0;
  b[coeff_matrix.get_n()-1] = -values[dimension_x-2]-s*sin(10 * M_PI * t);

  for (size_t i=1; i<coeff_matrix.get_n()-1; i++)
  {
    b[i] = -values[i+1];
  }

  linear_system_solve_triagonal(coeff_matrix, x, b);

  values[0] = 0;
  values[dimension_x-1] = sin(10 * M_PI * t);
  for (size_t i=0; i<coeff_matrix.get_n(); i++)
  {
    values[i+1] = x[i];
  }

  delete x;
  delete b;
}

//returns the error after one step
double ftcs_time_step(double* values, double t, double delta_t, double delta_x, size_t dimension_x)
{
  double* old_values = new double[dimension_x];
  for (size_t x=0; x<dimension_x; x++)
  {
    old_values[x] = values[x];
  }

  values[dimension_x-1] = sin(10 * M_PI * t);

  for (size_t x=1; x<dimension_x-1; x++)
  {
    values[x] += delta_t * (
        (old_values[x+1] + old_values[x-1] - 2*old_values[x])/(pow(delta_x,2))
        + (inhom_q(x*delta_x))
      );
  }

  double error = 0.0;

  for (size_t x=1; x<dimension_x-1; x++)
  {
    error += pow(values[x] - old_values[x],2);
  }
  delete old_values;

  return sqrt(error);
}

int main(int argc, char* argv[])
{
  ofstream outputFile;
	outputFile.open(argv[1]);
	outputFile << fixed << setprecision(10);
  cout << fixed << setprecision(20);

  size_t dimension_x = atoi(argv[2]);
  double delta_x = 1.0/dimension_x;

  double max_time = atof(argv[3]);
  double delta_t = atof(argv[4]);
  int max_time_steps = max_time / delta_t;

  double* values = new double[dimension_x];
  for (size_t x=0; x<dimension_x; x++)
  {
    values[x] = 0.0;
  }

/*  double error;
  for (int time_step=0; time_step < max_time_steps; time_step++)
  {
    error = ftcs_time_step(values, time_step * delta_t, delta_t, delta_x, dimension_x);
  }*/

  //begin: implicit calculation
  double s = delta_t/pow(delta_x,2);

  square_matrix coeff_matrix = square_matrix(dimension_x-2,0);
  coeff_matrix.set_value(0,0,-(1+2*s));
  coeff_matrix.set_value(0,1,s);

  for (size_t i=1; i<dimension_x-3; i++)
  {
    coeff_matrix.set_value(i,i,-(1+2*s));
    coeff_matrix.set_value(i,i-1,s);
    coeff_matrix.set_value(i,i+1,s);
  }

  coeff_matrix.set_value(dimension_x-3,dimension_x-4,s);
  coeff_matrix.set_value(dimension_x-3,dimension_x-3,-(1+2*s));

//coeff_matrix.print_to_stream(cout);

  cout << "Setup done" << endl << flush;
  for (int time_step=0; time_step < max_time_steps; time_step++)
  {
    btcs_time_step(values, coeff_matrix, time_step * delta_t, delta_t, delta_x, dimension_x);
  }
  //end: implicit calculation

  //begin: explicit calculation
  /*
  for (int time_step=0; time_step < max_time_steps; time_step++)
  {
    ftcs_time_step(values, time_step * delta_t, delta_t, delta_x, dimension_x);
  }
  */
  //end: explicit calculation

  //btcs_time_step(values, coeff_matrix, 0, delta_t, delta_x, dimension_x);
  //btcs_time_step(values, coeff_matrix, delta_t, delta_t, delta_x, dimension_x);
  //cout << values[0] << ";" << values[1] << ";" << values[2] << ";" << values[3] << ";" << values[4] << endl;
  for (size_t x=0; x<dimension_x; x++)
  {
    outputFile << x * delta_x << "\t" << values[x] << endl;
  }
  outputFile << flush;
}

//use s=dt/dx^2 <= 1/2 for s min => dt_min = 0.005 = 5e-3
