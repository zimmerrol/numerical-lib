#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "../../../lib/ode.h"

using namespace std;

//returns the error after one step
double ftcs_time_step(double* values, double delta_t, double delta_x, size_t dimension_x)
{
  double* old_values = new double[dimension_x];
  for (size_t x=0; x<dimension_x; x++)
  {
    old_values[x] = values[x];
  }

  for (size_t x=1; x<dimension_x-1; x++)
  {
    values[x] += delta_t * (
        (old_values[x+1] + old_values[x-1] - 2*old_values[x])/(pow(delta_x,2))
        + (sin(M_PI * x*delta_x))
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
	outputFile << fixed << setprecision(5);

  size_t dimension_x = atoi(argv[2]);
  double delta_x = 1.0/dimension_x;

  int max_time_steps = atoi(argv[3]);
  double delta_t = atof(argv[4]);

  double* values = new double[dimension_x];
  for (size_t x=0; x<dimension_x; x++)
  {
    values[x] = 0.0;
  }

  double error;
  for (int time_step=0; time_step < max_time_steps; time_step++)
  {
    error = ftcs_time_step(values, delta_t, delta_x, dimension_x);
    if (error < 0.00001)
    {
      cout << "Finished after " << time_step << " steps." << endl;
      break;
    }
  }

  for (size_t x=0; x<dimension_x; x++)
  {
    cout << x << "\t" << values[x] << endl;
  }
  cout << flush;
}

//use s=dt/dx^2 <= 1/2 for s min => dt_min = 0.005 = 5e-3
//pro zehnerprozent erhähen sich die rechnenschritte um den faktor 37
