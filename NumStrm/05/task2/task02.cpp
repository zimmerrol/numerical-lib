#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "../../../lib/ode.h"

using namespace std;

//returns the error after one step
double ftcs_time_step(double* values, double delta_t, double delta_r, size_t dimension_r)
{
  double* old_values = new double[dimension_r];
  for (size_t r=0; r<dimension_r; r++)
  {
    old_values[r] = values[r];
  }

  double first_derivation, second_derivation;

  values[0] += delta_t * (2*old_values[1]-2*old_values[0])/(pow(delta_r,2));
  for (size_t r=1; r<dimension_r-1; r++)
  {
    first_derivation = (old_values[r+1] - old_values[r-1])/(2*delta_r);
    second_derivation = (old_values[r+1] + old_values[r-1] - 2*old_values[r])/(pow(delta_r,2));
    values[r] += delta_t * (2.0*first_derivation/(r*delta_r) + second_derivation);
  }

  double error = 0.0;

  for (size_t r=1; r<dimension_r-1; r++)
  {
    error += pow(values[r] - old_values[r],2);
  }
  delete old_values;

  return sqrt(error);
}

int main(int argc, char* argv[])
{
  ofstream outputFile;
	outputFile.open(argv[1]);
	outputFile << fixed << setprecision(5);

  size_t dimension_r = (atoi)(argv[2]);
  double earth_radius = (atof)(argv[3]);
  double delta_r = 1/dimension_r;

  double delta_t = atof(argv[4]);
  double earth_water_gradient = (atof)(argv[5]);
  double kappa = (atof)(argv[6]);


  cout << "delta_t="<<delta_t << "\tearth_water_gradient=" << earth_water_gradient << "\tkappa=" << kappa << endl;

  double* values = new double[dimension_r];
  for (size_t r=0; r<dimension_r-1; r++)
  {
    values[r] = 6000.0;
  }

  values[dimension_r-1] = 0;

  double error;
  double derivation;
  for (int time_step=0; time_step < 1e6; time_step++)
  {
    error = ftcs_time_step(values, delta_t, delta_r, dimension_r);
    derivation = (values[dimension_r-1]-values[dimension_r-2])/delta_r;
    if (abs(derivation) < abs(earth_water_gradient))
    {
      cout << "Finished after " << time_step << " steps." << endl;
      cout << "\tThe earth seems to have an age of " << (time_step * delta_t)*kappa/pow(earth_radius,2) <<" time units" << endl;
      return 1;
    }
  }

  cout << "Did not find the age, sorry. Stopped at derivation of: " <<   (values[dimension_r-1]-values[dimension_r-2])/delta_r << endl;
}