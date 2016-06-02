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

  //if we assume T'[0] == 0 in the "usual" equation we get:
  values[0] += delta_t * (2.0*old_values[1]-2.0*old_values[0])/(pow(delta_r,2));

  //if we use the neuman bound.cond. according to the lecture like: T'[0] = 1/h*(-3T[0] +4T[1] -1T[2]) and say T'[0]==0 =>
  //values[0] = (-values[2]+4.0*values[1])/3.0;
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
	cout << fixed << setprecision(2);

  size_t dimension_r = (atoi)(argv[1]);
  double earth_radius = (atof)(argv[2]);
  double delta_r = 1.0/dimension_r;

  double delta_t = atof(argv[3]);
  double earth_water_gradient = (atof)(argv[4]);
  double kappa = (atof)(argv[5]);

  ofstream outputFile;
  outputFile.open(argv[6]);
  outputFile << fixed << setprecision(5);

  cout << "delta_r= "<< delta_r <<"\tdelta_t="<<delta_t << "\tearth_water_gradient=" << earth_water_gradient << "\tkappa=" << kappa << endl;

  double* values = new double[dimension_r];
  for (size_t r=0; r<dimension_r-1; r++)
  {
    values[r] = 6000.0;
  }

  values[dimension_r-1] = 0;

  double error;
  double derivation;
  int time_step_write_counter = 100;
  for (double t=0; t < 1e18; t+=delta_t)
  {
    if (isnan(values[dimension_r-1]) || isnan(values[dimension_r-2]))
    {
      cout << "nan detected for t=" << t << endl;
      return -1;
    }

if (time_step_write_counter == 100)
{
    for (size_t i=0; i<dimension_r; i++)
    {
      outputFile << "\t" << i << "\t" << values[i] << endl;
    }
time_step_write_counter = 0;
}
time_step_write_counter++;

    error = ftcs_time_step(values, delta_t, delta_r, dimension_r);
    derivation = (values[dimension_r-1]-values[dimension_r-2])/delta_r;
    if (abs(derivation) < abs(earth_water_gradient))
    {
      cout << "Finished after " << t/delta_t << " steps." << endl;
      cout << "\tThe earth seems to have an age of " << t/kappa*pow(earth_radius,2) / (60*60*24*365) <<" years" << endl;
      return 1;
    }
  }

  cout << "Did not find the age, sorry. Stopped at derivation of: " << (values[dimension_r-1]-values[dimension_r-2])/delta_r << endl;
}

/*
bash-4.3$ ./task02 100 6.4e6 2e-6 -211200 1e-6
delta_r= 0.01	delta_t=0.00	earth_water_gradient=-211200.00	kappa=0.00
Finished after 116.00 steps.
	The earth seems to have an age of 301329274.48 years

  bash-4.3$ ./task02 10000 6.4e6 2e-9 -211200 1e-6
delta_r= 0.00	delta_t=0.00	earth_water_gradient=-211200.00	kappa=0.00
Finished after 121473.00 steps.
	The earth seems to have an age of 315546301.37 years

*/
