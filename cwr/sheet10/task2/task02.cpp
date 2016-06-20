#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "../../../lib/ode.h"

using namespace std;

//use this to simulate the tsunami

void fcts_step(double* previous_x, double* x, double* beta, size_t dimension_x, double t, double T, double a0)
{
  double* new_x = new double[dimension_x];

//                                        0
  new_x[0] = t <= T ? a0*sin(2*M_PI/T*t) : 0;

  //                            -1
  for (size_t i=1; i<dimension_x; i++)
  {
      new_x[i] = 2*(1-pow(beta[i],2))*x[i]-previous_x[i] + pow(beta[i],2)*(x[i+1]+x[i-1]) + 1.0/4.0*(pow(beta[i+1],2) - pow(beta[i-1],2))*(x[i+1]-x[i-1]);
  }
  //new_x[dimension_x-1] = 0.0;
  if (t > T)
  {
      new_x[0] = new_x[0] + pow(beta[dimension_x-1],2)*(new_x[1]-new_x[0]);
  }
  new_x[dimension_x-1] = new_x[dimension_x-1] + pow(beta[dimension_x-1],2)*(new_x[dimension_x-1]-new_x[dimension_x-2]);
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
  cout << "#Usage: path, (x_length is fixed to 1000000), dimension_x, t_max, delta_t, a0, T" << endl;
  ofstream outputFile;
	outputFile.open(argv[1]);
	outputFile << fixed << setprecision(5);

  double x_length = 1000*1000;
  size_t dimension_x = atoi(argv[2]);
  double delta_x = x_length/dimension_x;
  double t_max = atof(argv[3]);
  double delta_t = atof(argv[4]);

  double a0 = atof(argv[5]);
  double T = atof(argv[6]);
  double omega = 2*M_PI/T;

  double x_reef_l = 300;
  double x_beach_l = 450;
  double x_beach_r = 550;
  double x_reef_r = 700;

  x_reef_l *= 1000/delta_x;
  x_beach_l *= 1000/delta_x;
  x_beach_r *= 1000/delta_x;
  x_reef_r *= 1000/delta_x;

  double h_ocean = 5000;
  double h_reef = 20;

  cout << "#x_length: " << x_length << "\tdimension_x" << dimension_x << "\tt_max: " << t_max << "\tdelta_t: " << delta_t << "\tA0: " << a0 << endl;

  double* beta = new double[dimension_x];
  //beta = u*delta_t/delta_x
  double h;
  for (size_t i=0; i<dimension_x; i++)
  {
      if (i < x_reef_l)
      {
         h = h_ocean;
      }
      else if (i < x_beach_l)
      {
        h = h_ocean + (h_reef - h_ocean)*pow(sin((M_PI*(i-x_reef_l))/(2*(x_beach_l-x_reef_l))),2);
      }
      else if (i < x_beach_r)
      {
        h = h_reef;
      }
      else if (i < x_reef_r)
      {
        h = h_reef - (h_reef - h_ocean)*pow(
            sin(
                (M_PI*(x_beach_r-i)) /
                (2*(x_beach_r-x_reef_r))
              )
          ,2);
      }
      else
      {
        h = h_ocean;
      }

      beta[i] = sqrt(h*9.81) * delta_t/delta_x;
      cout << i << "\t" << beta[i] << endl;
  }

  cout << "Setup finished. Starting simulation." << endl;

  double* x = new double[dimension_x];
  double* previous_x = new double[dimension_x];
  for (size_t i=0; i<dimension_x; i++)
  {
      x[i] = 0.0;
      previous_x[i] = 0.0;
  }

  for (double t=0; t<t_max; t+=delta_t)
  {
    if ((int)t % (int)(t_max / 1000.0) == 0)
    {
      cout << "t = " << t << endl;
      print(outputFile, x, delta_x, dimension_x);
    }
    outputFile << flush;
    fcts_step(x, previous_x, beta, dimension_x, t, T, a0);
  }
}

///usage: task02 res/out.dat 100 22000 1 0.5 2000
