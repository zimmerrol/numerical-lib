#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include "../lib/ode.h"

using namespace std;

double fX(double* args, double* params)
{
  return (998 * args[0]) + 1998 * args[1];
}

double fY(double* args, double* params)
{
  return (-999 * args[0]) - 1999 * args[1];
}


int main(int argc, char* argv[])
{
  double xExp = 1;
  double yExp = 0;

  ofstream outputFile;
  outputFile.open(argv[1]);
  outputFile << fixed << setprecision(9);

  double deltaT = atof(argv[2]);
  double* explicitValues = new double[2];
  explicitValues[0] = 1;
  explicitValues[1] = 0;
  numerical::odeFunction* functions = new numerical::odeFunction[2];
  functions[0] = &fX;
  functions[1] = &fY;

  for(double t =0;t <= 1;t+=deltaT)
  {
    outputFile << t << "\t" << explicitValues[0] << "\t" << explicitValues[1] << "\n";
    //numerical::stepEulerExplicit(2, functions, deltaT, explicitValues, NULL);
    numerical::step_rk4_explicit(2, functions, deltaT, explicitValues, NULL);
  }

}
