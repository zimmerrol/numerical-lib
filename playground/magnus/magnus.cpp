#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "../../lib/ode.h"

const double m = 0.5;
const double g = 9.81;
const double cw = 0.5;
const double cm = 0.5;
const double rho = 1.293;
const double A = 0.0380133;
const double r = 0.11;
const double wx = +4*M_PI*8;
const double wy = -2*M_PI*8;
const double wz = -2*M_PI*8;

using namespace std;

double fRX(double* args, double* params)
{
  return args[3];
}

double fRY(double* args, double* params)
{
  return args[4];
}

double fRZ(double* args, double* params)
{
  return args[5];
}

double fVX(double* args, double* params)
{
  return 1/m*(-cw*M_PI*r*r/2*sqrt(args[3]*args[3] + args[4]*args[4] + args[5]*args[5]) * args[3] + cm*rho*M_PI*r*r*r/2*(wy*args[5]-wz*args[4]));
}

double fVY(double* args, double* params)
{
  return 1/m*(-cw*M_PI*r*r/2*sqrt(args[3]*args[3] + args[4]*args[4] + args[5]*args[5]) * args[4] + cm*rho*M_PI*r*r*r/2*(wz*args[3]-wx*args[5]));
}

double fVZ(double* args, double* params)
{
  return 1/m*(-m*g-cw*M_PI*r*r/2*sqrt(args[3]*args[3] + args[4]*args[4] + args[5]*args[5]) * args[5] + cm*rho*M_PI*r*r*r/2*(wx*args[4]-wy*args[3]));
}


int main(int argc, char* argv[])
{

  ofstream outputFile;
  outputFile.open(argv[1]);
  outputFile << fixed << setprecision(9);

  double deltaT = 0.001;
  double v = 20.8333;
  double phi = M_PI/9;
  double theta =(409 * M_PI)/4500 ;

  double* explicitValues = new double[6];
  explicitValues[0] = 0;
  explicitValues[1] = 0;
  explicitValues[2] = 0;
  explicitValues[3] = cos(theta) * cos(phi) * v;
  explicitValues[4] = sin(theta) * cos(phi) * v;
  explicitValues[5] = sin(phi) * v;
  numerical::odeFunction* functions = new numerical::odeFunction[6];
  functions[0] = &fRX;
  functions[1] = &fRY;
  functions[2] = &fRZ;
  functions[3] = &fVX;
  functions[4] = &fVY;
  functions[5] = &fVZ;

  for(double t =0;t <= 2;t+=deltaT)
  {
    outputFile << t << "\t" << explicitValues[0] << "\t" << explicitValues[1] << "\t" << explicitValues[2] << "\n";
    if (explicitValues[2] < 0)
    {
       break;
    }
    //numerical::stepEulerExplicit(2, functions, deltaT, explicitValues, NULL);
    numerical::stepRK2Explicit(6, functions, deltaT, explicitValues, NULL);

  }

}
