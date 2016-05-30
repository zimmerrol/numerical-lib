#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "math.h"
#include "../../../lib/ode.h"

using namespace std;

const double epsilon = 1e-2;
const double density_gamma = 4.0/3.0;
const double k = 3.85e9;
const double G = 6.67384e-11;
const double r_max = 7e8;

double f_d_rho(double* args, double* params)
{
  const double r = params[0];
  const double rho = args[0];
  const double s = args[1];

  return s;
}

double f_d_s(double* args, double* params)
{
  const double r = params[0];
  const double rho = args[0];
  const double s = args[1];

  /*return -4*PI*G / (
                      k * pow(rho, density_gamma-3)
                    )
          - s/r * (
                    2+r*(density_gamma-2)*s/rho
                  );*/

    return (-4*M_PI*pow(rho,3-density_gamma)*G*pow(r,2) - k*pow(r,2)*density_gamma*(density_gamma-2.0)*pow(rho, -1)*pow(s,2) - k*2*r*density_gamma*s)/(pow(r,2)*density_gamma*k);

}

int main(int argc, char* argv[])
{
  ofstream outputFile;
  outputFile.open(argv[1]);
  outputFile << fixed << setprecision(17);

  double delta_r = atof(argv[2]);
  double rho0 = atof(argv[3]);

  double delta_rho0 = 1.0;
  double epsilon_rho0 = 0.01;

  double epsilon_rho = 2e-4;

  numerical::odeFunction* functions = new numerical::odeFunction[2];
  functions[0] = &f_d_rho;
  functions[1] = &f_d_s;

  double* rkValues = new double[2];
  rkValues[0] = rho0;
  rkValues[1] = -1/density_gamma * pow(rho0,2-density_gamma)*G*4/3*M_PI*rho0*epsilon; //-1/(density_gamma * k) * pow(rho0, 3-density_gamma)*G*PI*epsilon*4/3;

  double* params = new double[1];

  //determine real rho0
  double real_r = 0.0;
  while (true)
  {
    //real_r = r_max;
    //break;

    cout << "rho0 := " << rho0 << endl;
    rho0 += delta_rho0;

    rkValues[0] = rho0;
    rkValues[1] = -1/density_gamma * pow(rho0,2-density_gamma)*G*4/3*M_PI*rho0*epsilon; //-1/(density_gamma * k) * pow(rho0, 3-density_gamma)*G*PI*epsilon*4/3;

    for(double r=epsilon + delta_r; r<=r_max *2 ; r+=delta_r)
    {
      params[0] = r;

      numerical::step_rk4_explicit(2, functions, delta_r, rkValues, params);

      if (rkValues[0] < epsilon_rho)
      {
        cout << "\t" << "r = " << r << "\tdiff = " << r-r_max << "\t" << "rho = " << rkValues[0] << "\t delta_rho = " << delta_rho0 << endl;
        if (abs(r-r_max) < epsilon_rho0)
        {
          real_r = r;
        }

        break;
      }
    }

    if (real_r != 0.0)
    {
      break;
    }
    cout << "\t no result.\trho = " << rkValues[0] << endl;
  }

  cout << "found rho0 = " << rho0 << endl;

  rkValues[0] = rho0;
  rkValues[1] = -1/density_gamma * pow(rho0,2-density_gamma)*G*4/3*M_PI*rho0*epsilon; //-1/(density_gamma * k) * pow(rho0, 3-density_gamma)*G*PI*epsilon*4/3;
  outputFile << epsilon << "\t" << rkValues[0] << "\t" << rkValues[1] << "\n";//"\t" << xIm << "\t" << yIm  << "\n";

  int print_treshold = 100 * 1000/delta_r;
  int print = print_treshold;
  for(double r=epsilon + delta_r; r<=real_r; r+=delta_r)
  {
    params[0] = r;

    numerical::step_rk4_explicit(2, functions, delta_r, rkValues, params);

    if (std::isnan(rkValues[0]) || std::isnan(rkValues[1]))
    {
      cout << "nan at r=" << r <<", rho = " << rkValues[0] << ", rho' = " << rkValues[1] << endl;
      return -1;
    }

    //if (print == print_treshold)
    {

      outputFile << r << "\t" << rkValues[0] << "\t" << rkValues[1] << "\n" << flush;//"\t" << xIm << "\t" << yIm  << "\n";
      print = 0;
    }
    print++;
  }
}
