#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "math.h"
#include "../../lib/ode.h"

using namespace std;


/*
  structure:

  args contains the variables which will be integrated, so here omega, theta
    args[0] = omega
    args[1] = theta

  params contains the paramaters of the integration, so here t, q, f, big omega
    params[0] = t
    params[1] = q
    params[2] = f
    params[3] = big omega
*/

//returns dw/dt for the givven set of variables and parameters according to the structure above
double f_d_omega(double* args, double* params)
{
  const double omega = args[0];
  const double theta = args[1];

  const double t = params[0];
  const double q = params[1];
  const double f = params[2];
  const double big_omega = params[3];

  return -sin(theta) - q*omega + f*cos(big_omega * t);
}

//returns dtheta/dt for the givven set of variables and parameters according to the structure above
double f_d_theta(double* args, double* params)
{
  const double omega = args[0];

  return omega;
}

int main(int argc, char* argv[])
{
  cout << "expects the parameters:" << endl;

  cout << "\toutput file path, q, big omega, f, maximum time, time step size" << endl;
  //read input parameters
  char* file_path = argv[1];
  double q = atof(argv[2]);
  double big_omega = atof(argv[3]);
  double f = atof(argv[4]);
  double max_t = atof(argv[5]);
  double delta_t = atof(argv[6]);

  double start_print_t = max_t/2;


  //create output file stream
  ofstream outputFile;
  outputFile.open(file_path);
  //outputFile << fixed << setprecision(17);
  outputFile << "#t" << "\t" << "omega" << "\t" << "theta" << endl;


  //set the boundary conditions of the problem
  double theta_null = 0.0;
  double omega_null = big_omega;


  //create internal data for the RK4 integration implementation
  //set system of ODEs
  numerical::odeFunction* functions = new numerical::odeFunction[2];
  functions[0] = &f_d_omega;
  functions[1] = &f_d_theta;


  //set start values
  double* rkValues = new double[2];
  rkValues[0] = omega_null;
  rkValues[1] = theta_null;

  //set the other parameters of the problem for the evaluation of the ODEs
  double* rkParams = new double[4];
  //rkParams[0] will contain the current time t
  rkParams[1] = q;
  rkParams[2] = f;
  rkParams[3] = big_omega;

  double print_interval_t = 8.0/3.0 * M_PI;
  double last_printed_t = -print_interval_t;
  for (double t=0; t<=max_t; t+=delta_t)
  {


    //print only every 100th step for debug purposes
    //if (integrationSteps == 100)

    if (t - last_printed_t >= print_interval_t)
    {
      last_printed_t = t;
      //print the data
      outputFile << t << "\t" << rkValues[0] << "\t" << rkValues[1] << endl;
    }

    //set the time parameter
    rkParams[0] = t;

    //do one explicit RK4 integration step
    numerical::step_rk4_explicit(2, functions, delta_t, rkValues, rkParams);

    while (rkValues[1] > 2*M_PI)
    {
      rkValues[1] -= 2*M_PI;
    }

    while (rkValues[1] < -0)
    {
      rkValues[1] += 2*M_PI;
    }
  }
}
