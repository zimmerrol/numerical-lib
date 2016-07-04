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
  //print information for the usage
  cout << "Use this to solve the system of ODEs and to save the calculated values of theta, omega for the drivven, damped pendulum." << endl;
  cout << "Expects the follwing set of parameters:" << endl;

  cout << "\toutput file:\tThe path to the file in which the results of the bifurcation will be stored." << endl;
  cout << "\tq:\t\tThe damping parameter." << endl;
  cout << "\tOmega:\t\tThe frequency of the driving force of the pendulum." << endl;
  cout << "\tf:\tThe value of the strength of the force f." << endl;
  cout << "\tmaximum time:\tThe maximum time for the integration (should be high to reach a fixed point)" << endl;
  cout << "\tx:\t\t Indicates the time step size via delta_t = 2pi/(Omega * x)]" << endl << endl;
  cout << "\tremap theta alt. flag:\t (default value = 0) If this is set to 1, then theta will be remapped to [-pi,pi] instead of [0,2pi] (this is useful for phase diagrams)" << endl;
  cout << "\talt. output timing flag:\t (default value = 0) If this is set to 1, then not only the t_n but every 100th value will be saved (this is useful for phase diagrams)" << endl;

  cout << endl;

  //read input parameters
  char* file_path = argv[1];
  double q = atof(argv[2]);
  double big_omega = atof(argv[3]);
  double f = atof(argv[4]);
  double max_t = atof(argv[5]);
  double delta_t = 2.0/big_omega * M_PI / atof(argv[6]);// atof(argv[6]);

  bool use_alternative_theta_mapping = false;
  bool use_alternative_output_timing = false;
  if (argc == 8)
  {
      use_alternative_theta_mapping = atoi(argv[7]) == 1;
  }
  else if (argc == 9)
  {
    use_alternative_theta_mapping = atoi(argv[7]) == 1;
    use_alternative_output_timing = atoi(argv[8]) == 1;
  }

  //print entered parameters for the user again
  cout << "Entered parameters:" << endl;
  cout << "\toutput file:\t" << file_path << endl;
  cout << "\tq:\t\t." << q << endl;
  cout << "\tOmega:\t\t" << big_omega << endl;
  cout << "\tf:\t" << f << endl;
  cout << "\tmaximum time:\t" << max_t << endl;
  if (use_alternative_theta_mapping)
    cout << "\tAlternative mapping of theta has been enabled." << endl;

  //create output file stream
  ofstream outputFile;
  outputFile.open(file_path);
  outputFile << "#t" << "\t" << "omega" << "\t" << "theta" << endl;

  //set the boundary conditions of the problem
  double theta_null = 0.0;
  double omega_null = big_omega;

  //create internal data for the RK4 integration implementation
  //set system of ODEs
  numerical::odeFunction* functions = new numerical::odeFunction[2];
  functions[0] = &f_d_omega;
  functions[1] = &f_d_theta;

  //create storage for the variables. will contain omega and theta according to the comment in the top
  double* rkValues = new double[2];
  //set start values
  rkValues[0] = omega_null;
  rkValues[1] = theta_null;

  //set the other parameters of the problem for the evaluation of the ODEs
  double* rkParams = new double[4];
  //rkParams[0] will contain the current time t
  rkParams[1] = q;
  rkParams[2] = f;
  rkParams[3] = big_omega;

  //gives the difference in t between to t_n values: print_interval_t = t_{n+1}-t_n
  double print_interval_t = 2.0/big_omega * M_PI;
  //is beeing used to store the time for which the values haven been saved the last time
  //use this value for the initialization, so that the first iteration will be used
  double last_printed_t = -print_interval_t;

  //this helper index is being used of the use_alternative_output_timing has been enabled. now not only the t_ns but every 100th value will be saved.
  //this is helpful to draw phasediagrams. set this to 100 to print the first value set.
  int last_printed_index = 100;
  for (double t=0; t<=max_t; t+=delta_t)
  {
    //to be sure that t = t_n we have to test if the interval between this time and the last time hast been 2*pi/Omega = print_interval_t
    if ((!use_alternative_output_timing  && t - last_printed_t >= print_interval_t) || (use_alternative_output_timing && last_printed_index == 100))
    {
      //set this index to zero to start counting again
      last_printed_index = 0;

      //memorize the current time
      last_printed_t = t;

      //print the data
      outputFile << t << "\t" << rkValues[0] << "\t" << rkValues[1] << endl;
    }

    //increase this hlper index
    last_printed_index++;

    //RK4 integration

    //set the time parameter
    rkParams[0] = t;

    //do one explicit RK4 integration step
    numerical::step_rk4_explicit(2, functions, delta_t, rkValues, rkParams);

    //remap the values of theta
    if (!use_alternative_theta_mapping)
    {
      //default remapping of theta into the interval [0,2*pi]
      while (rkValues[1] > 2*M_PI)
      {
        rkValues[1] -= 2*M_PI;
      }

      while (rkValues[1] < -0)
      {
        rkValues[1] += 2*M_PI;
      }
    }
    else
    {
      //alternative remapping of theta into the interval [-pi,pi]
      while (rkValues[1] > M_PI)
      {
        rkValues[1] -= 2*M_PI;
      }

      while (rkValues[1] < -M_PI)
      {
        rkValues[1] += 2*M_PI;
      }
    }
  }
}
