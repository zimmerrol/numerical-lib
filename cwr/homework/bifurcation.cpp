#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "math.h"
#include <vector>
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


/*
for the creation of the bifurcation diagram the program has been called like this
  bifurcation data.dat 0.3333 0.75 1.35 1.45 0.0001 5000
*/

int main(int argc, char* argv[])
{
  //print information for the usage
  cout << "Use this to create a bifucation map for the drivven, damped pendulum." << endl;
  cout << "Expects the follwing set of parameters:" << endl;

  cout << "\toutput file:\tThe path to the file in which the results of the bifurcation will be stored." << endl;
  cout << "\tq:\t\tThe damping parameter." << endl;
  cout << "\tOmega:\t\tThe frequency of the driving force of the pendulum." << endl;
  cout << "\tstart f:\tThe start value for the iteration over the strength of the force f." << endl;
  cout << "\tend f:\t\tThe final value for the iteration over the strength of the force f." << endl;
  cout << "\tf step size:\tThe step size which is used in the iteration over the strength of the force f." << endl;
  cout << "\tmaximum time:\tThe maximum time for the integration (should be high to reach a fixed point)" << endl;
  cout << "\t[time step size is fixed to 2pi/(Omega * 5000)]" << endl << endl;

  //read input parameters
  char* file_path = argv[1];
  double q = atof(argv[2]);
  double big_omega = atof(argv[3]);
  double start_f = atof(argv[4]);
  double end_f = atof(argv[5]);
  double step_f = atof(argv[6]);
  double max_t = atof(argv[7]);

  //choose this delta_t to integrate until a valid t_n
  double delta_t = 2.0/big_omega * M_PI / 5000;

  //indicates the time after which the points will be used for the diagram. we need to calculate a lot of them before we can create the map to be sure,
  //that we have archived a fixed osicllation period. Usage of the last 100 t_n points is sufficient for the detection of the first period doubblings.
  double start_print_t = max_t - 2.0/big_omega * M_PI * 100;

  //print entered parameters for the user again
  cout << "Entered parameters:" << endl;
  cout << "\toutput file:\t" << file_path << endl;
  cout << "\tq:\t\t." << q << endl;
  cout << "\tOmega:\t\t" << big_omega << endl;
  cout << "\tstart f:\t" << start_f << endl;
  cout << "\tend f:\t\t." << end_f << endl;
  cout << "\tf step size:\t" << step_f << endl;
  cout << "\tmaximum time:\t" << max_t << endl << endl;

  //create output file stream
  ofstream outputFile;
  outputFile.open(file_path);

  //print the output data format
  outputFile << "#f" << "\t" << "theta" << endl;

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

  //set the other parameters of the problem for the evaluation of the ODEs
  double* rkParams = new double[4];

  //gives the difference in t between to t_n values: print_interval_t = t_{n+1}-t_n
  double print_interval_t = 2.0/big_omega * M_PI;

  //indicates the current period which has been detected below
  int detected_period = 1;
  for (double f=start_f; f<=end_f; f+=step_f)
  {
    //set the needed parameters
    //rkParams[0] will contain the current time t
    rkParams[1] = q;
    rkParams[2] = f;
    rkParams[3] = big_omega;

    //set start values
    rkValues[0] = omega_null;
    rkValues[1] = theta_null;

    //is beeing used to store the time for which the values haven been saved the last time
    //use this value for the initialization, so that the first iteration will be used
    double last_printed_t = -print_interval_t;

    //small paramter to compare double values
    //x is considered to be equal to y if |x-y| < epsilon is satisfied
    double epsilon = 1e-6;

    //stores all fixed points for this f iteration
    vector<double> calculated_values;

    for (double t=0; t<=max_t; t+=delta_t)
    {
      //only print/use the last 100 t_n for the map/diagram (see above)
      if (t > start_print_t)
      {
        //to be sure that t = t_n we have to test if the interval between this time and the last time hast been 2*pi/Omega = print_interval_t
        if (abs(t - last_printed_t + epsilon) >= print_interval_t)
        {
          //memorize the current time
          last_printed_t = t;

          //add the new value of theta (the fixed point) to the stack/vector
          calculated_values.push_back(rkValues[1]);
        }
      }

      //RK4 integration starts here
      //set the time parameter
      rkParams[0] = t;

      //do one explicit RK4 integration step
      numerical::step_rk4_explicit(2, functions, delta_t, rkValues, rkParams);

      //remap the values of theta into the interval [0,2*pi]
      while (rkValues[1] > 2*M_PI)
      {
        rkValues[1] -= 2*M_PI;
      }

      while (rkValues[1] < -0)
      {
        rkValues[1] += 2*M_PI;
      }
    }

    //detect the current period to detect the perdio doubling
    //loop over all possible period values for the current amount of recored values for this f iteration
    for (int period=1; period < (int)calculated_values.size()-1; period++)
    {
      bool found_period = true;
      //loop over all items and check if their period is "period" <=> calculated_values[i+period] is equal to calculated_values[i]
      for (int index=0; index < (int)calculated_values.size()-1-period; index++)
      {
        if (abs(calculated_values.at(index) - calculated_values.at(index+period)) > epsilon)
        {
          //this point is not periodically with the period "period"
          //mark this issue and break the loop to continue to try further period values
          found_period = false;
          break;
        }
      }

      //not all items seem to be periodically with the period "period" => try to use the next bigger period
      if (!found_period)
      {
        continue;
      }

      //all items seem to be periodically with the period "period" and this period has been altered compared to the previously detected period "detected_period"
      if (period != detected_period)
      {
        //print this information to the user and memorize the new period
        cout << "The period has been doubled (f = " << f <<  ")!\t" << detected_period << " => " << period << endl;
        detected_period = period;
      }

      break;
    }

    //detect multiple items now to keep the output file small
    for (size_t i=0; i<calculated_values.size(); i++)
    {
      //loop over all other items to detect multiples of item i
      for (size_t j=i+1; j<calculated_values.size(); j++)
      {
        //compare if the values seem to be equal (accoring to definition above)
        if (abs(calculated_values.at(i) - calculated_values.at(j)) < epsilon)
        {
          calculated_values.erase(calculated_values.begin()+j);
          j--;
        }
      }
    }

    //save the calculated theta values now
    for (size_t i=0; i<calculated_values.size(); i++)
    {
      outputFile << f << "\t" << calculated_values.at(i) << endl;
    }
  }
}
