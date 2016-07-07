#include "ode.h"
#include <cstddef>
#include <iostream>
#include <fstream>

namespace numerical{
  //calculates one RK4 ecplicit time step (default implementation of the algorithm as described in the pdf)
  //dimensions = amount if variables which have to be integrated
  //functions = array of odeFunctions (= pointer to a function which calculates the derivative of one variable)
  //timeStepWidth = the size of the time step delta_t
  //args = a list of variables for the functions which will ALL be integrated.
  //params = a list of parameters for the functions which will NOT be integrated.
  void step_rk4_explicit(size_t dimensions, odeFunction* functions, double timeStepWidth, double* args, double* params)
  {
    //create a working version of the argument list
    double* callingArguments = new double[dimensions];
    for (size_t i = 0;i < dimensions; i++)
    {
        callingArguments[i] = args[i];
    }

    //calculate the value k1 (according to the pdf-file)
    double* k1 = new double[dimensions];
    for (size_t i = 0;i < dimensions; i++)
    {
      k1[i] = (*(functions[i]))(callingArguments, params);
    }

    //use the variables of the step K1 for the calulcation of k2
    for (size_t i = 0;i < dimensions; i++)
    {
        callingArguments[i] = args[i] + timeStepWidth/2.0 * k1[i];
    }

    //calculate the value k2 (according to the pdf-file)
    double* k2 = new double[dimensions];
    for (size_t i = 0;i < dimensions; i++)
    {
      k2[i] = (*(functions[i]))(callingArguments, params);
    }

    //use the variables of the step K2 for the calulcation of k3
    for (size_t i = 0;i < dimensions; i++)
    {
        callingArguments[i] = args[i] + timeStepWidth/2.0 * k2[i];
    }

    //calculate the value k3 (according to the pdf-file)
    double* k3 = new double[dimensions];
    for (size_t i = 0;i < dimensions; i++)
    {
      k3[i] = (*(functions[i]))(callingArguments, params);
    }

    //use the variables of the step K3 for the calulcation of k4
    for (size_t i = 0;i < dimensions; i++)
    {
        callingArguments[i] = args[i] + timeStepWidth * k3[i];
    }

    //calculate the value k4 (according to the pdf-file)
    double* k4 = new double[dimensions];
    for (size_t i = 0;i < dimensions; i++)
    {
      k4[i] = (*(functions[i]))(callingArguments, params);
    }

    //add all the calculated k-values and integrate the variables according to the default RK4 scheme (see pdf file)
    for (size_t i = 0;i < dimensions; i++)
    {
      args[i] += timeStepWidth/6.0 * (1.0 * k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + 1.0 * k4[i]);
    }

    //tidy up & delete all created arrays
    delete callingArguments;
    delete k1;
    delete k2;
    delete k3;
    delete k4;
  }
}
