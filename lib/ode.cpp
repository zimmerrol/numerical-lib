#include "ode.h"
#include <cstddef>
#include <iostream>
#include <fstream>

namespace numerical{
  void step_euler_explicit(std::size_t dimensions, odeFunction* functions, double timeStepWidth, double* args, double* params)
  {
    //create copy of all submitted arguments
    double* origArgs = new double[dimensions];
    for(size_t i=0;i<dimensions;i++)
    {
      origArgs[i] = args[i];
    }

    for(size_t i=0;i<dimensions;i++)
    {
      args[i] += timeStepWidth * (*(functions[i]))(origArgs, params);
    }

    delete origArgs;
  }

  void step_leap_frog(size_t dimensions, odeFunction* functions, double timeStepWidth, double* args, double* previousArgs, double* params)
  {
    //create copy of all submitted arguments
    double* origArgs = new double[dimensions];
    for(size_t i=0;i<dimensions;i++)
    {
      origArgs[i] = args[i];
    }

    for(size_t i=0;i<dimensions;i++)
    {
      double prevI = previousArgs[i];
      previousArgs[i] = args[i];
      //x_i+1 = x_i-1 + 2* deltaT * x'_i
      args[i] = prevI + 2 * timeStepWidth * (*(functions[i]))(origArgs, params);
    }

    delete origArgs;
  }

  //first item in rkParams* is the the factor for the final summation
  //the other items are the factors for influence of the previous samples for the new relative position calculation
  void step_rk_general(size_t dimensions, odeFunction* functions, size_t order, double** rkParams, double timeStepWidth, double* args, double* params)
  {
    double** samplingPoints = new double*[order];

    for(size_t n = 0;n<order;n++)
    {
      //create copy of all submitted arguments
      samplingPoints[n] = new double[dimensions];
      for(size_t i=0;i<dimensions;i++)
      {
        samplingPoints[n][i] = args[i];

        double sum = 0.0;
        for (size_t prevIndex = 0;prevIndex < n; prevIndex++)
        {
          sum += rkParams[n][prevIndex+1] * samplingPoints[prevIndex][i];
        }
        samplingPoints[n][i] += timeStepWidth * sum;
      }

      //create backup
      double* origSamplingPoints = new double[dimensions];
      for(size_t i=0;i<dimensions;i++)
      {
        origSamplingPoints[i] = samplingPoints[n][i];
      }

      //calc new step
      for(size_t eq = 0;eq<dimensions;eq++)
      {
        samplingPoints[n][eq] += timeStepWidth * (*(functions[eq]))(origSamplingPoints, params);
      }
    }


    for(size_t i=0;i<dimensions;i++)
    {
      double value = 0.0;
      for(size_t k=0;k<order-1;k++)
      {
        value += rkParams[k][0] * (*(functions[k]))(samplingPoints[i], params);
      }

      args[i] += value;
    }

  }

  void step_rk2_explicit(size_t dimensions, odeFunction* functions, double timeStepWidth, double* args, double* params)
  {
    double* k1 = new double[dimensions];

    for (size_t i = 0;i < dimensions; i++)
    {

      k1[i] = args[i] + 0.5*timeStepWidth * (*(functions[i]))(args, params);
    }

    double* k2 = new double[dimensions];
    for (size_t i = 0;i < dimensions; i++)
    {
      k2[i] = args[i] + timeStepWidth * (*(functions[i]))(k1, params);
    }

    for (size_t i = 0;i < dimensions; i++)
    {
      args[i] = k2[i];
    }

    delete k1;
    delete k2;
  }

  void step_rk3_explicit(size_t dimensions, odeFunction* functions, double timeStepWidth, double* args, double* params)
  {
    double* callingArguments = new double[dimensions];
    for (size_t i = 0;i < dimensions; i++)
    {
        callingArguments[i] = args[i];
    }

    double* k1 = new double[dimensions];
    for (size_t i = 0;i < dimensions; i++)
    {
      k1[i] = (*(functions[i]))(callingArguments, params);
    }

    for (size_t i = 0;i < dimensions; i++)
    {
        callingArguments[i] = args[i] + timeStepWidth/2.0 * k1[i];
    }

    double* k2 = new double[dimensions];
    for (size_t i = 0;i < dimensions; i++)
    {
      k2[i] = (*(functions[i]))(callingArguments, params);
    }

    for (size_t i = 0;i < dimensions; i++)
    {
        callingArguments[i] = args[i] + timeStepWidth * (-k1[i] + 2*k2[i]);
    }

    double* k3 = new double[dimensions];
    for (size_t i = 0;i < dimensions; i++)
    {
      k3[i] = (*(functions[i]))(callingArguments, params);
    }

    for (size_t i = 0;i < dimensions; i++)
    {
      args[i] += timeStepWidth * (1.0/6.0 * k1[i] + 4.0/6.0 * k2[i] + 1.0/6.0 * k3[i]);
    }

    delete callingArguments;
    delete k1;
    delete k2;
    delete k3;
  }

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
