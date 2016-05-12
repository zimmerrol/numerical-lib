#include "ode.h"
#include <cstddef>
#include <iostream>
#include <fstream>

namespace numerical{


  void stepEulerExplicit(std::size_t dimensions, odeFunction* functions, double timeStepWidth, double* args, double* params)
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
  }

  void stepLeapFrog(size_t dimensions, odeFunction* functions, double timeStepWidth, double* args, double* previousArgs, double* params)
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
  void stepRKGeneral(size_t dimensions, odeFunction* functions, size_t order, double** rkParams, double timeStepWidth, double* args, double* params)
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

  void stepRK2Explicit(size_t dimensions, odeFunction* functions, double timeStepWidth, double* args, double* params)
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
  }

  void stepRK3Explicit(size_t dimensions, odeFunction* functions, double timeStepWidth, double* args, double* params)
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
  }

  void stepRK4Explicit(size_t dimensions, odeFunction* functions, double timeStepWidth, double* args, double* params)
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
        callingArguments[i] = args[i] + timeStepWidth/2.0 * k2[i];
    }

    double* k3 = new double[dimensions];
    for (size_t i = 0;i < dimensions; i++)
    {
      k3[i] = (*(functions[i]))(callingArguments, params);
    }

    for (size_t i = 0;i < dimensions; i++)
    {
        callingArguments[i] = args[i] + timeStepWidth * k3[i];
    }

    double* k4 = new double[dimensions];
    for (size_t i = 0;i < dimensions; i++)
    {
      k4[i] = (*(functions[i]))(callingArguments, params);
    }

    for (size_t i = 0;i < dimensions; i++)
    {
      args[i] += timeStepWidth/6.0 * (1.0 * k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + 1.0 * k4[i]);
    }
  }
}
