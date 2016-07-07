#include <cstddef>

namespace numerical
{
  //definition for functions which calculate the derivative of one variable in the RK4 integration method (see body for more information)
  typedef double (*odeFunction)(double* args, double* params);

  //default RK4 integration
  void step_rk4_explicit(size_t dimensions, odeFunction* functions, double timeStepWidth, double* args, double* params);
}
