#include <cstddef>

namespace numerical
{
  typedef double (*odeFunction)(double* args, double* params);

  void step_leap_frog(size_t dimensions, odeFunction* functions, double timeStepWidth, double* args, double* previousArgs, double* params);

  void step_euler_explicit(size_t dimensions, odeFunction* functions, double timeStepWidth, double* args, double* params);
  //void stepEulerImplicit(size_t dimensions, odeFunction* functions, double timeStepWidth, double* args, double* params);

  void step_rk_general(size_t dimensions, odeFunction* functions, size_t order, double** rkParams, double timeStepWidth, double* args, double* params);

  void step_rk2_explicit(size_t dimensions, odeFunction* functions, double timeStepWidth, double* args, double* params);
  void step_rk3_explicit(size_t dimensions, odeFunction* functions, double timeStepWidth, double* args, double* params);
  void step_rk4_explicit(size_t dimensions, odeFunction* functions, double timeStepWidth, double* args, double* params);
}
