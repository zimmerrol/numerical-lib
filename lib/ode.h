#include <cstddef>

namespace numerical
{
  typedef double (*odeFunction)(double* args, double* params);

  void stepLeapFrog(size_t dimensions, odeFunction* functions, double timeStepWidth, double* args, double* previousArgs, double* params);

  void stepEulerExplicit(size_t dimensions, odeFunction* functions, double timeStepWidth, double* args, double* params);
  //void stepEulerImplicit(size_t dimensions, odeFunction* functions, double timeStepWidth, double* args, double* params);

  void stepRKGeneral(size_t dimensions, odeFunction* functions, size_t order, double** rkParams, double timeStepWidth, double* args, double* params);

  void stepRK2Explicit(size_t dimensions, odeFunction* functions, double timeStepWidth, double* args, double* params);
  void stepRK3Explicit(size_t dimensions, odeFunction* functions, double timeStepWidth, double* args, double* params);
  void stepRK4Explicit(size_t dimensions, odeFunction* functions, double timeStepWidth, double* args, double* params);
}
