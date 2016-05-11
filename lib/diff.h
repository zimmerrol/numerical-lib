typedef double(*Function1D)(double x, double* args);

namespace numerical
{
  double differentiateNewtonGregroyForwards(Function1D func, double x, double h, double* args);
  double differentiateNewtonGregroyBackwards(Function1D func, double x, double h, double* args);
  double differentiateSterling(Function1D func, double x, double h, double* args);
}
