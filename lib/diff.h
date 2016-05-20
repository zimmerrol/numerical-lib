typedef double(*Function1D)(double x, double* args);

namespace numerical
{
  double differentiate_newton_gregroy_forwards(Function1D func, double x, double h, double* args);
  double differentiate_newton_gregroy_backwards(Function1D func, double x, double h, double* args);
  double differentiate_centered_difference(Function1D func, double x, double h, double* args);
}
