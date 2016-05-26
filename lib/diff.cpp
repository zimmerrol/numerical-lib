#include "diff.h"

using namespace std;

double numerical::differentiate_newton_gregroy_forwards(Function1D func, double x, double h, double* args)
{
  return ((*func)(x+h, args)-(*func)(x, args))/h;
}

double numerical::differentiate_newton_gregroy_backwards(Function1D func, double x, double h, double* args)
{
  return ((*func)(x, args)-(*func)(x-h, args))/h;
}

double numerical::differentiate_centered_difference(Function1D func, double x, double h, double* args)
{
  return ((*func)(x+h, args)-(*func)(x-h, args))/(2*h);
}
