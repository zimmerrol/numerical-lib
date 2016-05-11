#include "diff.h"

using namespace std;

double numerical::differentiateNewtonGregroyForwards(Function1D func, double x, double h, double* args)
{
  return ((*func)(x+h, args)-(*func)(x, args))/h;
}

double numerical::differentiateNewtonGregroyBackwards(Function1D func, double x, double h, double* args)
{
  return ((*func)(x, args)-(*func)(x-h, args))/h;
}

double numerical::differentiateSterling(Function1D func, double x, double h, double* args)
{
  return ((*func)(x+h, args)-(*func)(x-h, args))/(2*h);
}
