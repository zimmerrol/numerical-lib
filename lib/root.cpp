#include "root.h"
#include <vector>
#include <tuple>
#include <cmath>
#include <iostream>

using namespace std;

namespace numerical
{
  int countRoots(Function1D func, double start, double end, double deltaX, std::vector<std::tuple<double, double>>* rootPositions)
  {
    int result = 0;

    bool isGtn = (*func)(start) > 0;
    for (double x = start; x <= end; x+=deltaX)
    {
      if (((*func)(x) > 0) != isGtn)
      {
        isGtn = !isGtn;
        result++;
        (*rootPositions).push_back(make_tuple(x-deltaX, x));
      }
    }

    return result;
  }

  double findRootRegulaFalsi(Function1D func, double start, double end, double epsilon)
  {
    if (abs(start - end) <= epsilon)
    {

      return (start+end)/2;
    }

    double x = start + (end-start)/((*func)(end) - (*func)(start))*(*func)(start);

    if ((*func)(start) * (*func)(x) < 0)
    {
      return findRootRegulaFalsi(func, start, x, epsilon);
    }
    else
    {
      return findRootRegulaFalsi(func, x, end, epsilon);
    }
  }

  double findRootBisection(Function1D func, double start, double end, double epsilon)
  {
    double x = (start+end)/2;
    if (abs(start - end) < epsilon)
    {
      return x;
    }

    if ((*func)(start) * (*func)(x) < 0)
    {
      return findRootBisection(func, start, x, epsilon);
    }
    else
    {
      return findRootBisection(func, x, end, epsilon);
    }
  }

  //start = x_n-1
  //currentX = x_n
  double findRootSecantMethod(Function1D func, double previousX, double currentX, double epsilon)
  {

    double x = (previousX-currentX);
      cout << x << " " << flush;
    if (x < epsilon)
    {
      return (previousX+currentX)/2;
    }



    double xNew = previousX - (currentX - previousX)/((*func)(currentX) - (*func)(previousX))*(*func)(previousX);

    return findRootSecantMethod(func, currentX, xNew, epsilon);
  }

  double findRootNewtonRaphson(Function1D func, Function1D funcDiff, double previousX, double currentX, double epsilon)
  {
    double x = (previousX-currentX);
    if (x < epsilon)
    {
      return x/2;
    }

    double xNew = currentX + (*func)(currentX) / (*funcDiff)(currentX);

    return findRootNewtonRaphson(func, funcDiff, currentX, xNew, epsilon);
  }
}
