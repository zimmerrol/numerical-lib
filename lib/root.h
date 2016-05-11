#include <vector>
#include <tuple>

typedef double(*Function1D)(double x);

namespace numerical
{
  int countRoots(Function1D func, double start, double end, double stepWidth, std::vector<std::tuple<double, double>>* rootPositions);
  double findRootRegulaFalsi(Function1D func, double start, double end, double epsilon);
  double findRootBisection(Function1D func, double start, double end, double epsilon);
  double findRootSecantMethod(Function1D func, double start, double end, double epsilon);
  double findRootNewtonRaphson(Function1D func, double start, double end, double epsilon);
}
