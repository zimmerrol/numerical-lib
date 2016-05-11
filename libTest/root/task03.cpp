#include <math.h>
#include <cstdlib>
#include <cmath>
#include "../../lib/root.h"
#include <iostream>
#include <fstream>
#include <tuple>
#include <iomanip>

using namespace std;
using namespace numerical;

double func(double x)
{
  return tanh(5*x) - (5*exp(-x*x+6*x-9))/4.0;
  //return sin(x);
}

int main(int argc, char* argv[])
{
  //ofstream outputFile;
  //outputFile.open(argv[1]);

  double deltaX = atof(argv[1]);
  double epsilon = atof(argv[2]);

  std::vector<std::tuple<double, double>> result;

cout << std::setprecision(15);

  int res = countRoots(&func, -10, 10, deltaX, &result);

  for (int i=0;i<res;i++)
  {
    tuple<double, double> item = result.at(i);
    cout << "start=" << std::get<0>(item) << "\tend=" << std::get<1>(item);
    cout << flush;

    double exactRf = findRootRegulaFalsi(&func, std::get<0>(item),std::get<1>(item), epsilon);
    double exactBi = findRootBisection(&func, std::get<0>(item),std::get<1>(item), epsilon);
    double exactSe = findRootSecantMethod(&func, std::get<0>(item),std::get<1>(item), epsilon);

    cout <<  "\texactRf=" << exactRf << "\texactBi=" << exactBi  << "\texactSe=" << exactSe << "\n";


  }

cout << "#roots= " << res << "\n";

  for (double x=-10; x<10; x+=deltaX)
  {

  }

}
