#include <math.h>
#include <cstdlib>
#include <cmath>
#include "root.h"
#include <iostream>
#include <fstream>
#include <tuple>

using namespace std;

double func(double x)
{
  //return tanh(5*x) - (5*exp(-x*x+6*x-9))/4.0;
  return sin(x);
}

int main(int argc, char* argv[])
{
  //ofstream outputFile;
  //outputFile.open(argv[1]);

  double deltaX = atof(argv[1]);
  double epsilon = atof(argv[2]);

  std::vector<std::tuple<double, double>> result;

  int res = root::countRoots(&func, -10, 10, deltaX, &result);

  for (int i=0;i<res;i++)
  {
    tuple<double, double> item = result.at(i);
    double exact = root::findRootRegulaFalsi(&func, std::get<0>(item),std::get<1>(item), epsilon);
    cout << "start=" << std::get<0>(item) << "\tend=" << std::get<1>(item) << "\texact=" << exact << "\n";


  }

cout << res << "\n";

  for (double x=-10; x<10; x+=deltaX)
  {

  }

}
