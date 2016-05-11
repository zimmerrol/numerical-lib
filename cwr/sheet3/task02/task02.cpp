#include "../../../lib/diff.h"
#include <math.h>
#include <cstdlib>
#include <cmath>

#include <iostream>
#include <fstream>

using namespace std;

double funcDiff(double x)
{
  return 5/pow(cosh(5*x),2)-(5*exp(-x*x+6*x-9))/4.0 * (-2*x+6);
}

double func(double x)
{
  return tanh(5*x) - (5*exp(-x*x+6*x-9))/4.0;
}

int main(int argc, char* argv[])
{
  ofstream outputFile;
  outputFile.open(argv[1]);

  double deltaX = atof(argv[2]);
  double resultSt = 0.0;
  double resultFw = 0.0;
  double resultbw = 0.0;

  for (double x=-10; x<10; x+=deltaX)
  {
    resultSt += pow(funcDiff(x) - diff::differentiateSterling(&func, x, deltaX), 2);
    resultFw += pow(funcDiff(x) - diff::differentiateNewtonGregroyForwards(&func, x, deltaX), 2);
    resultbw += pow(funcDiff(x) - diff::differentiateNewtonGregroyBackwards(&func, x, deltaX), 2);
    outputFile << x << "\t" << funcDiff(x) << "\t" << diff::differentiateNewtonGregroyForwards(&func,x,deltaX) << "\t" << diff::differentiateNewtonGregroyBackwards(&func,x,deltaX) << "\t" << diff::differentiateSterling(&func, x, deltaX) << "\n";
  }

  resultSt /= (20.0/deltaX);
  resultFw /= (20.0/deltaX);
  resultbw /= (20.0/deltaX);

  cout << "MSD(Sterling)=" << resultSt << "\n";
  cout << "MSD(NGF)=" << resultFw << "\n";
  cout << "MSD(NGB)=" << resultbw << "\n";
}
