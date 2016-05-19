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

double func(double x, double* args)
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
    resultSt += pow(funcDiff(x) - numerical::differentiate_centered_difference(&func, x, deltaX, NULL), 2);
    resultFw += pow(funcDiff(x) - numerical::differentiate_newton_gregroy_forwards(&func, x, deltaX, NULL), 2);
    resultbw += pow(funcDiff(x) - numerical::differentiate_newton_gregroy_backwards(&func, x, deltaX, NULL), 2);
    outputFile << x << "\t" << funcDiff(x) << "\t" << numerical::differentiate_newton_gregroy_forwards(&func,x,deltaX, NULL) << "\t" << numerical::differentiate_newton_gregroy_backwards(&func,x,deltaX, NULL) << "\t" << numerical::differentiate_centered_difference(&func, x, deltaX, NULL) << "\n";
  }

  resultSt /= (20.0/deltaX);
  resultFw /= (20.0/deltaX);
  resultbw /= (20.0/deltaX);

  cout << "MSD(Sterling)=" << resultSt << "\n";
  cout << "MSD(NGF)=" << resultFw << "\n";
  cout << "MSD(NGB)=" << resultbw << "\n";
}
