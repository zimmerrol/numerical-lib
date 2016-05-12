#include "../../../lib/diff.h"
#include <math.h>
#include <cstdlib>
#include <cmath>

#include <iostream>
#include <fstream>

using namespace std;



double step(double x, double mu)
{
  return 4*mu*x*(1-x);
}

double func(double x, double* args)
{
  double mu = args[0];
  return step(x,mu);
}

int main(int argc, char* argv[])
{
  ofstream outputFile;
  outputFile.open(argv[1]);

  double deltaMu = atof(argv[2]);
  double deltaX = atof(argv[3]);

  double* args = new double[1];

  for (double mu=0; mu<1; mu+=deltaMu)
  {
    double lambda = 0.0;
    double x = 0.8;
    for (int n=0; n<1000; n++)
    {
      //calc x
      x = step(x, mu);
      args[0] = mu;
      lambda += log(abs(numerical::differentiateSterling(&func, x, deltaX, args)));
    }

    lambda /= (1000);

    outputFile << mu << "\t" << lambda << "\n";
  }
}
