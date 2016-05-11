#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>

using namespace std;

void stepExplicit(double *x, double *y, double deltaT)
{
  double origX = *x;
  *x += deltaT *( (998 * *x) + 1998 * *y);
  *y += deltaT * ((-999 * origX) - 1999 * *y);
}

void stepImplicit(double *x, double *y, double h)
{
  double origX = *x;
  double oneOverDet = 1.0/(1000*h*h + 1001*h+1);

  *x = oneOverDet * ((1999*h+1) * *x + (1998*h) * *y);
  *y = oneOverDet * ((-999*h)* origX + (1-998*h) * *y);
}

int main(int argc, char* argv[])
{
  double xExp = 1;
  double yExp = 0;

  ofstream outputFile;
  outputFile.open(argv[1]);
  outputFile << fixed << setprecision(9);

  double deltaT = atof(argv[2]);
  double xEx = 1,yEx = 0;
  double xIm = 1,yIm = 0;

  for(double t =0;t <= 1;t+=deltaT)
  {
    outputFile << t << "\t" << xEx << "\t" << yEx << "\t" << xIm << "\t" << yIm  << "\n";
    stepExplicit(&xEx, &yEx, deltaT);
    stepImplicit(&xIm, &yIm, deltaT);
  }

}
