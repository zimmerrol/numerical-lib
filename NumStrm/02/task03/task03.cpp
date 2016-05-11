#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>

using namespace std;

void stepExplicit(double *x, double *y, double *z, double deltaT, double sigma, double r, double b)
{
  double origX = *x;
  double origY = *y;
  double origZ = *z;

  *x += deltaT * (sigma * (origY-origX));
  *y += deltaT * (-origX*origZ + r*origX - origY);
  *z += deltaT * (origX*origY - b*origZ);
}

int main(int argc, char* argv[])
{
  double xExp = 1;
  double yExp = 0;

  ofstream outputFile;
  outputFile.open(argv[1]);
  outputFile << fixed << setprecision(9);

  double deltaT = 0.001;//atof(argv[2]);
  double sigma = 10;//;atof(argv[3]);
  double r = 28;//atof(argv[4]);
  double b = 8.0/3.0;//atof(argv[5]);

  double xEx = 1,yEx = 1, zEx = 1;
  for(double t =0;t <= 100;t+=deltaT)
  {
    outputFile << t << "\t" << xEx << "\t" << yEx << "\t" << zEx << "\n";//"\t" << xIm << "\t" << yIm  << "\n";
    stepExplicit(&xEx, &yEx, &zEx, deltaT, sigma, r, b);
    //stepImplicit(&xIm, &yIm, deltaT);
  }
}
