#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>

using namespace std;

const double k = 1;//;atof(argv[3]);
const double M = 1;//atof(argv[4]);

void stepExplicit(double *x, double *v, double deltaT)
{
  *v += -k/M* *x;
  *x += deltaT * *v;
}


void stepRK2(double *x, double *v, double deltaT)
{
  double k2v = *v + 0.5*deltaT * (-k/M* *x);
  double k2x = *x + 0.5*deltaT * *v;
  *v += deltaT * k2v;
  *x += deltaT * k2x;
}

void stepRK4(double *x, double *v, double deltaT)
{
  //double k1 =
  //*x += deltaT * ();
}

int main(int argc, char* argv[])
{
  double xExp = 1;
  double yExp = 0;

  ofstream outputFile;
  outputFile.open(argv[1]);
  outputFile << fixed << setprecision(9);

  double deltaT = 0.001;//atof(argv[2]);


  double xEx = 1, xRK2 = 1, xRK4 = 1;
  double vEx = 0, vRK2 = 0, vRK4 = 0;

  for(double t =0;t <= 5;t+=deltaT)
  {
    outputFile << t << "\t" << xEx << "\t" << vEx << "\t" << xRK2 << "\t" << vRK2 << "\n";//"\t" << xIm << "\t" << yIm  << "\n";
    stepExplicit(&xEx, &vEx, deltaT);
    stepRK2(&xRK2, &vRK2, deltaT);
    //stepImplicit(&xIm, &yIm, deltaT);
  }
}
