#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using namespace std;

const double g = -1.62;

void step(double* x, double* y, double* vx, double* vy, double deltaT, double c )
{
  *vy +=  deltaT * g - c* *vy * deltaT;
  *vx -= c* *vx * deltaT;
  *x += deltaT * *vx;
  *y += deltaT * *vy;
}

int main(int argc, char* argv[])
{
  double v0 = atof(argv[1]);
  double angle = atof(argv[2]);
  double deltaT = atof(argv[3]);

  double vy = v0*sin(angle);
  double vx = v0*cos(angle);
  double x = atof(argv[4]);
  double y = atof(argv[5]);
  double c = atof(argv[6]);

  ofstream outputFile;
  outputFile.open(argv[7]);

  double t = 0;
  do
  {
    outputFile << t << "\t"  << x << "\t" << y << "\t" << vx << "\t" << vy << "\n";
    step(&x, &y, &vx, &vy, deltaT, c);
    t += deltaT;
  }
  while ((vy < 0 && y > 0) || vy > 0);

  outputFile.flush();
  outputFile.close();
}
