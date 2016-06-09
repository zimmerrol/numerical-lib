#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "../../../lib/square_matrix.h"
#include "../../../lib/linear_system.h"

using namespace std;
using namespace numerical;

double analytic_solution(double x, double D, double v)
{
  return (exp(v/D*x)-1)/(exp(v/D)-1);
}

//returns the error after one step
void step(size_t n, double* T, double* a, double*b, double* c, double*d)
{
  c[1] /= b[1];
   d[1] /= b[1];

   for(size_t i = 2; i < n; i++)
   {
     c[i] /= (b[i] - a[i] * c[i-1] );
     d[i] = (d[i] - a[i] * d[i-1])/(b[i] - a[i] * c[i-1]);
   }

   for(size_t i = n-1; i>= 1 ; i--)
   {
     T[i] = d[i] - c[i]*T[i+1];
   }
}

void integrate_exim(double* grid, double &kappa, double &v,double &dx, double &dt, size_t N)
{
  double* a = new double[N];
  double* b = new double[N];
  double* c = new double[N];
  double* d = new double[N];

  double alpha = 1/dt + 2*kappa/dx/dx;
  double beta = -kappa/dx/dx;
  //stay with the convention to label a, b, c starting from 1 to avoid confusion
  // with subscribts of grid points

  for(size_t i = 2; i < N-1; i++)
  {
    b[i] = alpha;
    a[i] = beta; c[i] = beta;
    d[i] = v/2/dx*(grid[i-1] - grid[i+1]) + grid[i]/dt;
  }

  b[1] = alpha;
  d[1] = v/2/dx*(grid[0] - grid[2]) + grid[1]/dt ;

  d[N-1] = v/2/dx*(grid[N-2] - grid[N]) + kappa*grid[N]/dx/dx + grid[N-1]/dt;

  c[1] = beta; a[N-1] = beta;

  step(N-1, a, b, c, d, grid);
  delete a; delete b; delete c; delete d;
}

int main(int argc, char* argv[])
{
  cout << "Expects:" << endl;
  cout << "\t" << "outputPath\tdimension_x\tmax_time\tdelta_t\tv\tD" << endl;
  ofstream outputFile;
	outputFile.open(argv[1]);
	outputFile << fixed << setprecision(5);

  size_t dimension_x = (atoi)(argv[2]);
  dimension_x++;

  double delta_x = 1.0/(dimension_x-1);

  int max_time = (int)atof(argv[3]);
  double delta_t = atof(argv[4]);

  double v = atof(argv[5]);
  double D = atof(argv[6]);

  cout << "delta_t="<<delta_t << "\tmax_time=" << max_time << "\tdelta_x=" << delta_x << "\tv=" << v << "\tD=" << D << endl;

  double* values = new double[dimension_x];
  for (size_t i=0; i<dimension_x-1; i++)
  {
    values[i] = i*delta_x;
  }

  for(double t=0; t<max_time; t+=delta_t)
  {
    integrate_exim(values, D,v, delta_x, delta_t, dimension_x+1);
  }

  double diff = 0.0;
  for (size_t i=0; i<dimension_x; i++)
  {
    outputFile << i*delta_x << "\t" << values[i] << "\t" << analytic_solution(i*delta_x, D, v) << endl;
    diff += pow(analytic_solution(i*delta_x, D, v) - values[i],2);
  }

  cout << "MSE = " << sqrt(diff) << endl;
  cout << flush;
}
