#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "../../../lib/square_matrix.h"
#include "../../../lib/linear_system.h"

using namespace std;
using namespace numerical;

double analytic_solution(double x, double D, double v)
{
  return (exp(v/D*x)-1)/(exp(v/D)-1);
}

void integrate(vector<double> &values, square_matrix &coeff_matrix, double D, double v, double delta_x, double delta_t)
{
  size_t n = coeff_matrix.get_n();
  double T_left = 0.0;
  double T_right = 1.0;

  vector<double> x(n);
  vector<double> b(n);

//b starts at T_1 and yields till T_N-1
  b.at(0) = T_left*D*delta_t + v/2.0*delta_t*delta_x*(T_left-values.at(2)) + pow(delta_x,2)*values.at(1);
  b.at(n-1) = T_right*D*delta_t + v/2.0*delta_t*delta_x*(values.at(n) - T_right) + pow(delta_x, 2)*values.at(n+1);

  for (size_t i=1; i<n-1; i++)
  {
    b.at(i) = v/2.0*delta_t*delta_x*(values.at(i)-values.at(i+2)) + pow(delta_x,2)*values.at(i+1);
  }

  linear_system_solve_tridiagonal(coeff_matrix, &x[0], &b[0]);

  values.at(0) = T_left;
  values.at(n+1) = T_right;
  for (size_t i=0; i<n; i++)
  {
    values.at(i+1) = x.at(i);
  }
}

int main(int argc, char* argv[])
{
  cout << "Expects:" << endl;
  cout << "\t" << "outputPath\tdimension_x\tmax_time\tdelta_t\tv\tD" << endl;
  ofstream outputFile;
	outputFile.open("res/out.dat");
	outputFile << fixed << setprecision(5);

  size_t dimension_x = 100;

  double delta_x = 1.0/(dimension_x);

  int max_time = 20;
  double delta_t = 3.85e-2;

  double v = 1;
  double D = 0.01;

  cout << "delta_t="<<delta_t << "\tmax_time=" << max_time << "\tdelta_x=" << delta_x << "\tv=" << v << "\tD=" << D << endl;

  vector<double> values(dimension_x);
  for (size_t i=0; i<dimension_x; i++)
  {
    values.at(i) = i*delta_x;
  }

  double alpha = delta_x*delta_x + 2*D*delta_t;
  double beta = D*delta_t;

  square_matrix coeff_matrix = square_matrix(dimension_x-2,0);
  coeff_matrix.set_value(0,0,alpha);
  coeff_matrix.set_value(0,1,-beta);

  for (size_t i=1; i<dimension_x-3; i++)
  {
    coeff_matrix.set_value(i,i,alpha);
    coeff_matrix.set_value(i,i-1,-beta);
    coeff_matrix.set_value(i,i+1,-beta);
  }

  coeff_matrix.set_value(dimension_x-3,dimension_x-4,-beta);
  coeff_matrix.set_value(dimension_x-3,dimension_x-3,alpha);


  for(double t=0; t<max_time; t+=delta_t)
  {
    integrate(values, coeff_matrix, D, v, delta_x, delta_t);
  }

  double diff = 0.0;
  for (size_t i=0; i<dimension_x; i++)
  {
    outputFile << i*delta_x << "\t" << values[i] << "\t" << analytic_solution(i*delta_x, D, v) << endl;
    diff += pow(analytic_solution(i*delta_x, D, v) - values[i],2);
  }

  cout << "MSE = " << sqrt(diff)/dimension_x << endl;
  cout << flush;

  system("gnuplot plot.plt");
}
