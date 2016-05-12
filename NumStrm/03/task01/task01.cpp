#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <complex>

using namespace std;

complex<double> dx(complex<double> x, complex<double> lambda)
{
  return lambda * x;
}

int main(int argc, char* argv[])
{
  ofstream euler_ofstream;
  euler_ofstream.open(argv[1]);
  euler_ofstream << fixed << setprecision(9);

  ofstream ab2_ofstream;
  ab2_ofstream.open(argv[2]);
  ab2_ofstream << fixed << setprecision(9);

  ofstream ab3_ofstream;
  ab3_ofstream.open(argv[3]);
  ab3_ofstream << fixed << setprecision(9);

  double PI = acos(-1);


  //euler
  for (double theta = 0.000;theta < 2*PI;theta += 0.001)
  {
      complex<double> h_times_lambda = polar(1.0,theta) - 1.0;
      euler_ofstream << real(h_times_lambda) << "\t" << imag(h_times_lambda) << "\n";
  }

  //ab2
  for (double theta = 0.000;theta < 2*PI;theta += 0.001)
  {
      complex<double> epsilon = polar(1.0,theta);
      complex<double> h_times_lambda = 2.0*(epsilon*epsilon-epsilon) / (3.0*epsilon - 1.0);
      ab2_ofstream << real(h_times_lambda) << "\t" << imag(h_times_lambda) << "\n";
  }

  //ab3
  for (double theta = 0.000;theta < 2*PI;theta += 0.001)
  {
      complex<double> epsilon = polar(1.0,theta);
      complex<double> h_times_lambda = 12.0*(epsilon*epsilon*epsilon-epsilon*epsilon) / (23.0*epsilon*epsilon - 16.0 * epsilon + 5.0);
      ab3_ofstream << real(h_times_lambda) << "\t" << imag(h_times_lambda) << "\n";
  }
}
