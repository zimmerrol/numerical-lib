#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

#include "fft_fn.h"

using namespace std;

int main(int argc, char* argv[])
{
  cout << std::setprecision(4);

  size_t n = pow(2,10);

  double delta_x = 2*M_PI/n;
  float* data_ft = new float[n+1];

  for (size_t i=0; i<n; i++)
  {
    data_ft[i] = sin(i*delta_x);
  }

  realft(data_ft, n+1, 1);

  for (size_t i=2; i<n+1; i+=2)
  {
  //  cout << -(M_PI)+i*delta_x << "\t" << 2.0/n * data_ft[i] << "\t" << 2.0/n * data_ft[i+1] << endl;
  }

  float* data_ft_derivative = new float[n+1];
  data_ft_derivative[0] = data_ft[0];
  data_ft_derivative[1] = data_ft[1];
  data_ft_derivative[2] = data_ft[2];

  for (size_t i=3; i<n; i+=2)
  {
      data_ft_derivative[i] = data_ft[i+1]*((i-1)/2);
      data_ft_derivative[i+1] = -data_ft[i]*((i-1)/2);
  }
  realft(data_ft_derivative, n+1,-1);

  for (size_t i=2; i<n+1; i+=2)
  {
    cout << -(M_PI)+i*delta_x << "\t" << 2.0/n * data_ft_derivative[i] << "\t" << 2.0/n * data_ft_derivative[i+1] << endl;
  }
}
