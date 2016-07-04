#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <iostream>
#include <valarray>

using namespace std;

const double PI = 3.141592653589793238460;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

// Cooley–Tukey FFT (in-place, divide-and-conquer)
// Higher memory requirements and redundancy although more intuitive
void fft(CArray& x)
{
    const size_t N = x.size();
    if (N <= 1) return;

    // divide
    CArray even = x[std::slice(0, N/2, 2)];
    CArray  odd = x[std::slice(1, N/2, 2)];

    // conquer
    fft(even);
    fft(odd);

    // combine
    for (size_t k = 0; k < N/2; ++k)
    {
        Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
        x[k    ] = even[k] + t;
        x[k+N/2] = even[k] - t;
    }
}

// inverse fft (in-place)
void ifft(CArray& x)
{
    // conjugate the complex numbers
    x = x.apply(std::conj);

    // forward fft
    fft( x );

    // conjugate the complex numbers again
    x = x.apply(std::conj);

    // scale the numbers
    x /= x.size();
}

int main(int argc, char* argv[])
{
  size_t n = 256;
  double max_t = atof(argv[1]);
  double delta_t = 1e-3;
  double D = 0.01;
  double max_x = 2*PI;

  Complex functionData[n];

  double delta_x = max_x/n;

  for (size_t i=0; i<n; i++)
  {
    functionData[i] = sin(i*delta_x);
  }

  Complex imaginary(0,1);

  // forward fft
  CArray data(functionData, n);

  //fft
  //1. abl
  //ifft
  //non-lin
  //fft
  //euler-step
  //ifft

int j=0;
  for (double t=0; t<max_t; t+=delta_t)
  {
    j++;
    data[0] = 0;
    data[n-1] = 0;
    CArray derivative_data(n);
    CArray second_derivation_data(n);

    //fft
    fft(data);

    //calc (f')^
    for (size_t i = 0; i < n; ++i)
    {
      Complex k(2.0* i / n -1,0);
      derivative_data[i] = imaginary*k*data[i];
      second_derivation_data = -k*k*data[i];
    }

    // inverse fft to calc f'
    ifft(data);
    ifft(derivative_data);

    CArray non_linearity_data(n);
    for (size_t i = 0; i < n; ++i)
    {
      non_linearity_data[i] = derivative_data[i] * data[i];
    }

    //get transformed non-linearity
    fft(non_linearity_data);
    fft(data);

    for (size_t i = 0; i < n; ++i)
    {
        data[i] += delta_t * (D*second_derivation_data[i] - non_linearity_data[i]);
    }

    ifft(data);
  }

  for (size_t i = 0; i < n; ++i)
  {
     std::cout << i*delta_x << "\t" << data[i].real() << std::endl;
  }

  return 0;
}
