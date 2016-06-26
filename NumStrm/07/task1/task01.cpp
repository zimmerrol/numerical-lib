#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <iostream>
#include <valarray>

#include "fft_fn.h"

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
  size_t n = pow(2,10);
  Complex functionData[n];

  double delta_x = 2*PI/n;

  for (size_t i=0; i<n; i++)
  {
    functionData[i] = sin(i*delta_x);
  }


  CArray data(functionData, n);

  // forward fft
  fft(data);

  Complex imaginary(0,1);

  CArray d1data(functionData, n);
  CArray d2data(functionData, n);

  for (int i = 0; i < n; ++i)
  {
    Complex k(2.0* i / n -1,0);
    d1data[i] = -imaginary*k*data[i];
    d2data[i] = (-imaginary*k)*(-imaginary*k)*data[i];
  }

  // inverse fft
  ifft(data);
  ifft(d1data);
  ifft(d2data);

  for (int i = 0; i < n; ++i)
  {
     std::cout << i*delta_x << "\t" << data[i].real() << "\t" << d1data[i].real() << "\t" << d2data[i].real() << std::endl;
  }
  return 0;
}
