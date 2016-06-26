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
  size_t n = 64;
  double max_t = 4.2;
  double delta_t = 1e-3;
  double D = 1.0;

  Complex functionData[n];

  double delta_x = 1.0/n;

  for (size_t i=0; i<n; i++)
  {
    functionData[i] = sin(2*PI*i*delta_x);
  }

  // forward fft
  CArray data(functionData, n);


  for (double t=0; t<max_t; t+=delta_t)
  {
    data[0] = 0;
    data[n] = 0;

    fft(data);
    Complex second_derivation;
    Complex k;
    for (size_t i=0; i<n; i++)
    {
      k = (Complex)(2.0* i / n -1);
      second_derivation = - k*k*data[i];
      data[i] += D*delta_t*second_derivation;
    }

    ifft(data);
  }

  for (size_t i = 0; i < n; ++i)
  {
     std::cout << i*delta_x << "\t" << data[i].real() << std::endl;
  }

  return 0;
}
