#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <complex>
#include <cmath>
#include "../../../lib/diff.h"

using namespace std;

double real_diff_square(double x)
{
  return 2*x;
}

double real_diff_cube(double x)
{
  return 3*x*x;
}

double square(double x, double* args)
{
  return x*x;
}

double cube(double x, double* args)
{
  return x*x*x;
}

int main(int argc, char* argv[])
{
  ofstream output;
  output.open(argv[1]);
  output << fixed << setprecision(35);

  double delta= atof(argv[2]);

  for (double h = delta;h< 1;h+=delta)
  {
    double st_error_square = 0.0;
    double st_error_cube = 0.0;

    double fw_error_square = 0.0;
    double fw_error_cube = 0.0;
    for (double x = -2.0; x <2.0;x+=h)
    {
      double st_dx_square = numerical::differentiate_centered_difference(&square, x, h, NULL);
      double st_dx_cube = numerical::differentiate_centered_difference(&cube, x, h, NULL);
      st_error_square += pow((st_dx_square - real_diff_square(x)),2);
      st_error_cube += pow((st_dx_cube - real_diff_cube(x)),2);

      double fw_dx_square = numerical::differentiate_newton_gregroy_forwards(&square, x, h, NULL);
      double fw_dx_cube = numerical::differentiate_newton_gregroy_forwards(&cube, x, h, NULL);
      fw_error_square += pow((fw_dx_square - real_diff_square(x)),2);
      fw_error_cube += pow((fw_dx_cube - real_diff_cube(x)),2);
    }

    st_error_square /= delta / 2.0;
    st_error_cube *= delta / 2.0;

    fw_error_square *= delta / 2.0;
    fw_error_cube *= delta / 2.0;

    output << h << "\t" << st_error_square << "\t" << fw_error_square<< "\t" << st_error_cube << "\t" << fw_error_cube << "\n";
  }

}
