#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <complex>
#include <cmath>
#include "../../../lib/diff.h"

using namespace std;

double real_diff_func(double x)
{
  return 5*cos(5*x);
}

double func(double x, double* args)
{
  return sin(5*x);
}

int main(int argc, char* argv[])
{
  ofstream output_diff;
  output_diff.open(argv[1]);
  output_diff << fixed << setprecision(35);

  ofstream output;
  output.open(argv[2]);
  output << fixed << setprecision(35);

  double delta = 4*M_PI/100;

  double st_error = 0.0;
  double fw_error = 0.0;
  for (double x = 0; x <4*M_PI;x+=delta)
  {
    double st_dx = numerical::differentiate_centered_difference(&func, x, delta, NULL);
    double fw_dx= numerical::differentiate_newton_gregroy_forwards(&func, x, delta, NULL);
    st_error += pow((st_dx - real_diff_func(x)),2);
    fw_error += pow((fw_dx - real_diff_func(x)),2);

    output << x << "\t" << real_diff_func(x) << "\t" << st_dx << "\t" << fw_dx << "\n";

  }

  st_error /= delta/ 2.0;
  fw_error *= delta / 2.0;

  output_diff << delta << "\t" << st_error << "\t" << fw_error << "\n";
}
