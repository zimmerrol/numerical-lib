#include "square_matrix.h"
#include <iostream>
#include <fstream>

namespace numerical
{
  //solves M*x = b for M €tridiagMatrix(...)
  //implementation of the thomas algortihm
  void linear_system_solve_triagonal(square_matrix& matrix, double*x, double* b)
  {
    double* c = new double[matrix.get_n()];
    double* d = new double[matrix.get_n()];

    c[0] = matrix.get_value(0,1) / matrix.get_value(0,0);
    d[0] = b[0] / matrix.get_value(0,0);
    for (size_t i=1; i<matrix.get_n();i++)
    {
        c[i] = matrix.get_value(i,i+1) / (matrix.get_value(i,i) - c[i-1]*matrix.get_value(i,i-1));
        d[i] = (b[i] - d[i-1]*matrix.get_value(i-1,i)) / (matrix.get_value(i,i) - c[i-1]*matrix.get_value(i,i-1));
    }

    x[matrix.get_n()-1] = d[matrix.get_n()-1];

    for (size_t i = matrix.get_n()-1; i>0; i--)
    {
      x[i] = d[i] - c[i]*x[i+1];
    }
    x[0] = d[0] - c[0]*x[1];

    delete c;
    delete d;
  }
}
