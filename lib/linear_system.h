#include "square_matrix.h"

#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H
namespace numerical
{
  void linear_system_solve_tridiagonal(size_t n, double* diag_left, double* diag, double* diag_right, double*x, double* b);
  void linear_system_solve_tridiagonal(square_matrix& matrix, double*x, double* b);
  double linear_system_solve_gauss_seidel_step(square_matrix& coeff_matrix, double*x, double* b);
  double linear_system_solve_sor_step(square_matrix& coeff_matrix, double*x, double* b, double alpha);

  bool linear_system_solve_gauss_seidel(square_matrix& coeff_matrix, double*x, double* b, double error_threshold);
  bool linear_system_solve_sor(square_matrix& coeff_matrix, double*x, double* b, double alpha, double error_threshold);
}
#endif
