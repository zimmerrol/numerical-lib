#include "square_matrix.h"
#include <vector>

using namespace std;

#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H
namespace numerical
{
	size_t linear_system_solve_sor(square_matrix& coeff_matrix, vector<double> &x, vector<double> &b, double alpha, double error_threshold);
	double linear_system_solve_sor_step(square_matrix& coeff_matrix, vector<double> &x, vector<double> &b, double alpha);

	void multiply_matrix_vector(matrix& matrix, vector<double> &x, vector<double> &result);
}
#endif
