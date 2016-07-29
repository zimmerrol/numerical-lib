#include "linear_system.h"
#include "square_matrix.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

namespace numerical
{
	//solves M*x = b for M €tridiagMatrix(...)
	//implementation of the thomas algortihm
	void linear_system_solve_tridiagonal(square_matrix& matrix, double*x, double* b)
	{
		double* c = new double[matrix.get_n()];
		double* d = new double[matrix.get_n()];

		c[0] = matrix.get_value(0, 1) / matrix.get_value(0, 0);
		d[0] = b[0] / matrix.get_value(0, 0);
		for (size_t i = 1; i<matrix.get_n(); i++)
		{
			c[i] = matrix.get_value(i, i + 1) / (matrix.get_value(i, i) - c[i - 1] * matrix.get_value(i, i - 1));
			d[i] = (b[i] - d[i - 1] * matrix.get_value(i - 1, i)) / (matrix.get_value(i, i) - c[i - 1] * matrix.get_value(i, i - 1));
		}

		x[matrix.get_n() - 1] = d[matrix.get_n() - 1];

		for (size_t i = matrix.get_n() - 1; i>0; i--)
		{
			x[i] = d[i] - c[i] * x[i + 1];
		}
		x[0] = d[0] - c[0] * x[1];

		delete c;
		delete d;
	}

	//solves M*x = b for M €tridiagMatrix(...)
	//implementation of the thomas algortihm
	void linear_system_solve_tridiagonal(size_t n, double* diag_left, double* diag, double* diag_right, double*x, double* b)
	{
		double* c = new double[n];
		double* d = new double[n];

		c[0] = diag_right[0] / diag[0];
		d[0] = b[0] / diag[0];
		for (size_t i = 1; i<n; i++)
		{
			c[i] = diag_right[i] / (diag[i] - c[i - 1] * diag_left[i]);
			d[i] = (b[i] - d[i - 1] * diag_right[i - 1]) / (diag[i] - c[i - 1] * diag_left[i]);
		}

		x[n - 1] = d[n - 1];

		for (size_t i = n - 1; i>0; i--)
		{
			x[i] = d[i] - c[i] * x[i + 1];
		}
		x[0] = d[0] - c[0] * x[1];

		delete c;
		delete d;
	}

	double linear_system_solve_gauss_seidel_step(square_matrix& coeff_matrix, double* x, double* b)
	{
		//this is just a special case of SOR with alpha = 1.0
		return linear_system_solve_sor_step(coeff_matrix, x, b, 1.0);
	}

	double linear_system_solve_sor_step(square_matrix& coeff_matrix, double*x, double* b, double alpha)
	{
		double error = 0.0;
		double new_value;
		double sum;
		for (size_t k = 0; k<coeff_matrix.get_n(); k++)
		{
			sum = 0.0;
			//splitted the one loop into two seperate loops for better performance
			for (size_t i = 0; i<k; i++)
			{
				sum += coeff_matrix.get_value(k, i)*x[i];
			}
			for (size_t i = k + 1; i<coeff_matrix.get_n(); i++)
			{
				sum += coeff_matrix.get_value(k, i)*x[i];
			}
			new_value = (1.0 - alpha)*x[k] + alpha / coeff_matrix.get_value(k, k)*(b[k] - sum);

			error += pow(new_value - x[k], 2);
			x[k] = new_value;
		}

		return error;
	}


	bool linear_system_solve_gauss_seidel(square_matrix& coeff_matrix, double*x, double* b, double error_threshold)
	{
		error_threshold = pow(error_threshold, 2);
		double error = -1.0;
		double old_error;
		while (true)
		{
			old_error = error;
			error = linear_system_solve_gauss_seidel_step(coeff_matrix, x, b);

			if (error < error_threshold)
			{
				return true;
			}
			else if (error > old_error)
			{
				return false;
			}
		}
	}


	bool linear_system_solve_sor(square_matrix& coeff_matrix, double*x, double* b, double alpha, double error_threshold)
	{
		error_threshold = pow(error_threshold, 2);
		double error = -1.0;
		double old_error;
		while (true)
		{
			old_error = error;
			error = linear_system_solve_sor_step(coeff_matrix, x, b, alpha);

			if (error < error_threshold)
			{
				return true;
			}
			else if (error > old_error)
			{
				return false;
			}
		}
	}


	bool linear_system_solve_sor2(square_matrix& coeff_matrix, vector<double> &x, vector<double> &b, double alpha, double error_threshold)
	{
		double error = 1000000000;
		double old_error;
		while (true)
		{
			old_error = error;
			error = linear_system_solve_sor_step2(coeff_matrix, x, b, alpha);

			if (error < error_threshold)
			{
				return true;
			}
			else if (error > old_error)
			{
				return false;
			}
		}
	}


	double linear_system_solve_sor_step2(square_matrix& coeff_matrix, vector<double> &x, vector<double> &b, double alpha)
	{
		double error = 0.0;
		double new_value;
		double sum;
		for (size_t k = 0; k<coeff_matrix.get_n(); k++)
		{
			sum = 0.0;
			//splitted the one loop into two seperate loops for better performance
			for (size_t i = 0; i<k; i++)
			{
				sum += coeff_matrix.get_value(k, i)*x.at(i);
			}
			for (size_t i = k + 1; i<coeff_matrix.get_n(); i++)
			{
				sum += coeff_matrix.get_value(k, i)*x.at(i);
			}
			new_value = (1.0 - alpha)*x.at(k) + alpha / coeff_matrix.get_value(k, k)*(b.at(k) - sum);

			error += pow(new_value - x.at(k), 2);
			x.at(k) = new_value;
		}

		return sqrt(error);
	}



}
