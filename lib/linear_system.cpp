#include "linear_system.h"
#include "square_matrix.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <float.h>

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


	void multiply_matrix_vector(matrix& matrix, vector<double> &xVec, vector<double> &result)
	{
		result.clear();
		double val;
		for (size_t y = 0; y<matrix.get_n(); y++)
		{
			val = 0.0;
			for (size_t x = 0; x<matrix.get_m(); x++)
			{
				val += matrix.get_value(y, x) * xVec.at(x);
			}
			result.push_back(val);
		}
	}

	//calculates one SOR step the for system coeff_matrix * x = b with the relaxation parameter alpha
	double linear_system_solve_sor_step(square_matrix& coeff_matrix, vector<double> &x, vector<double> &b, double alpha)
	{
		//the new value of the current item
		double new_value;
		//a helper variable to store a needed sum of items in the loop
		double sum;

		//loop over all items and calculate the next elemt of x
		for (size_t k = 0; k<coeff_matrix.get_n(); k++)
		{
			sum = 0.0;

			//splitted the one loop into two seperate loops for better performance.
			//one of the loops is already using the new values while the other one is using the old one
			for (size_t i = 0; i<k; i++)
			{
				sum += coeff_matrix.get_value(k, i)*x.at(i);
			}
			for (size_t i = k + 1; i<coeff_matrix.get_n(); i++)
			{
				sum += coeff_matrix.get_value(k, i)*x.at(i);
			}
			new_value = (1.0 - alpha)*x.at(k) + alpha / coeff_matrix.get_value(k, k)*(b.at(k) - sum);

			//update the new x value
			x.at(k) = new_value;
		}

		//calculate the residue now
		double residue = 0.0;
		vector<double> testVector;
		multiply_matrix_vector(coeff_matrix, x, testVector);
		for (size_t k = 0; k<coeff_matrix.get_n(); k++)
		{
			residue += pow(testVector.at(k) - b.at(k), 2);
		}

		//return the length of the residue vector
		return sqrt(residue);
	}

	//tries to solve the linear equation coeff_matrix * x = b for x iterative with the SOR method untill either:
	//1) the residues are lower than the specified error_threshold
	//2) the residues of 50 following iterations have been increasing => propably no convergence for this system
	//returns the needed iterations

	size_t linear_system_solve_sor(square_matrix& coeff_matrix, vector<double> &x, vector<double> &b, double alpha, double error_threshold)
	{
		//contains the value of the residues of the current iteration. start with a high number because of the if-comparison below
		double error = DBL_MAX;
		//contains the value of the residues of the previous iteration. start with a high number because of the if-comparison below
		double old_error;
		//counter of the needed iterations
		size_t iteration = 0;

		//contains how often the residues of two following iterations have been increasing
		size_t residue_increased_counter = 0;

		//loop, until a solution has been found
		while (true)
		{
			//swap the old and the current error value
			old_error = error;
			//increase the number of needed iterations
			iteration++;

			//calculate the next SOR step
			error = linear_system_solve_sor_step(coeff_matrix, x, b, alpha);


			//	cout << "iteration: " << iteration << "\terror:" << error << "\n";


			//check for exit condition 1)
			if (error < error_threshold)
			{
				return iteration;
			}
			else if (error > old_error)
			{
				//exit condition 2
				if (residue_increased_counter > 50){
					return iteration;
				}
				residue_increased_counter++;
			}
			else
			{
				//reset the counter for the second exit condition.
				residue_increased_counter=0;
			}
		}
	}
}
