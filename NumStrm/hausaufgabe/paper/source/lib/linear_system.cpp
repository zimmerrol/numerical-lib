#include "linear_system.h"
#include "square_matrix.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

namespace numerical
{
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
		double error = 1000000000;
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
