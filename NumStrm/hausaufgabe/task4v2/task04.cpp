#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include "../../../lib/linear_system.h"
#include "../../../lib/square_matrix.h"

using namespace std;
using namespace numerical;

//custom type which discribes a 2d grid. use this instead of arrays for access violation safety
typedef vector<vector<double > > grid_t;

/*calculates the variable beta. expects the parameters:
	plus: indicates whether beta plus or beta minus will be calculated
	delta_t, delta_x, delta_y, pe: the time-step-size, the step size in x/y direction and the Péclet number
	x,y: the x/y coordinate of the point on which beta is desired
*/
double beta(bool plus, double delta_t, double delta_x, double delta_y, double pe, double x, double y)
{
	if (plus)
	{
		return delta_t*(-1.0 / pow(delta_x, 2) + (pe*(M_PI*sin(2 * M_PI*x*delta_x)*cos(M_PI*y*delta_y))) / (2.0*delta_x));
	}
	else
	{
		return delta_t*(-1.0 / pow(delta_x, 2) - (pe*(M_PI*sin(2 * M_PI*x*delta_x)*cos(M_PI*y*delta_y))) / (2.0*delta_x));
	}
}

/*calculates the variable ganna. expects the parameters:
	plus: indicates whether gamma plus or gamma minus will be calculated
	delta_t, delta_x, delta_y, pe: the time-step-size, the step size in x/y direction and the Péclet number
	x,y: the x/y coordinate of the point on which gamma is desired
*/
double gamma(bool plus, double delta_t, double delta_x, double delta_y, double pe, double x, double y)
{
	if (plus)
	{
		return delta_t*(-1.0 / pow(delta_y, 2) + (pe*(-2.0*M_PI*cos(2 * M_PI*x*delta_x)*sin(M_PI*y*delta_y))) / (2.0*delta_y));
	}
	else
	{
		return delta_t*(-1.0 / pow(delta_y, 2) - (pe*(-2.0*M_PI*cos(2 * M_PI*x*delta_x)*sin(M_PI*y*delta_y))) / (2.0*delta_y));
	}
}

int main(int argc, char* argv[])
{
	//print information for the usage
  cout << "Use this program to calculate the temperature with the BTCS-scheme in the rectangle." << endl;
  cout << "Expects the follwing set of parameters:" << endl;

  cout << "\toutput file:\tThe path to the file in which the results will be stored." << endl;
  cout << "\tN_x:\t\t.The amount of grid points in the x-direction." << endl;
  cout << "\tN_y\t\tThe amount of grid points in the y-direction." << endl;
  cout << "\tPe:\t\tThe Péclet-Number" << endl;
  cout << "\tMax t\t\tThe maximum time until the system will be simulated." << endl;
  cout << "\tDelta t\t\tThe time step size for the integration." << endl;

  //read entered parameters

  //create output stream
  ofstream outputFile;
	outputFile.open(argv[1]);
  outputFile << fixed << setprecision(5);

	//the amount of grid points in the x direction
  const size_t dimension_x = (atoi)(argv[2])+1;
	//the distance between to grid points
  const double delta_x = 1.0/(dimension_x-1);

	//the amount of grid points in the x direction
  const size_t dimension_y = (atoi)(argv[3])+1;
	//the distance between to grid points
  const double delta_y = 1.0/(dimension_y-1);

	//Péclet-number
  const double pe = atof(argv[4]);
	//maximum time until which the system will be simulated
  const double max_t = atof(argv[5]);
	//the time step size
  const double delta_t = atof(argv[6]);

  cout << "Entered parameters:" << endl;
  cout << "\tNx: " << dimension_x << "\tNy: " << dimension_y << "\tPe: " << pe << "\tMax t: " << max_t << "\tDelta t: " << delta_t << endl;

	//the length of the value vector
	size_t n = dimension_x*dimension_y;

	//create a new square matrix (nXn) to store the matrix M
	square_matrix coeff_matrix(n, 0);

	//create const value for the variable alpha
  const double alpha = delta_t * (2.0 / pow(delta_x, 2) + 2.0 / pow(delta_y, 2)) + 1.0;

	//now set the values of the three matrices

	//unit matrix
	//main diagonal for the dirichlet boundaries again
	for (size_t i = 0; i<dimension_x; i++)
	{
		coeff_matrix.set_value(i, i, 1.0);
	}

	//Matrix M
	//contains the biggest part of the BTCS-scheme
	for (size_t i = dimension_x+1; i<n - dimension_x-1; i++)
	{
		//ingore this diagonals, as they affect the boundaries
		if (i % dimension_x == 0 || i % dimension_x == dimension_x-1)
		{
			continue;
		}

		//i % dimension_x gives the x coordinate of the field-value on which the current matrix's item operates
		//i / dimension_x gives the y coordinate of the field-value on which the current matrix's item operates
		coeff_matrix.set_value(i, i - dimension_x, gamma(false, delta_t, delta_x, delta_y, pe, i % dimension_x, i/dimension_x));
		coeff_matrix.set_value(i, i - 1, beta(false, delta_t, delta_x, delta_y, pe, i % dimension_x, i/dimension_x));
		coeff_matrix.set_value(i, i, alpha);
		coeff_matrix.set_value(i, i + 1, beta(true, delta_t, delta_x, delta_y, pe, i % dimension_x, i/dimension_x));
		coeff_matrix.set_value(i, i + dimension_x, gamma(true, delta_t, delta_x, delta_y, pe, i % dimension_x, i/dimension_x));
	}

	//now set the repeating line which take care of the neumann boundaries
	//first left boundary
	coeff_matrix.set_value(dimension_x, dimension_x+2, 1.0/3.0);
	coeff_matrix.set_value(dimension_x, dimension_x+1, -4.0/3.0);
	coeff_matrix.set_value(dimension_x, dimension_x, 1);
	for (size_t i=2; i<dimension_x-1; i++)
	{
		//right
		coeff_matrix.set_value(i*dimension_x-1, i*dimension_x-3, 1.0/3.0);
		coeff_matrix.set_value(i*dimension_x-1, i*dimension_x-2, -4.0/3.0);
		coeff_matrix.set_value(i*dimension_x-1, i*dimension_x-1, 1);

		//left
		coeff_matrix.set_value(i*dimension_x, i*dimension_x, 1);
		coeff_matrix.set_value(i*dimension_x, i*dimension_x+1, -4.0/3.0);
		coeff_matrix.set_value(i*dimension_x, i*dimension_x+2, 1.0/3.0);
	}

	//last right boundary
	coeff_matrix.set_value((dimension_x)*(dimension_x-1)-1, (dimension_x)*(dimension_x-1)-3, 1.0/3.0);
	coeff_matrix.set_value((dimension_x)*(dimension_x-1)-1, (dimension_x)*(dimension_x-1)-2, -4.0/3.0);
	coeff_matrix.set_value((dimension_x)*(dimension_x-1)-1,(dimension_x)*(dimension_x-1)-1, 1);

	//unit matrix
	//main diagonal for the dirichlet boundaries again
	for (size_t i = n - dimension_x; i<n; i++)
	{
		coeff_matrix.set_value(i, i, 1.0);
	}


	//create two vectors to store the temperature and to calculate the new one
	vector<double> values;
	vector<double> old_values;

	//fill them with the starting temperature
	for (size_t y = 0; y<dimension_y; y++)
	{
		for (size_t x = 0; x<dimension_x; x++)
		{
			old_values.push_back(y*delta_y);
			values.push_back(y*delta_y);
		}
	}

	//start to measure the runtime now
	clock_t start_time = clock();

	//loop over the time, until we have reached the maximum time
	for (double t = 0.0; t<max_t; t += delta_t)
	{
		//set the right hand side of the equation (set some values to zero for the neuman boundaries)
		old_values.at(dimension_x) = 0;
		for (size_t i=2; i<dimension_x-1; i++)
		{
			old_values.at(i*dimension_x-1) = 0;
			old_values.at(i*dimension_x) = 0;
		}
		old_values.at((dimension_x)*(dimension_x-1)-1) = 0;

		//invert the matrix M/solve the linear system to get the new temperature(value)
		cout << "t: " << t << "\titerations: " << linear_system_solve_sor(coeff_matrix, values, old_values, 1.455, 1e-4) << endl;

		//copy the new values into the old_values vector for the next iteration
		for (size_t i = 0; i<n; i++)
		{
			old_values.at(i) = values.at(i);
		}
	}

	//measure the end time
	clock_t end_time = clock();

	//print the runtime
	cout << "Integration finished after " << (end_time - start_time) / CLOCKS_PER_SEC << " seconds.";

	//print the final grid's values in the givven output file in a matrix form
	for (size_t i = 0; i<n; i++)
	{
		outputFile << values.at(i) << "\t";
		if (i > 0 && (i+1) % dimension_x == 0)
		{
			outputFile << endl;
		}
	}

	return 1;
}
