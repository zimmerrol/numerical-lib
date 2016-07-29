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


typedef vector<vector<double > > grid_t;

double beta(bool plus, double delta_x, double delta_y, double pe, double x, double y)
{
  if (plus)
  {
    return -1.0/pow(delta_x,2) + (pe*(M_PI*sin(2*M_PI*x*delta_x)*cos(M_PI*y*delta_y)))/(2.0*delta_x);
  }
  else
  {
    return -1.0/pow(delta_x,2) - (pe*(M_PI*sin(2*M_PI*x*delta_x)*cos(M_PI*y*delta_y)))/(2.0*delta_x);
  }
}

double gamma(bool plus, double delta_x, double delta_y, double pe, double x, double y)
{
  if (plus)
  {
    return -1.0/pow(delta_y,2) + (pe*(-2.0*M_PI*cos(2*M_PI*x*delta_x)*sin(M_PI*y*delta_y)))/(2.0*delta_y);
  }
  else
  {
    return -1.0/pow(delta_y,2) - (pe*(-2.0*M_PI*cos(2*M_PI*x*delta_x)*sin(M_PI*y*delta_y)))/(2.0*delta_y);
  }
}

int main(int argc, char* argv[])
{
	//print information for the usage
	cout << "" << endl;
	cout << "Expects the follwing set of parameters:" << endl;

	cout << "\toutput file:\tThe path to the file in which the results will be stored." << endl;
	cout << "\tN_x:\t\t." << endl;
	cout << "\tN_y\t\t" << endl;
	cout << "\tPe:\t\t" << endl;
	cout << "\tMax t\t\t" << endl;
	cout << "\tDelta t\t\t" << endl;

	ofstream outputFile;
	outputFile.open(argv[1]);
	outputFile << fixed << setprecision(5);

	size_t dimension_x = (atoi)(argv[2]) + 1;
	double delta_x = 1.0 / (dimension_x - 1);

	size_t dimension_y = (atoi)(argv[3]) + 1;
	double delta_y = 1.0 / (dimension_y - 1);

	double pe = atof(argv[4]);

	double max_t = atof(argv[5]);
	double delta_t = atof(argv[6]);

	cout << dimension_x << "\t" << dimension_y << "\t" << pe << "\t" << max_t << "\t" << delta_t << endl;

	//ftcs_time_step(grid_t values, grid_t v0x, grid_t v0y, double delta_t, double delta_x, double delta_y, double pe)

	size_t n = dimension_x*dimension_y;
	square_matrix coeff_matrix(n, 0);

	//Matrix A
	/*for (size_t i=1; i<dimension_y-1; i++)
	{
	//top left
	coeff_matrix.set_value(i, i, 1.0);
	coeff_matrix.set_value(i, i + dimension_y, -4.0/3.0);
	coeff_matrix.set_value(i, i+ dimension_y*2, 1.0/3.0);


	//bottom right
	coeff_matrix.set_value(i + (n-dimension_x+1), i, 1.0);
	coeff_matrix.set_value(i + (n-dimension_x+1), i + (n-dimension_x+1) + dimension_y, -4.0/3.0);
	coeff_matrix.set_value(i + (n-dimension_x+1), i + (n-dimension_x+1) + dimension_y*2, 1.0/3.0);
	}

	coeff_matrix.set_value(dimension_y+1, dimension_y+1, 1.0);
	for (size_t i=dimension_y+1; i<n-dimension_y+1; i++)
	{
	coeff_matrix.set_value(i, i-1, -1.0/pow(delta_y,2) - pe/(2.0*delta_y) * );
	coeff_matrix.set_value(i, i, 2.0/pow(delta_x,2) + 2.0/pow(delta_y,2) + 1.0/delta_t);
	coeff_matrix.set_value(i, i+1, 1.0);
	}*/


	double alpha = delta_t * (2.0 / pow(delta_x, 2) + 2.0 / pow(delta_y, 2)) + 1.0;

	for (size_t i = 0; i<dimension_x; i++)
	{
		coeff_matrix.set_value(i, i, 1.0);
	}

  cout << "b\n";

	coeff_matrix.set_value(0, 1, -4.0 / 3.0);
	coeff_matrix.set_value(0, 2, 1.0 / 3.0);

  cout << "c\n";
	for (size_t i = dimension_x; i<n; i += dimension_x)
	{
		coeff_matrix.set_value(0, i - 2, 1.0 / 3.0);
		coeff_matrix.set_value(0, i - 1, -4.0 / 3.0);
		coeff_matrix.set_value(0, i + 1, -4.0 / 3.0);
		coeff_matrix.set_value(0, i + 2, 1.0 / 3.0);
	}

  cout << "d\n";

	coeff_matrix.set_value(0, dimension_x - 1, -4.0 / 3.0);
	coeff_matrix.set_value(0, dimension_x - 2, 1.0 / 3.0);

  cout << "e\n";
	for (size_t i = dimension_x; i<n - dimension_x; i++)
	{
		coeff_matrix.set_value(i, i - dimension_x, gamma(false, delta_x, delta_y, pe, i / dimension_x, i%dimension_x));
		coeff_matrix.set_value(i, i - 1, beta(false, delta_x, delta_y, pe, i / dimension_x, i%dimension_x));
		coeff_matrix.set_value(i, i, alpha);
		coeff_matrix.set_value(i, i + 1, beta(true, delta_x, delta_y, pe, i / dimension_x, i%dimension_x));
		coeff_matrix.set_value(i, i + dimension_x, gamma(true, delta_x, delta_y, pe, i / dimension_x, i%dimension_x));
	}
  cout << "f\n";

	for (size_t i = n - dimension_x; i<n; i++)
	{
		coeff_matrix.set_value(i, i, 1.0);
	}
  cout << "g\n";

	//grid_t v0x;
	//grid_t v0y;
	//grid_t values;
	vector<double> values;
	vector<double> old_values;
	for (size_t y = 0; y<dimension_y; y++)
	{
		//v0x.push_back(vector<double>());
		//v0y.push_back(vector<double>());
		//values.push_back(vector<double>());

		for (size_t x = 0; x<dimension_x; x++)
		{
			//v0x.at(x).push_back(M_PI*sin(2*M_PI*x*delta_x)*cos(M_PI*y*delta_y));
			//v0y.at(x).push_back(-2.0*M_PI*cos(2*M_PI*x*delta_x)*sin(M_PI*y*delta_y));
			old_values.push_back(y*delta_y);
			values.push_back(0.0);
		}
	}

	for (double t = 0.0; t<max_t; t += delta_t)
	{
		linear_system_solve_sor2(coeff_matrix, values, old_values, 1.0, 0.0001);

		for (size_t i = 0; i<n; i++)
		{

			old_values.at(i) = values.at(i);
		}

		//ftcs_time_step(values, v0x, v0y, sources, delta_t, delta_x, delta_y, pe) ;

	}

	for (size_t i = 0; i<n; i++)
	{
		outputFile << values.at(i) << "\t";
		if (i > 0 && (i+1) % dimension_x == 0)
		{
			outputFile << endl;
		}
	}

	/*for (size_t y=0; y<dimension_y; y++)
	{
	for (size_t x=0; x<dimension_x; x++)
	{
	outputFile << values.at(x + y*dimension_y) << "\t";//x*delta_x << "\t" << y*delta_y << "\t" << values.at(x).at(y) << endl;
	}
	outputFile << "\n";
	}*/

	return 1;
}
