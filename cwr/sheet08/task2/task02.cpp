#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "../../../lib/ode.h"

using namespace std;

const double x_position_electrode_1 = 0.4;
const double y_position_electrode_1 = 0.6;
const double r_electrode_1 = 0.2;
const double v_electrode_1 = -1.0;

const double x_position_electrode_2 = 1.1;
const double y_position_electrode_2 = 0.4;
const double r_electrode_2 = 0.3;
const double v_electrode_2 = 1.0;

const double length_x = 1.6;
const double length_y = 1.0;


bool point_has_fixed_border_constraint(int x, int y, int x_dimension, int y_dimension, double delta_x, double delta_y)
{
  //outer border values
  if (x == 0 || x == x_dimension || y == 0 || y == y_dimension )
  {
    return true;
  }

  //get electrode position in the grid and check for collision
  if (sqrt(pow(x_position_electrode_1 - x * delta_x,2) + pow(y_position_electrode_1 - x * delta_y,2)) < r_electrode_1)
  {
    return true;
  }

  if (sqrt(pow(x_position_electrode_2 - x * delta_x,2) + pow(y_position_electrode_2 - x * delta_y,2)) < r_electrode_2)
  {
    return true;
  }

  //nothing detected
  return false;
}

void gauss_seidel_step(double** x_values, double* b_values, int x_dimension, int y_dimension, double delta_x, double delta_y)
{
  //skip the values at the border, as they are fixed
  for (int x = 0; x<x_dimension; x++)
  {
    for (int y = 0; y<y_dimension; y++)
    {
      if (point_has_fixed_border_constraint(x,y,x_dimension, y_dimension, delta_x, delta_y))
      {
        continue;
      }
      //laplace u = 1/h^2 * (u_i+1,j + u_i-1,j + u_i,j+0 + u_i,j-1 - 4*u_i,j) for delta_x = delta_y
      //x_values[x][y] = 1/pow(delta_x,2) * (x_values[x+1][y] + x_values[x-1][y] - 2*x_values[x][y]) + 1/pow(delta_y,2) * (x_values[x][y+1] + x_values[x][y-1] - 2*x_values[x][y]);
      x_values[x][y] = (-b_values[x+y] + (x_values[x+1][y] + x_values[x-1][y])/pow(delta_x,2) + (x_values[x][y+1] + x_values[x][y-1])/pow(delta_y,2)) / (pow(delta_x,2) + pow(delta_y,2));
    }
  }
}

double calculate_error(double** x_values, double* b_values, int x_dimension, int y_dimension, double delta_x, double delta_y)
{
  double squared_error = 0;

  for (int x = 0; x<x_dimension; x++)
  {
    for (int y = 0; y<y_dimension; y++)
    {
      squared_error += pow(b_values[x+y] - 1/pow(delta_x,2) * (x_values[x+1][y] + x_values[x-1][y] - 2*x_values[x][y]) - 1/pow(delta_y,2) * (x_values[x][y+1] + x_values[x][y-1] - 2*x_values[x][y]) ,2);
    }
  }

  return sqrt(squared_error);
}


int main(int argc, char* argv[])
{
	ofstream outputFile;
	outputFile.open(argv[1]);
	outputFile << fixed << setprecision(5);

  double epsilon = atof(argv[2]);

  int dimension_x = atoi(argv[3]);
  int dimension_y = atoi(argv[4]);

  double** values = new double*[dimension_x];
  for (int x = 0; x<dimension_x; x++)
  {
    values[x] = new double[dimension_y];
  }

  double* b = new double[dimension_x * dimension_y];
  for (int i=0; i < dimension_x * dimension_y; i++)
  {
    b[i] = 0;
  }

  //set border constraints
  //top. bottom
  for (int x=0; x<dimension_x; x++)
  {
    values[x][0] = 0;
    values[x][dimension_y-1] = 0;
  }

  //left, right
  for (int y=0; y<dimension_y; y++)
  {
    values[0][y] = 0;
    values[dimension_x-1][y] = 0;
  }

  //electrode 1

  //electrode 2

}
