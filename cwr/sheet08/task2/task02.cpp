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
  if (x == 0 || x == x_dimension-1 || y == 0 || y == y_dimension-1 )
  {
    return true;
  }

  //get electrode position in the grid and check for collision
  if (sqrt(pow(x_position_electrode_1 - x * delta_x,2) + pow(y_position_electrode_1 - y * delta_y,2)) < r_electrode_1)
  {
    return true;
  }

  if (sqrt(pow(x_position_electrode_2 - x * delta_x,2) + pow(y_position_electrode_2 - y * delta_y,2)) < r_electrode_2)
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
      x_values[x][y] = (-b_values[x+y*x_dimension] + (x_values[x+1][y] + x_values[x-1][y])/pow(delta_x,2) + (x_values[x][y+1] + x_values[x][y-1])/pow(delta_y,2)) / ( 2* (1/pow(delta_x,2) + 1/pow(delta_y,2)));
    }
  }
}

double calculate_error(double** x_values, double* b_values, int x_dimension, int y_dimension, double delta_x, double delta_y)
{
  double squared_error = 0;

  for (int x = 1; x<x_dimension-1; x++)
  {
    for (int y = 1; y<y_dimension-1; y++)
    {
      squared_error += pow(b_values[x_dimension*y + x] - 1/pow(delta_x,2) * (x_values[x+1][y] + x_values[x-1][y] - 2*x_values[x][y]) - 1/pow(delta_y,2) * (x_values[x][y+1] + x_values[x][y-1] - 2*x_values[x][y]) ,2);
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

  double delta_x = length_x / dimension_x;
  double delta_y = length_y / dimension_y;

  double** values = new double*[dimension_x];

  for (int x = 0; x<dimension_x; x++)
  {
    values[x] = new double[dimension_y];
  }

  for (int x=0; x<dimension_x; x++)
  {
    for (int y=0; y<dimension_y; y++)
    {
      values[x][y] = (rand() % 2000) / 1000.0 - 1.0;
    }
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
  for (int x = 0; x<dimension_x; x++)
  {
    for (int y=0; y<dimension_y; y++)
    {
      if (sqrt(pow(x_position_electrode_1 - x * delta_x,2) + pow(y_position_electrode_1 - y * delta_y,2)) < r_electrode_1)
      {
        values[x][y] = v_electrode_1;
      }
    }
  }

  //electrode 2
  for (int x = 0; x<dimension_x; x++)
  {
    for (int y=0; y<dimension_y; y++)
    {
      if (sqrt(pow(x_position_electrode_2 - x * delta_x,2) + pow(y_position_electrode_2 - y * delta_y,2)) < r_electrode_2)
      {
        values[x][y] = v_electrode_2;
      }
    }
  }

  int step = 0;
  int aa = 0;
  double err = calculate_error(values,  b, dimension_x, dimension_y, delta_x, delta_y);
  while (err > epsilon)
  {
    aa++;
    if (aa == 25)
    {
      cout << "error: " << err << endl;
      aa = 0;
    }
    gauss_seidel_step(values,  b, dimension_x, dimension_y, delta_x, delta_y);
    err = calculate_error(values,  b, dimension_x, dimension_y, delta_x, delta_y);

    step++;

    if (step == 25*200)
    {
      break;
    }
  }

  //calc E
  for (int x = 0; x<dimension_x-1; x++)
  {
    for (int y = 0; y<dimension_y-1; y++)
    {
      outputFile << x * delta_x << "\t" << y * delta_y << "\t" << (values[x+1][y] - values[x][y]) << "\t" << (values[x][y+1] - values[x][y]) << endl;
      //outputFile << x << "\t" << y << "\t" << values[x][y] << endl;
    }
  }


  cout << "final error=" <<  err << endl;


/*  cout << "field:\n";
  for (int x = 0; x<dimension_x; x++)
  {
    for (int y = 0; y<dimension_y; y++)
    {
      cout << values[x][y] << "\t";
    }
    cout << endl;
  }*/

}
