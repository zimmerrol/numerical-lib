#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "../../../lib/ode.h"

using namespace std;

const double x_position_electrode_1 = 0.4;
const double y_position_electrode_1 = 0.6;
const double size_electrode_1 = 0.2;
const double v_electrode_1 = -1.0;

const double x_position_electrode_2 = 1.1;
const double y_position_electrode_2 = 0.4;
const double size_electrode_2 = 0.3;
const double v_electrode_2 = 1.0;

const double length_x = 1.6;
const double length_y = 1.0;

class point2d
{
public:
    double x,y;
    point2d(){}
    point2d(double x, double y)
    {
      this-> x = x;
      this-> y = y;
    }
};

int jordan_cross_product_test(point2d a, point2d b, point2d c)
{
  if (a.y == b.y && b.y == c.y)
  {
    if ((b.x <= a.x && a.x <= c.x) || (c.x <= a.x && a.x <= b.x))
    {
      return 0;
    }
    else
    {
      return 1;
    }
  }

  if (b.y > c.y)
  {
    double tmp_x = b.x, tmp_y = b.y;
    b.x = c.x;
    b.y = c.y;
    c.x = tmp_x;
    c.y = tmp_y;
  }

  if (a.y == b.y || a.x == b.x)
  {
    return 0;
  }

  if (a.y <= b.y || a.y > c.y)
  {
    return 1;
  }

  double delta = (b.x-a.x) * (c.y-a.y) - (b.y-a.y) * (c.x-a.x);
  if (delta > 0)
  {
    return -1;
  }
  else if (delta < 0)
  {
    return 1;
  }
  else{
    return 0;
  }
}

int jordan_point_in_polygon_test(size_t n, point2d* points, point2d test_point)
{
  int t = -1;

  for (int i=0;i<n-1;i++)
  {
    t *= jordan_cross_product_test(test_point, points[i], points[i+1]);
  }

  t *= jordan_cross_product_test(test_point, points[n-1], points[0]);

  return t;
}

bool point_is_in_hexagon(double x, double y, double middle_x, double middle_y, double side_length)
{
  point2d* points = new point2d[3];

  points[0] = point2d(middle_x, middle_y);
  points[1] = point2d(middle_x - side_length, middle_y);
  points[2] = point2d(middle_x - side_length/2, middle_y + sqrt(3)/2 * side_length);

  if (jordan_point_in_polygon_test(3,points, point2d(x,y)))
  {
    cout << x << "\t" << y << "\t" << middle_x << "\t" << middle_y << endl;
    return true;
  }

  points[0] = point2d(middle_x, middle_y);
  points[1] = point2d(middle_x - side_length/2, middle_y + sqrt(3)/2 * side_length);
  points[2] = point2d(middle_x + side_length/2, middle_y + sqrt(3)/2 * side_length);

  if (jordan_point_in_polygon_test(3,points, point2d(x,y)))
  {
    cout << "2\t";
    return true;
  }

  points[0] = point2d(middle_x, middle_y);
  points[1] = point2d(middle_x + side_length/2, middle_y + sqrt(3)/2 * side_length);
  points[2] = point2d(middle_x + side_length, middle_y);

  if (jordan_point_in_polygon_test(3,points, point2d(x,y)))
  {
    cout << "3\t";
    return true;
  }

  points[0] = point2d(middle_x, middle_y);
  points[1] = point2d(middle_x + side_length, middle_y);
  points[2] = point2d(middle_x + side_length/2, middle_y - sqrt(3)/2 * side_length);

  if (jordan_point_in_polygon_test(3,points, point2d(x,y)))
  {
    cout << "4\t";
    return true;
  }

  points[0] = point2d(middle_x, middle_y);
  points[1] = point2d(middle_x - side_length/2, middle_y - sqrt(3)/2 * side_length);
  points[2] = point2d(middle_x + side_length/2, middle_y - sqrt(3)/2 * side_length);

  if (jordan_point_in_polygon_test(3,points, point2d(x,y)))
  {
    cout << "5\t";
    return true;
  }

  points[0] = point2d(middle_x, middle_y);
  points[1] = point2d(middle_x - side_length/2, middle_y - sqrt(3)/2 * side_length);
  points[2] = point2d(middle_x - side_length, middle_y );

  if (jordan_point_in_polygon_test(3,points, point2d(x,y)))
  {
    cout << "6\t";
    return true;
  }

  return false;
}

bool point_is_in_triangle(double x, double y, double middle_x, double middle_y, double rotation_angle, double side_length)
{
    //at first: calculate relative coordinates to the midle of the triangle, then rotate them
    point2d a_raw = point2d(0.0,sqrt(3)/2.0 * side_length/2.0);
    point2d b_raw = point2d(-side_length/2, -sqrt(3)/2.0 * side_length/2.0);
    point2d c_raw = point2d(side_length/2, -sqrt(3)/2.0 * side_length/2.0);

    //now turn all points
    point2d a = point2d(middle_x + cos(rotation_angle) * a_raw.x + sin(rotation_angle) * a_raw.y, middle_y -sin(rotation_angle) * a_raw.x + cos(rotation_angle) * a_raw.y);
    point2d b = point2d(middle_x + cos(rotation_angle) * b_raw.x + sin(rotation_angle) * b_raw.y, middle_y -sin(rotation_angle) * b_raw.x + cos(rotation_angle) * b_raw.y);
    point2d c = point2d(middle_x + cos(rotation_angle) * c_raw.x + sin(rotation_angle) * c_raw.y, middle_y -sin(rotation_angle) * c_raw.x + cos(rotation_angle) * c_raw.y);

    point2d points[] = {a,b,c};

    int result = jordan_point_in_polygon_test(3, points, point2d(x,y));
    return result >= 0;
}

bool point_has_fixed_border_constraint(int x, int y, int x_dimension, int y_dimension, double delta_x, double delta_y)
{
  //outer border values
  if (x == 0 || x == x_dimension-1 || y == 0 || y == y_dimension-1 )
  {
    return true;
  }

  //get electrode position in the grid and check for collision
  if (point_is_in_hexagon(x * delta_x, y*delta_y, x_position_electrode_2, y_position_electrode_2, size_electrode_2))
  {
    return true;
  }

  if (point_is_in_triangle(x * delta_x, y*delta_y, x_position_electrode_2, y_position_electrode_2, 0.0, size_electrode_2))
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

void sor_step(double** x_values, double* b_values, int x_dimension, int y_dimension, double delta_x, double delta_y, double s)
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
      x_values[x][y] = (1-s) * x_values[x][y] + s*((-b_values[x+y*x_dimension] + (x_values[x+1][y] + x_values[x-1][y])/pow(delta_x,2) + (x_values[x][y+1] + x_values[x][y-1])/pow(delta_y,2)) / ( 2* (1/pow(delta_x,2) + 1/pow(delta_y,2))));
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
      if (point_has_fixed_border_constraint(x,y,x_dimension, y_dimension, delta_x, delta_y))
      {
        continue;
      }

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

  double s = atof(argv[5]);

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
      if (point_is_in_hexagon(x * delta_x, y*delta_y, x_position_electrode_2, y_position_electrode_2, size_electrode_2))
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
      if (point_is_in_triangle(x * delta_x, y*delta_y, x_position_electrode_2, y_position_electrode_2, 0.0, size_electrode_2))
      {
        values[x][y] = v_electrode_2;
      }
    }
  }

  int totalIterations = 0;
  int lastCoutIteration = 0;
  double err = calculate_error(values,  b, dimension_x, dimension_y, delta_x, delta_y);
  while (err > epsilon)
  {
    if (totalIterations - lastCoutIteration > 25)
    {
      cout <<"iteration: " << totalIterations  << "\terror: " << err << endl;
      lastCoutIteration = totalIterations;
    }
    //gauss_seidel_step(values,  b, dimension_x, dimension_y, delta_x, delta_y);
    sor_step(values,  b, dimension_x, dimension_y, delta_x, delta_y,s);
    err = calculate_error(values,  b, dimension_x, dimension_y, delta_x, delta_y);

    totalIterations++;
  }

  //calc E
  for (int x = 0; x<dimension_x-1; x++)
  {
    for (int y = 0; y<dimension_y-1; y++)
    {
      //outputFile << x * delta_x << "\t" << y * delta_y << "\t" << (values[x+1][y] - values[x][y]) << "\t" << (values[x][y+1] - values[x][y]) << endl;
      outputFile << x << "\t" << y << "\t" << values[x][y] << endl;
    }
  }

  cout <<"final iteration: " << totalIterations  << "\tfinal error: " << err << endl;
}
