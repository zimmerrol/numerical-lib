#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

using namespace std;

typedef vector<vector<double> > grid_t;

double calculate_non_linear_term(grid_t &grid, const double x, const double y, const double delta_x, const double delta_y)
{
  return grid.at(x).at(y) * (grid.at(x+1).at(y) - grid.at(x-1).at(y))/(2.0*delta_x) +
            grid.at(x).at(y) * (grid.at(x).at(y+1) - grid.at(x).at(y-1))/(2.0*delta_y);
}

double calculate_visc_term(grid_t &grid, const double x, const double y, const double delta_x, const double delta_y, const double re)
{
  return 1.0/re * ((grid.at(x+1).at(y) + grid.at(x-1).at(y) - 2.0*grid.at(x).at(y))/(pow(delta_x,2)) +
                          (grid.at(x).at(y+1) + grid.at(x).at(y-1) - 2.0*grid.at(x).at(y))/(pow(delta_y,2)));
}

void integrate(grid_t &velocity_x, grid_t &velocity_y, grid_t &old_velocity_x, grid_t &old_velocity_y, grid_t &pressure, const double delta_x, const double delta_y, const double re, const double delta_t)
{
  const size_t dimension_x = velocity_x.size();
  const size_t dimension_y = velocity_y.size();

  double non_linear_term_x, non_linear_term_y;
  double visc_term_x, visc_term_y;

  double old_non_linear_term_x, old_non_linear_term_y;
  double old_visc_term_x, old_visc_term_y;

  grid_t new_velocity_x, new_velocity_y;

  double velocity_star_x, velocity_star_y;

  for (size_t x=0; x<dimension_x; x++)
  {
    new_velocity_x.push_back(vector<double>(dimension_y));

    for (size_t y=0; y<dimension_y; y++)
    {
      /*non_linear_term_x = velocity_x.at(x).at(y) * (velocity_x.at(x+1).at(y) - velocity_x.at(x-1).at(y))/(2.0*delta_x) +
                            velocity_y.at(x).at(y) * (velocity_x.at(x).at(y+1) - velocity_x.at(x).at(y-1))/(2.0*delta_y);
      non_linear_term_y = velocity_x.at(x).at(y) * (velocity_y.at(x+1).at(y) - velocity_y.at(x-1).at(y))/(2.0*delta_x) +
                            velocity_y.at(x).at(y) * (velocity_y.at(x).at(y+1) - velocity_y.at(x).at(y-1))/(2.0*delta_y);*/
      non_linear_term_x = calculate_non_linear_term(velocity_x, x, y, delta_x, delta_y);
      non_linear_term_y = calculate_non_linear_term(velocity_y, x, y, delta_x, delta_y);
      visc_term_x = calculate_visc_term(velocity_x, x, y, delta_x, delta_y, re);
      visc_term_y = calculate_visc_term(velocity_y, x, y, delta_x, delta_y, re);

      old_non_linear_term_x = calculate_non_linear_term(old_velocity_x, x, y, delta_x, delta_y);
      old_non_linear_term_y = calculate_non_linear_term(old_velocity_y, x, y, delta_x, delta_y);
      old_visc_term_x = calculate_visc_term(old_velocity_x, x, y, delta_x, delta_y, re);
      old_visc_term_y = calculate_visc_term(old_velocity_y, x, y, delta_x, delta_y, re);

      /*visc_term_x = 1.0/re * ((velocity_x.at(x+1).at(y) + velocity_x.at(x-1).at(y) - 2.0*velocity_x.at(x).at(y))/(pow(delta_x,2)) +
                              (velocity_x.at(x).at(y+1) + velocity_x.at(x).at(y-1) - 2.0*velocity_x.at(x).at(y))/(pow(delta_y,2)));

      visc_term_y = 1.0/re * ((velocity_y.at(x+1).at(y) + velocity_y.at(x-1).at(y) - 2.0*velocity_y.at(x).at(y))/(pow(delta_x,2)) +
                              (velocity_y.at(x).at(y+1) + velocity_y.at(x).at(y-1) - 2.0*velocity_y.at(x).at(y))/(pow(delta_y,2)));*/


      velocity_star_x = velocity_x.at(x).at(y)  + delta_t * 1.0/2.0*(3.0*(-non_linear_term_x+visc_term_x) - (-old_non_linear_term_x+old_visc_term_x));
      velocity_star_y = velocity_y.at(x).at(y)  + delta_t * 1.0/2.0*(3.0*(-non_linear_term_y+visc_term_y) - (-old_non_linear_term_y+old_visc_term_y));

      //now solve laplace(p^n+1) = div(v*)/delta_t
      //gauss seidel


      //calculate now grad(p)
      //calculate then v^n+1 = u*-delta_t*grad(p)

      //save the calculated value in the new_velocity_x/y vectors
    }
  }

  //copy velocity_x/y to old_velocity_x/y
  //copy new_velocity_x/y to velocity_x/y
  for (size_t x=0; x<dimension_x; x++)
  {
    for (size_t y=0; y<dimension_y; y++)
    {
      old_velocity_x.at(x).at(y) = velocity_x.at(x).at(y);
      old_velocity_y.at(x).at(y) = velocity_y.at(x).at(y);

      velocity_x.at(x).at(y) = new_velocity_x.at(x).at(y);
      velocity_y.at(x).at(y) = new_velocity_y.at(x).at(y);
    }
  }

}

int main(int argc, char* argv[])
{
  cout << "Expects:" << endl;
  cout << "\t" << "outputPath\tdimension_x\tdimension_y\tmax_time\tdelta_t\tRe" << endl;
  ofstream outputFile;
	outputFile.open("res/out.dat");
	outputFile << fixed << setprecision(5);

  const double width = 1.0;
  const double height = 1.0;
  const double delta_x = 1.0/(atof(argv[2]));
  const double delta_y = 1.0/(atof(argv[3]));
  const float max_time = atof(argv[4]);
  const double delta_t = atof(argv[5]);
  const double Re = atof(argv[6]);

  grid_t velocity_x, velocity_y, old_velocity_x, old_velocity_y, pressure;

  for (double t=0; t<max_time; t+=delta_t)
  {
    integrate(velocity_x, velocity_y, old_velocity_x, old_velocity_y, pressure, delta_x, delta_y, Re, delta_t);
  }

  return 0;
}
