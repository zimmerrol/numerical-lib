#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <iomanip>
#include <cmath>

using namespace std;

//custom type which discribes a 2d grid. use this instead of arrays for access violation safety
typedef vector<vector<double > > grid_t;


//integrates the ODE with the FTCS scheme. expects the parameters:
//  values: a grid_t object which contains current distribution of the temperature T
//  v0x, v0y: a grid_t object which contains the x/y component of the velocity v0
//  delta_t, delta_x, delta_y, pe: the time-step-size, the step size in x/y direction and the Péclet number
void ftcs_time_step(grid_t& values, grid_t v0x, grid_t v0y, const double delta_t, const double delta_x, const double delta_y, const double pe)
{
  //get the dimensions of the grid
  size_t dimension_x = values.size();
  size_t dimension_y = values.at(0).size();

  //create copy of the grid with the same dimensions
  //this will be used for the calculation of the derivatives
  grid_t old_values;
  for (size_t x=0; x<dimension_x; x++)
  {
    //add a new row
    old_values.push_back(vector<double>());
    for (size_t y=0; y<dimension_y; y++)
    {
      //copy value
      old_values.at(x).push_back(values.at(x).at(y));
    }
  }

  //loop over the entire inner grid (ignore the 4 outer sides of the rectangle) and calculate the new value using the FTCS-scheme
  for (size_t x=1; x<dimension_x-1; x++)
  {
    for (size_t y=1; y<dimension_y-1; y++)
    {
      values.at(x).at(y) += delta_t*
                (
                  1.0/pow(delta_x,2)*(old_values.at(x+1).at(y) - 2.0*old_values.at(x).at(y) + old_values.at(x-1).at(y))
                + 1.0/pow(delta_y,2)*(old_values.at(x).at(y+1) - 2.0*old_values.at(x).at(y) + old_values.at(x).at(y-1))
                  -pe*
                    (
                      v0x.at(x).at(y)/(2.0*delta_x)*(old_values.at(x+1).at(y)-old_values.at(x-1).at(y))
                    + v0y.at(x).at(y)/(2.0*delta_y)*(old_values.at(x).at(y+1)-old_values.at(x).at(y-1))
                    )
                );
    }
  }

  //now take care of the four boundary conditions

  //left/right boundaries (Neumann boundaries)
  for (size_t y=0; y<dimension_y; y++)
  {
    values.at(0).at(y) = 1.0/3.0*(4.0*values.at(0+1).at(y) - values.at(0+2).at(y));
    values.at(dimension_x-1).at(y) = -2.0/3.0*(-2.0*values.at(dimension_x-2).at(y) + 1.0/2.0*values.at(dimension_x-3).at(y));
  }

  //bottom/top boundaries (Dirichlet boundaries)
  for (size_t x=0; x<dimension_x; x++)
  {
    values.at(x).at(0) = 0.0;
    values.at(x).at(dimension_y-1) = 1.0;
  }

}

int main(int argc, char* argv[])
{
  //print information for the usage
  cout << "Use this program to calculate the temperature with the FTCS-scheme in the rectangle." << endl;
  cout << "Expects the follwing set of parameters:" << endl;

  cout << "\toutput file:\tThe path to the file in which the results will be stored." << endl;
  cout << "\tN_x:\t\tThe amount of grid points in the x-direction." << endl;
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
  cout << "\tNx:" << dimension_x << "\tNy:" << dimension_y << "\tPe:" << pe << "\tMax t" << max_t << "\tDelta t:" << delta_t << endl;

  //to store the x component of the velocity v0
  grid_t v0x;
  //to store the y component of the velocity v0
  grid_t v0y;
  //to store the current temperature
  grid_t values;

  //create the new grids
  for (size_t x=0; x<dimension_x; x++)
  {
    //add a new row
    v0x.push_back(vector<double>());
    v0y.push_back(vector<double>());
    values.push_back(vector<double>());

    for (size_t y=0; y<dimension_y; y++)
    {
      //add the v0 values
      v0x.at(x).push_back(M_PI*sin(2*M_PI*x*delta_x)*cos(M_PI*y*delta_y));
      v0y.at(x).push_back(-2.0*M_PI*cos(2*M_PI*x*delta_x)*sin(M_PI*y*delta_y));

      //add the start temperature of the rectangle
      values.at(x).push_back(y*delta_y);
    }
  }

  //loop over the time, until we have reached the maximum time
  for (double t=0.0; t<max_t; t+=delta_t)
  {
    ftcs_time_step(values, v0x, v0y, delta_t, delta_x, delta_y, pe);
  }

  //print the final grid's values in the givven output file in a matrix form
  for (size_t y=0; y<dimension_y; y++)
  {
    for (size_t x=0; x<dimension_x; x++)
    {
      outputFile << values.at(x).at(y) << "\t";
    }
    outputFile << "\n";
  }

  return 1;
}
