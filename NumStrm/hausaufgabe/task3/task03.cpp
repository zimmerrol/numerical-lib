#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include "../../../lib/ode.h"

using namespace std;

typedef vector<vector<double > > grid_t;


//returns the error
double ftcs_time_step(grid_t& values, grid_t v0x, grid_t v0y, grid_t source, const double delta_t, const double delta_x, const double delta_y, const double pe)
{
  size_t dimension_x = values.size();
  size_t dimension_y = values.at(0).size();

  grid_t old_values;
  for (size_t x=0; x<dimension_x; x++)
  {
    old_values.push_back(vector<double>());
    for (size_t y=0; y<dimension_y; y++)
    {
      old_values.at(x).push_back(values.at(x).at(y));
    }
  }

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
                    + source.at(x).at(y)
                );
    }
  }

  //left/right boundaries
  for (size_t y=0; y<dimension_y; y++)
  {
    values.at(0).at(y) = 1.0/3.0*(4.0*values.at(0+1).at(y) - values.at(0+2).at(y));
    values.at(dimension_x-1).at(y) = -2.0/3.0*(-2.0*values.at(dimension_x-2).at(y) + 1.0/2.0*values.at(dimension_x-3).at(y));
  }

  //bottom/top boundaries
  for (size_t x=0; x<dimension_x; x++)
  {
    values.at(x).at(0) = 0.0;
    values.at(x).at(dimension_y-1) = 1.0;
  }

  double error = 0.0;
  for (size_t x=0; x<dimension_x; x++)
  {
    for (size_t y=0; y<dimension_y; y++)
    {
      error += pow(
                  values.at(x).at(y) - (cos(M_PI*delta_x*x)*sin(M_PI*delta_y*y)+y*delta_y)
                ,2);
    }
  }

  return sqrt(error);
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

  size_t dimension_x = (atoi)(argv[2])+1;
  double delta_x = 1.0/(dimension_x-1);

  size_t dimension_y = (atoi)(argv[3])+1;
  double delta_y = 1.0/(dimension_y-1);

  double pe = atof(argv[4]);

  double max_t = atof(argv[5]);
  double delta_t = atof(argv[6]);

  cout << dimension_x << "\t" << dimension_y << "\t" << pe << "\t" << max_t << "\t" << delta_t << endl;

  //ftcs_time_step(grid_t values, grid_t v0x, grid_t v0y, double delta_t, double delta_x, double delta_y, double pe)

  grid_t v0x;
  grid_t v0y;
  grid_t values;
  grid_t sources;
  for (size_t x=0; x<dimension_x; x++)
  {
    v0x.push_back(vector<double>());
    v0y.push_back(vector<double>());
    values.push_back(vector<double>());
    sources.push_back(vector<double>());

    for (size_t y=0; y<dimension_y; y++)
    {
      v0x.at(x).push_back(M_PI*sin(2*M_PI*x*delta_x)*cos(M_PI*y*delta_y));
      v0y.at(x).push_back(-2.0*M_PI*cos(2*M_PI*x*delta_x)*sin(M_PI*y*delta_y));
      values.at(x).push_back(y*delta_y);
      sources.at(x).push_back(
        pow(M_PI,2)*(cos(M_PI*delta_x*x)*sin(M_PI*delta_y*y)*2
        - pe*pow(cos(x*delta_x*M_PI),3)*sin(2*M_PI*delta_y*y))
        -pe*2*M_PI*cos(2*M_PI*delta_x*x)*sin(M_PI*y*delta_y)
      );
    }
  }

int aa = 0;
  for (double t=0.0; t<max_t; t+=delta_t)
  {
    double err = ftcs_time_step(values, v0x, v0y, sources, delta_t, delta_x, delta_y, pe) ;
    aa++;
    if (aa % 100 == 0)
    {
        cout << "time:\t" << t << "\terror:\t" << err << "\n" << flush;
        aa = 0;
    }
  }

  for (size_t y=0; y<dimension_y; y++)
  {
    for (size_t x=0; x<dimension_x; x++)
    {
      outputFile << sources.at(x).at(y) << "\t";//values.at(x).at(y) << "\t";//x*delta_x << "\t" << y*delta_y << "\t" << values.at(x).at(y) << endl;
    }
    outputFile << "\n";
  }

  return 1;
}
