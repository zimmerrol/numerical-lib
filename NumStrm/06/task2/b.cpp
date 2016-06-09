#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

void integrate_FTCS(vector<double> &grid, double &kappa,
  double &v,double &dx, double &dt)
{
  size_t N = grid.size();
  vector<double> new_grid(N) ;
  new_grid[N-1] = 1;
  new_grid[0] = 0.0;
  for(size_t i = 1; i < N-1; i++)
  {
      new_grid[i] = grid[i] + dt/dx * ((grid[i-1]-grid[i+1])/2.0*v
                  + kappa*(grid[i+1] - 2*grid[i] + grid[i-1])/dx);
  }
  grid = new_grid;
}

void make_plot_file(vector<double> & grid, double dx)
{
  ofstream out; out.open("data.dat");
  for(size_t i = 0; i < grid.size(); i++)
  {
    out << i*dx << "\t" << grid[i] << endl;
  }
  out.close();
}

int tridiag(int n, double * a,double * b, double * c, double * d, vector<double> &T)
{

  c[1] /= b[1];
  d[1] /= b[1];

  for(size_t i = 2; i < n; i++)
  {
    c[i] /= (b[i] - a[i] * c[i-1] );
  //  if(isnan(d[i-1])) cout << "d[" << i-1 <<"] is nan" << endl;
    d[i] = (d[i] - a[i] * d[i-1])/(b[i] - a[i] * c[i-1]);
  //  if(isnan(d[i])) cout << "d[" << i <<"] is nan" << endl;
    //cout << d[i] << "\t" << c[i] << endl;
  }

//  T[n] = 1;

  for(size_t i = n-1; i>= 1 ; i--)
  {
    T[i] = d[i] - c[i]*T[i+1];
    /*if( isnan(T[i]) )
      {
          cout << "nan at " << i << endl;
      }*/
  }

  return 42;
}



void integrate_exim(vector<double> &grid, double &kappa, //######################################
  double &v,double &dx, double &dt)
{

  int N = grid.size();

  double* a = new double[N+1];
  double* b = new double[N+1];
  double* c = new double[N+1];
  double* d = new double[N+1];

  double alpha = 1/dt + 2*kappa/dx/dx;
  double beta = -kappa/dx/dx;
  //stay with the convention to label a, b, c starting from 1 to avoid confusion
  // with subscribts of grid points

  for(size_t i = 2; i < N-1; i++)
  {
    b[i] = alpha;
    a[i] = beta; c[i] = beta;
    d[i] = v/2/dx*(grid.at(i-1) - grid.at(i+1)) + grid.at(i)/dt;
  }

  b[1] = alpha;
  d[1] = v/2/dx*(grid[0] - grid[2]) + grid[1]/dt ;
  //cout << d[1] << endl;

  d[N-1] = v/2/dx*(grid.at(N-2) - grid[N]) + kappa*grid.at(N)/dx/dx + grid.at(N-1)/dt;

  //cout << "d[N] = " << d[N]<< endl;

  c[1] = beta; a[N-1] = beta;

  tridiag(N-1, a, b, c, d, grid);
  delete a; delete b; delete c; delete d;
}


int main()
{
  size_t N = 200;
  double D = 0.1;
  double v = 1.0;
  double dt = 1e-2;
  double dx = 1.0/(N-1);
  vector<double> grid;

  for(size_t i = 0; i < N ; i++)
  {
      grid.push_back( (double) i/N);
  }

  for(double t = 0; t<20 ; t+=dt)
  {
    integrate_exim(grid, D, v, dx, dt);
  }

  make_plot_file(grid,dx);
  return 0;
}
