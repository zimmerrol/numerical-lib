#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "../../../lib/ode.h"

const double R_EARTH = 100;
const double R_MOON = 1000;
const double R_S = 1.06;
const double OMEGA = 2 * M_PI / (27.332 * 24);
const double ALPHA = 0.1 * M_PI;

class vector2d
{
  public:
    double x,y;
    vector2d();
    vector2d(double x, double y);
};

vector2d::vector2d()
{
  this->x = 0;
  this->y = 0;
}

vector2d::vector2d(double x, double y)
{
  this->x = x;
  this->y = y;
}

using namespace std;

const double v_sound = 343.0;
const double m = 3000.0;

vector2d gravity_pro_mass(vector2d target_position, vector2d source_position, double source_gravity_const)
{
  double squared_distance = pow(target_position.x - source_position.x,2)+pow(target_position.y - source_position.y,2);
  double distance = sqrt(squared_distance);

  vector2d result = vector2d(source_gravity_const / squared_distance * (source_position.x - target_position.x) / distance,
    source_gravity_const / squared_distance * (source_position.y - target_position.y) / distance);
  return result;
}

double dv_x(double* args, double* params)
{
  vector2d satelite = vector2d(args[0], args[1]);
  vector2d earth = vector2d(R_EARTH * cos(OMEGA * params[0]),R_EARTH * sin(OMEGA * params[0]));
  vector2d moon = vector2d(R_MOON * cos(OMEGA * params[0]),R_MOON * sin(OMEGA * params[0]));

  return gravity_pro_mass(satelite, earth, 20.0).x + gravity_pro_mass(satelite, moon, 20.0/81.3).x;
}

double dv_y(double* args, double* params)
{
  vector2d satelite = vector2d(args[0], args[1]);
  vector2d earth = vector2d(R_EARTH * cos(OMEGA * params[0]),R_EARTH * sin(OMEGA * params[0]));
  vector2d moon = vector2d(R_MOON * cos(OMEGA * params[0]),R_MOON * sin(OMEGA * params[0]));

  return gravity_pro_mass(satelite, earth, 20.0).y + gravity_pro_mass(satelite, moon, 20.0/81.3).y;
}

double dr_x(double* args, double* params)
{
  return args[2];
}

double dr_y(double* args, double* params)
{
  return args[3];
}

int main(int argc, char* argv[])
{
  //vector2d satelite_start = vector2d(R_S * cos(ALPHA) + R_EARTH * cos(0),R_S * sin(ALPHA) + R_EARTH * sin(0));

  for (double t = 0;t < 10;t+=0.01)
  {
    
  }
}
