#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "../../../lib/ode.h"

using namespace std;

const double R_EARTH = 100;
const double R_MOON = 1000;
const double R_S = 1.06;
const double OMEGA = 2 * atan(1) * 4 / (27.332 * 24);
double ALPHA;
double phase_shift;

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

vector2d gravity_pro_mass(vector2d target_position, vector2d source_position, double source_gravity_const)
{
  double squared_distance = pow(target_position.x - source_position.x,2)+pow(target_position.y - source_position.y,2);
  double distance = sqrt(squared_distance);

  vector2d result = vector2d(
		source_gravity_const / squared_distance * (source_position.x - target_position.x) / distance,
		source_gravity_const / squared_distance * (source_position.y - target_position.y) / distance
	  );

 // cout << "dist=" << distance << "\t" << (source_position.x - target_position.x)/ distance << "\t" << (source_position.y - target_position.y) / distance <<   "\n";
 
  return result;
}

double dv_x(double* args, double* params)
{
	double t = params[0];
  vector2d satelite = vector2d(args[0], args[1]);
  vector2d earth = vector2d(R_EARTH * cos(OMEGA * t), R_EARTH * sin(OMEGA * t));
  vector2d moon = vector2d(R_MOON * cos(OMEGA * t + phase_shift*atan(1) * 4), R_MOON * sin(OMEGA * t + phase_shift*atan(1) * 4));

  double acc = gravity_pro_mass(satelite, earth, 20.0).x + gravity_pro_mass(satelite, moon, 20.0/81.3).x;

  //cout << "x''=" << acc << "\n";
  return acc;
}

double dv_y(double* args, double* params)
{
	double t = params[0];
  vector2d satelite = vector2d(args[0], args[1]);
  vector2d earth = vector2d(R_EARTH * cos(OMEGA * t), R_EARTH * sin(OMEGA * t));
  vector2d moon = vector2d(R_MOON * cos(OMEGA * t + phase_shift*atan(1) * 4), R_MOON * sin(OMEGA * t + phase_shift*atan(1) * 4));

  double acc = gravity_pro_mass(satelite, earth, 20.0).y +gravity_pro_mass(satelite, moon, 20.0 / 81.3).y;
 // cout << "y''="<<  acc << "\n";
  return acc;

  
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
	ofstream outputFile;
	outputFile.open(argv[1]);
	outputFile << fixed << setprecision(5);

	double deltaT = atof(argv[2]);
	double t_max = atof(argv[3]);
	double v_0 = atof(argv[4]);

	ALPHA = atan(1) * 4 * atof(argv[5]);
	phase_shift = atof(argv[6]);

	double  v_init = sqrt(20.0 / R_S);


	double* explicitValues = new double[4];
	explicitValues[0] = R_S * cos(ALPHA) + R_EARTH * cos(0);
	explicitValues[1] = R_S * sin(ALPHA) + R_EARTH * sin(0);
	explicitValues[2] = -(v_0 + v_init)* sin(ALPHA);
	explicitValues[3] = (v_0 + v_init)* cos(ALPHA);

	numerical::odeFunction* functions = new numerical::odeFunction[4];
	functions[0] = &dr_x;
	functions[1] = &dr_y;
	functions[2] = &dv_x;
	functions[3] = &dv_y;

	double* params = new double[1];
	double delta_t_save = 1;
	double last_t = -2;
	for (double t = 0; t < t_max; t += deltaT)
	{
		params[0] = t;//;
		//params[1] = R_EARTH * sin(OMEGA * t);
		//params[2] = R_MOON * cos(OMEGA * t);
		//params[3] = R_MOON * sin(OMEGA * t);

		if (sqrt(pow(explicitValues[0] - R_MOON * cos(OMEGA * t + phase_shift*atan(1) * 4),2) + pow(explicitValues[1] - R_MOON * sin(OMEGA * t + phase_shift*atan(1) * 4), 2)) < 0.5)
		{
			cout << "The satelite hit the moon at t = " << t << " hours!\n";
			return 0;
		}

		if (t - last_t >= 0.1)
		{
			cout << "time = " << t << "\tdistance = " << sqrt(pow(explicitValues[0] - R_MOON * cos(OMEGA * t + phase_shift*atan(1) * 4), 2) + pow(explicitValues[1] - R_MOON * sin(OMEGA * t + phase_shift*atan(1) * 4), 2)) <<  "\n";
			outputFile << t << "\t";
			outputFile << explicitValues[0] << "\t" << explicitValues[1] << "\t";
			outputFile << explicitValues[2] << "\t" << explicitValues[3] << "\t";
			outputFile << R_EARTH * cos(OMEGA * t) << "\t" << R_EARTH * sin(OMEGA * t) << "\t";
			outputFile << R_MOON * cos(OMEGA * t + phase_shift*atan(1) * 4) << "\t" << R_MOON * sin(OMEGA * t + phase_shift*atan(1) * 4) << "\n";
			last_t = t;
		}


		numerical::stepRK4Explicit(4, functions, deltaT, explicitValues, params);
	}
}

//hit for: task02.exe res/out.dat 0.001 733.5 1.95 0.4 0.4677