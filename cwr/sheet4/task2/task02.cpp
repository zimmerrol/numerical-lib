#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "../../../lib/ode.h"

using namespace std;

const double v_sound = 343.0;
const double m = 3000.0;

double f_c_w(double* args, double* params)
{
  double M = args[1] / v_sound;
  return 0.227+0.120 * (pow(M,2) - 1.18 * M + 0.205)/( pow(M,3) - 3 * M + 2.05);
}

double f_F_friction(double* args, double* params)
{
  return pow(args[1],2) * f_c_w(args, params);
}

double f_r(double* args, double* params)
{
  return args[1];
}

double f_v(double* args, double* params)
{
  return args[2];
}

double f_a(double* args, double* params)
{
  double F = params[0];
  return 1/m * (F-f_F_friction(args, params));
}

void leapFrogStep(double &x, double &v, double deltaT, double f)
{
  double M = v/v_sound;
  double fr = (0.227+0.120 * (pow(M,2) - 1.18 * M + 0.205)/( pow(M,3) - 3 * M + 2.05)) * v *v;
  double a = (f-fr)/m;

  v += deltaT * a;
  x += v * deltaT;
}

void calcLFforFixedF(char* argv[])
{
  ofstream outputFile;
  outputFile.open(argv[2]);
  outputFile << fixed << setprecision(9);

  double deltaT = atof(argv[3]);

  double F = atof(argv[4]);

  double x = 0;
  double v = F/m * deltaT;

  for(double t=0; t<=500; t+=deltaT)
  {
    outputFile << t << "\t" << x << "\t" << v << "\n";//"\t" << xIm << "\t" << yIm  << "\n";
    leapFrogStep(x,v, deltaT, F);
  }
}

void calcLFforFixedF2(char* argv[])
{
  ofstream outputFile;
  outputFile.open(argv[2]);
  outputFile << fixed << setprecision(9);

  double deltaT = atof(argv[3]);

  cout << deltaT << "\n";

  double F = atof(argv[4]);
  double params[] = {F};

  numerical::odeFunction* functions = new numerical::odeFunction[2];
  functions[0] = &f_r;
  functions[1] = &f_a;

  double* previousValues = new double[2];
  previousValues[0] = 0;
  previousValues[1] = 0;
  double* values = new double[2];
  values[0] = 0;
  values[1] = F/m*deltaT;

  outputFile << 0 << "\t" << previousValues[0] << "\t" << previousValues[1] << "\n";
  for(double t=deltaT; t<=500; t+=deltaT)
  {
    outputFile << t << "\t" << values[0] << "\t" << values[1] << "\n";
    //numerical::stepLeapFrog(2, functions, deltaT, values, previousValues, params);
    double origX = values[0];
    double origV = values[1];
    double oldX = previousValues[0];
    double oldV = previousValues[1];

    double M = origV/v_sound;

    previousValues[0] = values[0];
    previousValues[1] = values[1];

    values[0] = oldX + 2*deltaT * origV;
    values[1] = oldV + 2*deltaT * ((F-(0.227+0.120 * (pow(M,2) - 1.18 * M + 0.205)/( pow(M,3) - 3 * M + 2.05)) * origV * origV)/m);
  }
}

void findMaxSpeed(char* argv[])
{
  ofstream outputFile;
  outputFile.open(argv[2]);
  outputFile << fixed << setprecision(9);

  double deltaT = atof(argv[3]);

  double F = atof(argv[4]);
  double params[] = {F};

  numerical::odeFunction* functions = new numerical::odeFunction[2];
  functions[0] = &f_r;
  functions[1] = &f_a;

  while (true)
  {
    double* rkValues = new double[2];
    rkValues[0] = 0;
    rkValues[1] = 0;

    for(double t=0; t<=500; t+=deltaT)
    {
      numerical::stepRK2Explicit(2, functions, deltaT, rkValues, params);
    }

    if (rkValues[1] > v_sound)
    {
      //calc again and save the results for the plot
      rkValues[0] = 0;
      rkValues[1] = 0;

      for(double t=0; t<=500; t+=deltaT)
      {
        outputFile << t << "\t" << rkValues[0] << "\t" << rkValues[1] << "\t" << rkValues[2] << "\n";//"\t" << xIm << "\t" << yIm  << "\n";
        numerical::stepRK2Explicit(2, functions, deltaT, rkValues, params);
      }

      break;
    }
    else
    {
      F += 10;
      cout << "trying: " << F << "\n";
    }
  }
}


void calcRK2forFixedF(char* argv[])
{
  ofstream outputFile;
  outputFile.open(argv[2]);
  outputFile << fixed << setprecision(9);

  double deltaT = atof(argv[3]);

  double F = atof(argv[4]);
  double params[] = {F};

  numerical::odeFunction* functions = new numerical::odeFunction[2];
  functions[0] = &f_r;
  functions[1] = &f_a;

  double* rkValues = new double[2];
  rkValues[0] = 0;
  rkValues[1] = 0;

  for(double t=0; t<=500; t+=deltaT)
  {
    outputFile << t << "\t" << rkValues[0] << "\t" << rkValues[1] << "\t" << rkValues[2] << "\n";//"\t" << xIm << "\t" << yIm  << "\n";
    numerical::stepRK2Explicit(2, functions, deltaT, rkValues, params);
  }
}

void calcMaxSpeed(char* argv[])
{
  ofstream outputFile;
  outputFile.open(argv[2]);
  outputFile << fixed << setprecision(9);

  double deltaT = atof(argv[3]);

  double F = atof(argv[4]);
  double params[] = {F};
  double F_final = atof(argv[5]);

  numerical::odeFunction* functions = new numerical::odeFunction[2];
  functions[0] = &f_r;
  functions[1] = &f_a;

  while (F < F_final)
  {
    double* rkValues = new double[2];
    rkValues[0] = 0;
    rkValues[1] = 0;

    for(double t=0; t<=500; t+=deltaT)
    {
      numerical::stepRK2Explicit(2, functions, deltaT, rkValues, params);
    }

    outputFile << F << "\t" << rkValues[1] << "\n";

    F += 10;
  }
}

int main(int argc, char* argv[])
{
  int mode = atoi(argv[1]);
  if (mode == 0)
  {
    calcRK2forFixedF(argv);
  }
  else if (mode == 1)
  {
    findMaxSpeed(argv);
  }
  else if (mode == 2)
  {
    calcMaxSpeed(argv);
  }
  else if (mode == 3)
  {
    calcLFforFixedF2(argv);
  }
  else{
      calcLFforFixedF(argv);
  }
}
