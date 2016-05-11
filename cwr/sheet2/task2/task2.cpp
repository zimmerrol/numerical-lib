#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

double step(double mu, double x)
{
  return 4*mu*x*(1-x);
}

void calculateSeriesAndWriteData(double x, double mu, int count, char* path)
{
  ofstream outputFile;
  outputFile.open(path);

  for (int i=0;i<=count; i++)
  {
    outputFile << i << "\t" << x << "\n";
    x = step(mu, x);
  }

  outputFile.flush();
  outputFile.close();
}

void calculateFixPointsAndWriteData(double x, double muStep, int count, char* path)
{
  ofstream outputFile;
  outputFile.open(path);

  double xTmp = x;

  for (double mu = 0;mu<=1;mu+=muStep)
  {
    xTmp = x;
    for (int i=0;i<=count; i++)
    {
      xTmp = step(mu, xTmp);
    }

    //no write all fixpoints
    for (int i=0;i<=count; i++)
    {
      outputFile << mu << "\t" << xTmp << "\n";
      xTmp = step(mu, xTmp);
    }
  }

  outputFile.flush();
  outputFile.close();
}

int main(int argc, char* argv[])
{
  if (argc < 2)
  {
    cout << "Wrong parameter count." << "\n";
    return -1;
  }

  int mode = atoi(argv[1]);
  if (mode == 0)
  {
    if (argc != 6)
    {
      cout << "Wrong parameter count." << "\n";
      return -1;
    }
    else
    {
      double x = atof(argv[2]);
      double mu = atof(argv[3]);
      int count = atoi(argv[4]);

      calculateSeriesAndWriteData(x, mu, count, argv[5]);
    }
  }
  else if(mode == 1)
  {
    if (argc != 6)
    {
      cout << "Wrong parameter count." << "\n";
      return -1;
    }
    else
    {
      double x = atof(argv[2]);
      double muStep = atof(argv[3]);
      int count = atoi(argv[4]);

      calculateFixPointsAndWriteData(x, muStep, count, argv[5]);
    }
  }




}
/*ofstream outputFile;
outputFile.open(argv[1]);
outputFile.flush();
outputFile.close();*/
