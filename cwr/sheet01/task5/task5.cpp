#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <cstddef>
#include <iomanip>

using namespace std;

double faculty(int n)
{
  if (n == 0)
  {
    return 1;
  }
  return n*faculty(n-1);
}

double csin(double x)
{
  double epsilon = std::numeric_limits<double>::min();

  double result = 0.0;
  double change = 1.0;
  int k = 0;

  while (change > epsilon || change < -epsilon)
  {
    change = pow(-1,k)*pow(x,1+2*k)/faculty(1+2*k);
    result += change;
    k++;
  }

  return result;
}

int main(int argc, char* argv[])
{
  if (argc != 4)
  {
    cout << "No output/range specified!" << endl;
    return -1;
  }

  ofstream outputFile;
  outputFile.open(argv[1]);

  outputFile << fixed << setprecision(9);

  for (double x = atof(argv[2]);x<atof(argv[3]);x+=0.001)
  {
    outputFile << x << "\t" << csin(x) << endl;
  }
  outputFile << atof(argv[3]) << "\t" << csin(atof(argv[3])) << endl;

  outputFile.flush();
  outputFile.close();

  cout << "Data written to: " << argv[1] << endl;
}
