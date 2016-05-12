#include <iostream>

using namespace std;

double sumUp(int n)
{
  double result = 0.0;
  for (int i=0;i<=n;i++)
  {
    result += i*i*i;
  }
  return result;
}

double sumDown(int n)
{
  double result = 0.0;
  for (int i=n;i>=0;i--)
  {
    result += i*i*i;
  }
  return result;
}

int main(int argc, char* argv[])
{
 cout << sumUp(200) << "\t" << sumDown(200) << endl;
 cout << sumUp(200) - sumDown(200) << endl;
}
