#include <iostream>
using namespace std;

double calc(int n)
{
  return n*(1.0/3.0)-n/3.0;
}

int main(int argc, char* argv[])
{
  cout << fixed << calc(7) << ";" << calc(8) << endl << scientific;
  cout <<  calc(7) << ";" << calc(8) << endl;

}
