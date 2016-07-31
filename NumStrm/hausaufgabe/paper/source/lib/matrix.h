#include <fstream>

#ifndef MATRIX_H
#define MATRIX_H
namespace numerical{
  class matrix
  {
    public:
      //creates a new nXm matrix
      matrix(size_t n, size_t m);
      //creates a new nXm matrix and initializes all values with defaultValue
      matrix(size_t n, size_t m, double defaultValue);
      //deconstructor
      ~matrix();

      //returns the number of rows/columns
      size_t get_n();
      size_t get_m();

      //returns the value in row n, column m
      double get_value(size_t n, size_t m);
      //sets the value in row n, column m
      void set_value(size_t n, size_t m, double value);
    private:
      double** data;
      size_t n;
      size_t m;
  };
}
#endif
