#include "matrix.h"
#include <iostream>

namespace numerical{
  //creates a new nXm matrix
  matrix::matrix(size_t n, size_t m)
  {
    this->n = n;
    this->m = m;
    data = new double*[m];
    for (size_t i=0;i<m;i++)
    {
      data[i] = new double[n];
    }
  }

  //creates a new nXm matrix and initializes all values with defaultValue
  matrix::matrix(size_t n, size_t m, double defaultValue)
  {
    this->n = n;
    this->m = m;
    data = new double*[m];
    for (size_t i=0;i<m;i++)
    {
      data[i] = new double[n];
      for (size_t j=0;j<n;j++)
      {
        data[i][j] = defaultValue;
      }
    }
  }

  //deconstructor
  matrix::~matrix()
  {
    delete data;
  }

  //sets the value in the nth row, mth column
  void matrix::set_value(size_t n, size_t m, double value)
  {
    data[n][m] = value;
  }

  //returns the value in the nth row, mth column
  double matrix::get_value(size_t n, size_t m)
  {
    return data[n][m];
  }

  //returns the number of rows
  size_t matrix::get_n()
  {
      return this->n;
  }

  //returns the number of columns
  size_t matrix::get_m()
  {
    return this->m;
  }
}
