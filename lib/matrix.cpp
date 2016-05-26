#include "matrix.h"
#include <iostream>

namespace numerical{
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

  matrix::~matrix()
  {
    delete data;
  }

  // n=row, m=column
  //data = [rowIndex][columnIndex]
  void matrix::set_value(size_t n, size_t m, double value)
  {
    data[n][m] = value;
  }

  // n=row, m=column
  //data = [rowIndex][columnIndex]
  double matrix::get_value(size_t n, size_t m)
  {
    return data[n][m];
  }

  size_t matrix::get_n()
  {
      return this->n;
  }

  size_t matrix::get_m()
  {
    return this->m;
  }

  void matrix::multiply_row(size_t target_m, double scalar)
  {
    for (size_t i=0; i<this->m; i++)
    {
      this->data[target_m][i] *= scalar;
    }
  }

  void matrix::multiply_column(size_t target_n, double scalar)
  {
    for (size_t i=0; i<this->n; i++)
    {
      this->data[i][target_n] *= scalar;
    }
  }

  void matrix::swap_rows(size_t target_m, size_t source_m)
  {
    double buffer;
    for (size_t i=0;i<this->n;i++)
    {
      buffer = this->data[target_m][i];
      this->data[target_m][i] = this->data[source_m][i];
      this->data[source_m][i] = buffer;
    }
  }

  void matrix::replace_row(size_t target_m, double* source)
  {
    for (size_t i=0;i<this->n;i++)
    {
      this->data[target_m][i] = source[i];
    }
  }

  void matrix::swap_column(size_t target_n, size_t source_n)
  {
    double buffer;
    for (size_t i=0;i<this->n;i++)
    {
      buffer = this->data[i][target_n];
      this->data[i][target_n] = this->data[i][source_n];
      this->data[i][source_n] = buffer;
    }
  }

  void matrix::replace_column(size_t target_n, double* source)
  {
    for (size_t i=0;i<this->n;i++)
    {
      this->data[i][target_n] = source[i];
    }
  }

  void matrix::print_to_stream(std::ostream& stream) const
  {
    for (size_t y = 0;y<n;y++)
    {
      for (size_t x = 0;x<n;x++)
      {
        stream << this->data[y][x] << "\t";
      }
      stream << std::endl;
    }
    stream << std::flush;
  }

  std::ostream& operator<<(std::ostream& stream, const matrix& matrix)
  {
      matrix.print_to_stream(stream);
      return stream;
   }
}
