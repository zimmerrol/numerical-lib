#include "matrix.h"
#include <iostream>

namespace numerical{
  matrix::matrix(int m, int n)
  {
    this->n = n;
    this->m = m;
    data = new double*[m];
    for (int i=0;i<m;i++)
    {
      data[i] = new double[n];
    }
  }

  matrix::matrix(int m, int n, double defaultValue)
  {
    this->n = n;
    this->m = m;
    data = new double*[m];
    for (int i=0;i<m;i++)
    {
      data[i] = new double[n];
      for (int j=0;j<n;j++)
      {
        data[i][j] = defaultValue;
      }
    }
  }

  matrix::~matrix()
  {
    delete data;
  }

  void matrix::swap_rows(int targetm, int sourcem)
  {
    double buffer;
    for (int i=0;i<this->n;i++)
    {
      buffer = this->data[targetm][i];
      this->data[targetm][i] = this->data[sourcem][i];
      this->data[sourcem][i] = buffer;
    }
  }

  void matrix::replace_row(int targetm, double* source)
  {
    for (int i=0;i<this->n;i++)
    {
      this->data[targetm][i] = source[i];
    }
  }

  void matrix::swap_column(int targetn, int sourcen)
  {
    double buffer;
    for (int i=0;i<this->n;i++)
    {
      buffer = this->data[i][targetn];
      this->data[i][targetn] = this->data[i][sourcen];
      this->data[i][sourcen] = buffer;
    }
  }

  void matrix::replace_column(int targetn, double* source)
  {
    for (int i=0;i<this->n;i++)
    {
      this->data[i][targetn] = source[i];
    }
  }

  void matrix::print_to_stream(std::ostream& stream) const
  {
    for (int y = 0;y<m;y++)
    {
      for (int x = 0;x<n;x++)
      {
        stream << this->data[x][y] << "\t";
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
