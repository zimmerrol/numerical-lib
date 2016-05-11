#include "matrix.h"
#include <iostream>

namespace numerical{
  Matrix::Matrix(int m, int n)
  {
    this->n = n;
    this->m = m;
    data = new double*[m];
    for (int i=0;i<m;i++)
    {
      data[i] = new double[n];
    }
  }

  Matrix::Matrix(int m, int n, double defaultValue)
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

  Matrix::~Matrix()
  {
    delete data;
  }

  void Matrix::swapRows(int targetm, int sourcem)
  {
    double buffer;
    for (int i=0;i<this->n;i++)
    {
      buffer = this->data[targetm][i];
      this->data[targetm][i] = this->data[sourcem][i];
      this->data[sourcem][i] = buffer;
    }
  }

  void Matrix::replaceRow(int targetm, double* source)
  {
    for (int i=0;i<this->n;i++)
    {
      this->data[targetm][i] = source[i];
    }
  }

  void Matrix::swapColumn(int targetn, int sourcen)
  {
    double buffer;
    for (int i=0;i<this->n;i++)
    {
      buffer = this->data[i][targetn];
      this->data[i][targetn] = this->data[i][sourcen];
      this->data[i][sourcen] = buffer;
    }
  }

  void Matrix::replaceColumn(int targetn, double* source)
  {
    for (int i=0;i<this->n;i++)
    {
      this->data[i][targetn] = source[i];
    }
  }

  void Matrix::printToCout()
  {
    for (int y = 0;y<m;y++)
    {
      for (int x = 0;x<n;x++)
      {
        std::cout << this->data[x][y] << "\t";
      }
      std::cout << "\n";
    }
  }
}
