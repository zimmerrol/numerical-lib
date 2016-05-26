#include "matrix.h"

#ifndef SQUAREMATRIX_H
#define SQUAREMATRIX_H
namespace numerical
{
  class square_matrix :  public matrix
  {
    public:
      square_matrix(int n) : matrix(n,n) {}
      square_matrix(int n, double defaultValue) : matrix(n,n, defaultValue) {};

      double determinant();
  };
}
#endif
