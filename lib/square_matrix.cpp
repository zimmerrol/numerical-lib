#include "square_matrix.h"

namespace numerical
{
  //constructors: see header file

  double square_matrix::determinant()
  {

  }

  std::ostream& operator<<(std::ostream& stream, const square_matrix& matrix)
  {
      matrix.print_to_stream(stream);
      return stream;
   }
}
