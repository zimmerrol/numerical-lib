#include <cmath>
#include "vector2d.h"

namespace numerical{
  vector2d::vector2d()
  {
    this->x = 0;
    this->y = 0;
  }

  vector2d::vector2d(double x, double y)
  {
    this->x = x;
    this->y = y;
  }

  vector2d::vector2d(const vector2d &value)
  {
    this->x = value.x;
    this->y = value.y;
  }

  vector2d::~vector2d()
  {

  }

  double vector2d::length()
  {
    return sqrt(x*x+y*y);
  }

  void vector2d::add(vector2d value)
  {
    x += value.x;
    y += value.y;
  }

  void vector2d::substract(vector2d value)
  {
    x -= value.x;
    y -= value.y;
  }

  void vector2d::scalar_multiply(double value)
  {
    x *= value;
    y *= value;
  }

  double vector2d::inner_product(vector2d value)
  {
    return value.x*x+value.y*y;
  }

  void vector2d::print_to_stream(std::ostream& stream) const
  {
    stream << "{" << x << ";" << y << "}";
    stream << std::flush;
  }

  std::ostream& operator<<(std::ostream& stream, const vector2d& vector)
  {
      vector.print_to_stream(stream);
      return stream;
   }
}
