#include <cmath>
#include "vector3d.h"

namespace numerical{
  vector3d::vector3d()
  {
    this->x = 0;
    this->y = 0;
    this->z = 0;
  }

  vector3d::vector3d(double x, double y, double z)
  {
    this->x = x;
    this->y = y;
    this->z = z;
  }

  vector3d::vector3d(const vector3d &value)
  {
    this->x = value.x;
    this->y = value.y;
    this->z = value.z;
  }

  vector3d::~vector3d()
  {

  }

  double vector3d::length()
  {
    return sqrt(product(*this));
  }

  void vector3d::add(vector3d value)
  {
    x += value.x;
    y += value.y;
    z += value.z;
  }

  void vector3d::substract(vector3d value)
  {
    x -= value.x;
    y -= value.y;
    z -= value.z;
  }

  void vector3d::scalar_multiply(double value)
  {
    x *= value;
    y *= value;
    z *= value;
  }

  vector3d vector3d::cross_product(vector3d value)
  {
    return vector3d(y*value.z-z*value.y, z*value.x-x*value.z, x*value.y-y*value.x);
  }

  double vector3d::inner_product(vector3d value)
  {
    return value.x*x+value.y*y+value.z*z;
  }

  void vector3d::print_to_stream(std::ostream& stream) const
  {
    stream << "{" << x << ";" << y << ";" << z << "}";
    stream << std::flush;
  }

  std::ostream& operator<<(std::ostream& stream, const vector3d& vector)
  {
      vector.print_to_stream(stream);
      return stream;
   }
}
