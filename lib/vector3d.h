#include <fstream>

namespace numerical{
  class vector3d
  {
    public:
      vector3d();
      vector3d(const vector3d &value);
      vector3d(double x, double y, double z);
      ~vector3d();

      double length();
      void add(vector3d value);
      void substract(vector3d value);
      void scalar_multiply(double value);
      vector3d cross_product(vector3d value);
      double inner_product(vector3d value);
      void print_to_stream(std::ostream& stream) const;
    private:
      double x;
      double y;
      double z;
  };
}
