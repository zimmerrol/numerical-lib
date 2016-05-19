#include <fstream>

namespace numerical{
  class vector2d
  {
    public:
      vector2d();
      vector2d(const vector2d &value);
      vector2d(double x, double y);
      ~vector2d();

      double length();
      void add(vector2d value);
      void substract(vector2d value);
      void scalar_multiply(double value);
      double inner_product(vector2d value);
      void print_to_stream(std::ostream& stream) const;
    private:
      double x;
      double y;
  };
}
