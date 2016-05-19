#include <fstream>

namespace numerical{
  class matrix
  {
    public:
      matrix(int n, int m);
      matrix(int n, int m, double defaultValue);
      ~matrix();

      void swap_rows(int targetm, int sourcem);
      void replace_row(int targetm, double* source);
      void swap_column(int targetn, int sourcen);
      void replace_column(int targetn, double* source);
      double determinant();
      void print_to_stream(std::ostream& stream) const;
    private:
      double** data;
      int n;
      int m;
  };
}
