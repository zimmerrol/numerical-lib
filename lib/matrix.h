namespace numerical{
  class Matrix
  {
    public:
      Matrix(int n, int m);
      Matrix(int n, int m, double defaultValue);
      ~Matrix();
      
      void swapRows(int targetm, int sourcem);
      void replaceRow(int targetm, double* source);
      void swapColumn(int targetn, int sourcen);
      void replaceColumn(int targetn, double* source);
      double determinant();
      void printToCout();
    private:
      double** data;
      int n;
      int m;
  };


}
