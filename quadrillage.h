#include <eigen3/Eigen/Dense>
#include <iostream>

class Mesh
{
 private:
  Eigen::MatrixXd meshGrille;
 public:
  Mesh();
  ~Mesh();
  int getRows();
  int getCols();
  void print();
  void set(int row, int col, double value);
  void resize(int newRows, int newCols);
};
