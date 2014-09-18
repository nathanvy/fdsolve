#include "quadrillage.h"

int Mesh::getRows()
{
  return meshGrille.rows();
}
int Mesh::getCols()
{
  return meshGrille.cols();
}

void Mesh::resize(int newRows, int newCols)
{
  //relatif à l'élément en haut, à gauche
  meshGrille.conservativeResize(newRows, newCols);
}

void Mesh::set(int row, int col, double value)
{
  meshGrille(row, col) = value;
}
void Mesh::print()
{
  std::cout << meshGrille << std::endl;
}

Mesh::Mesh()
{
}
Mesh::~Mesh()
{
}
