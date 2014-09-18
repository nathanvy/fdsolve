#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>

//#include "quadrillage.h"

class ReaderWriter
{
 private:
  std::string filename;
 public:
  ReaderWriter();
  //ReaderWriter(std::string filename, Mesh* mesh);
  ~ReaderWriter();
  //void writeFile(std::string filename, Mesh* mesh);
};
