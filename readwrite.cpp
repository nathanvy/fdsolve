#include "readwrite.h"

ReaderWriter::ReaderWriter()
{
  std::cout << "default rw const" << std::endl;
}

ReaderWriter::~ReaderWriter()
{
  std::cout << "default rw dest" << std::endl;
}

/*ReaderWriter::ReaderWriter(std::string filename, Mesh* mesh)
{
  //get columns of matrix
  int cols = 0;
  int rows = 0;
  int i=0;
  int j=0;

  std::fstream infile(filename);
  std::string line;
  if(infile.good()){
    getline(infile, line);
    cols = std::count(line.begin(), line.end(), ',');
  }
  
  //rewind
  infile.clear();
  infile.seekg(0, std::ios::beg );
  line = "";

  //get rows of matrix
  while(getline(infile, line)){
    if( !line.empty() ){
      rows++;
    }
  }
  mesh->resize(rows,cols);

  //rewind again
  infile.clear();
  infile.seekg(0, std::ios::beg );
  line = "";

  while (infile) {
    if(!getline( infile, line)) {
      std::cout << "1" << std::endl;
      break;
    }    
   
    std::istringstream ss(line);
    
    while(ss){
      if(!getline( ss, line, ',' )) {
	std::cout << "2" << std::endl;
	break;
      }
      mesh->set(i,j, std::stod(line) );
      j++;
    }
    j=0;
    i++;
  }

  //mesh->print();
}

void ReaderWriter::writeFile(std::string filename, Mesh* mesh)
{
  //write to gnuplot-compatible format
}
*/
