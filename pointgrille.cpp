#include "pointgrille.h"

TermVec::TermVec()
{
  cont = 0.0;
  elanX = 0.0;
  elanY = 0.0;
  elanZ = 0.0;
  nrg = 0.0;
}

TermVec::TermVec(double c, double x, double y, double z, double e)
{
  cont = c;
  elanX = x;
  elanY = y;
  elanZ = z;
  nrg = e;
}

TermVec::~TermVec()
{
}

PointGrille::PointGrille()
{
}

PointGrille::~PointGrille()
{
}

void PointGrille::initICs(){
  //
}
