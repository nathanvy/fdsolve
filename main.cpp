// fdsolve - ©2014 Nathan Van Ymeren
//solutionneur des équations Navier-Stokes
//Navier-Stokes solver
//GPLv3

#include <fstream>
#include <iostream>
#include <exception>
//#include <eigen3/Eigen/Dense>
#include <cmath>

//#include "quadrillage.h"
#include "pointgrille.h"
#include "readwrite.h"

double timeStep() {
  //Courant-Friedrichs-Lewy
  return 0.1;
}

bool checkConvergence(){
  return false;
}

double mac() {
  return 0;
}

void checkContinuity() {
}

int main(int argc, char* argv[])
{
  /*
  if( argc != 5 )
    {
      std::cout << "Usage: " << argv[0] << " <continuity data> <x-momentum data> <y-momentum data> <energy data>" << std::endl;

      //tuer le processus toute suite
      return 1;
    }
  */
  
  //std::string cont = argv[1];
  //Mesh* p_contMesh = new Mesh();
  //ReaderWriter* p_rw = new ReaderWriter(cont, p_contMesh);
  
  //grid stuff, to-do:  move this to ini file
  const double machInf = 4.0;
  const double nx = 70;
  const double ny = 70;
  const double length = 0.00001;
  const double gamma = 1.4;
  const double aInf = 340.28;
  const double pInf = 101325.0;
  const double tempInf = 288.16;
  const double temp0 = 288.16;
  const double Pr = 0.71;
  const double R = 287;
  const double tempRatio = 1.0;
  const double mu0 = 1.7894 * pow((double)10, (double)-5);
  const double Courant = 0.6;
 
  int iterations;
  const int maxiterations = 10000;
  
  double Cv = R / (gamma - 1);
  double Cp = gamma * Cv;
  double rhoInf = pInf / ( R*tempInf);
  double eInf = Cv*tempInf;
  double muInf = mu0 * pow(tempInf/temp0, 1.5) * (temp0+110)/(tempInf+110);
  double ReL = rhoInf * ( machInf*aInf ) * length / muInf;
  double kInf = muInf * Cp / Pr;

  double deltaX = length / (nx-1);
  double deltaY = 5 * ( 5*length/ReL ) / (ny-1);
  double deltaT;

  // mesh[i][j] is now mesh[i*ny + j]
  // e.g. p_mesh[i*ny+j].U.cont = 0;
  PointGrille* p_mesh = new PointGrille[(int)(nx*ny)];
  
  // start iterating
  for(iterations = 0; iterations >= maxiterations; iterations++) {
    if(iterations == 0){
      p_mesh->initICs();
    }
    //calculate time step deltaT
    deltaT = timeStep();
    mac();
    
    if( checkConvergence() ) {
      break;
    }

  }

  checkContinuity();
  
  //write data to gnuplot somewhere here  

  return 0;
}
