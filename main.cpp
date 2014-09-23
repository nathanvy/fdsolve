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

//grid stuff, to-do:  move this to ini file
// also make them global, fuck
const double machInf = 4.0;
const double nx = 70;
const double ny = 70;
const double length = 0.00001;
const double Gamma = 1.4;
const double aInf = 340.28;
const double pInf = 101325.0;
const double tempInf = 288.16;
const double temp0 = 288.16;
const double Pr = 0.71;
const double R = 287;
const double tempRatio = 1.0;
const double mu0 = 1.7894 * pow((double)10, (double)-5);
const double Courant = 0.6;

double timeStep(PointGrille* pg, double dx, double dy, int nx, int ny) {
  //Courant-Friedrichs-Levy
  double vPrimeMax = 0;
  double cflmin = 0;
  double u, v, a, mu, vprime, cfl;
  
  // mesh[i][j] is now mesh[i*ny + j]
  // e.g. p_mesh[i*ny+j].U.cont = 0;

  //see Anderson p 457
  for(int i=0;i<nx;i++){
    for(int j=0;j<ny;j++){
      u = pg[i*ny+j].U.elanX / pg[i*ny+j].U.cont;
      v = pg[i*ny+j].U.elanY / pg[i*ny+j].U.cont;
      a = sqrt(Gamma * R * pg[i*ny+j].T);
      mu = pg[i*ny+j].mu;
      vprime = (4/3)*mu*(Gamma*mu/Pr)/pg[i*ny+j].U.cont;
      
      if ( vprime > vPrimeMax ) {
        vPrimeMax = vprime;
      }
    }
  }

  for(int i=0;i<nx;i++){
    for(int j=0;j<ny;j++){
      cfl = 1/( (abs(u)/dx) + 
                (abs(v)/dy) + 
                a*sqrt( pow(dx,-2)+ pow(dy,-2) ) +
                2*vPrimeMax*( pow(dx,-2)+ pow(dy,-2) )
                ); //yuck
      
      if( (cfl < cflmin) || (cflmin == 0) ) {
        cflmin = cfl;
      }
    }
  }
    
  return Courant*cflmin;
}

bool checkConvergence(){
  return false;
}

double mac() {
  //calculate delU/delt at t,i,j from fwd differencing in space
  //predict U at time t+delta t
  //get predicted delU/delt at time t+delta t
  //run that through eqn 10.6 to get corrected, 2nd-order-accurate values
  //repeat until residuals drop below desired threshhold

  //after each predictor corrector, call decode(U);
  return 0;
}

void decode(PointGrille* pg) {
  //e.g. rho = U1
  // u = rho u / rho = U2/U1
  // cv = r/gamma-1
  // cp = gamma*cv
  // ...
  // gives rho, u, v, and e
  // lets us determine remaining unknowns T, p, mu, k
  double u, v, e;
  u = pg[i*ny+j].U.elanX / pg[i*ny+j].U.cont;
  v = pg[i*ny+j].U.elanY / pg[i*ny+j].U.cont;
  e = (pg[i*ny+j].U.nrg / pg[i*ny+j].U.cont) - 0.5*(pow(u,2)+pow(v,2));

  //set the actual flow variables to close the system
  pg[i*ny+j].T = e / (R/(Gamma - 1)); // T=e/Cv
  pg[i*ny+j].p = pg[i*ny+j].U.cont * R * pg[i*ny+j].T; //p = rho R T
  pg[i*ny+j].mu = mu0 * pow(pg[i*ny+j].T / temp0, 1.5) * (temp0+110) / (pg[i*ny+j].T + 110); //Sutherland's
  pg[i*ny+j].k = pg[i*ny+j].mu * (Gamma* R/(Gamma - 1) ) / Pr;
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
 
  int iterations;
  const int maxiterations = 10000;
  
  double Cv = R / (Gamma - 1);
  double Cp = Gamma * Cv;
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
    deltaT = timeStep(p_mesh, deltaX, deltaY, (int)nx, (int)ny);
    mac();
    
    if( checkConvergence() ) {
      break;
    }

  }

  checkContinuity();
  
  //write data to gnuplot somewhere here  

  return 0;
}
