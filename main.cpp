// fdsolve - ©2014 Nathan Van Ymeren
//solutionneur des équations Navier-Stokes
//Navier-Stokes solver
//GPLv3

/*
  This program currently only solves the case of supersonic flow over a flate plate using the MacCormack Technique, a problem for which the author believes a time-marching technique is well-posed.  Future versions may or may not solve different geometries without needing to be recompiled.  The value to the author is the chance to get down and dirty with CFD, not to worry about the vagarities of parsing a shitload of input data
*/

#include <fstream>
#include <iostream>
#include <exception>
//#include <eigen3/Eigen/Dense>
#include <cmath>

//#include "quadrillage.h"
#include "pointgrille.h"
#include "readwrite.h"

//grid stuff, to-do:  move this to ini file
const double machInf = 4.0;
const int nx = 70;
const int ny = 70;
const double length = 0.00001; //length of the plate
const double Gamma = 1.4;
const double aInf = 340.28;
const double pInf = 101325.0;
const double tempInf = 288.16;
const double temp0 = 288.16;
const double Pr = 0.71; //Prandtl
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

void applyBCs(PointGrille* pg){
}

void decode(PointGrille* pg, int i, int j) {
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

double tauxx(PointGrille* pg, int i, int j){
  return 0;
}
double tauxy(PointGrille* pg, int i, int j){
  return 0;
}
double tauyy(PointGrille* pg, int i, int j){
  return 0;
}

double qx(PointGrille* pg, int i, int j){
  return 0;
}
double qy(PointGrille* pg, int i, int j){
  return 0;
}

bool checkContinuity() {
  return true;
}

void mac(PointGrille* pg, PointGrille* scratch, double dt, double dx, double dy, int nx, int ny) {
  for(int i=0;i<nx;i++){
    for(int j=0;i<ny;j++){
      //assume we know U at the current time step.
      decode(pg, i, j);
      //store current solution vector in scratch[] for use later
      scratch[i*ny+j].U.cont = pg[i*ny+j].U.cont;
      scratch[i*ny+j].U.elanX = pg[i*ny+j].U.elanX;
      scratch[i*ny+j].U.elanY = pg[i*ny+j].U.elanY;
      scratch[i*ny+j].U.nrg = pg[i*ny+j].U.nrg;
      
      // get E and F vectors since we need those for predicting
      pg[i*ny+j].E.cont = pg[i*ny+j].U.elanX; //rho u
      pg[i*ny+j].E.elanX = pow(pg[i*ny+j].U.elanX, 2) + pg[i*ny+j].p - tauxx(pg, i , j); //rho u^2 +p - tauxx
      pg[i*ny+j].E.elanY = ( pg[i*ny+j].U.elanX*pg[i*ny+j].U.elanY/pg[i*ny+j].U.cont ) - tauxy(pg, i, j); //rho u v - tauxy
      pg[i*ny+j].E.nrg = ( (pg[i*ny+j].U.nrg + pg[i*ny+j].p)*pg[i*ny+j].U.elanX/pg[i*ny+j].U.cont ) - ( (pg[i*ny+j].U.elanX/pg[i*ny+j].U.cont)*tauxx(pg, i, j) ) - ( (pg[i*ny+j].U.elanY/pg[i*ny+j].U.cont)*tauxy(pg, i, j) ) + qx(pg, i , j); //(Et+p)u - u tauxx - v tauxy + qx
      
      pg[i*ny+j].F.cont = pg[i*ny+j].U.elanY; //rho v
      pg[i*ny+j].F.elanX = ( pg[i*ny+j].U.elanX*pg[i*ny+j].U.elanY/pg[i*ny+j].U.cont ) - tauxy(pg, i, j); //rho u v - tauxy
      pg[i*ny+j].F.elanY = pow(pg[i*ny+j].U.elanY, 2) + pg[i*ny+j].p - tauxx(pg, i , j); //rho v^2 +p - tauxx
      pg[i*ny+j].F.nrg = ( (pg[i*ny+j].U.nrg + pg[i*ny+j].p)*(pg[i*ny+j].U.elanY/pg[i*ny+j].U.cont) ) - ( (pg[i*ny+j].U.elanX/pg[i*ny+j].U.cont)*tauxy(pg, i, j) ) - ( (pg[i*ny+j].U.elanY/pg[i*ny+j].U.cont)*tauyy(pg, i, j)) + qy(pg, i , j); //(Et+p)v - u tauxy - v tauyy + qy      
    }
  }
  
  //now, for interior points only
  //predict values of U at time t+delta:t
  for(int i=1;i<nx-1;i++){
    for(int j=1;j<nx-1;j++){
      pg[i*ny+j].U.cont = pg[i*ny+j].U.cont -
	(dt/dx)*( pg[(i+1)*ny+j].E.cont - pg[i*ny+j].E.cont ) -
	(dt/dy)*( pg[i*ny+(j+1)].F.cont - pg[i*ny+j].F.cont );
      pg[i*ny+j].U.elanX = pg[i*ny+j].U.elanX -
	(dt/dx)*( pg[(i+1)*ny+j].E.elanX - pg[i*ny+j].E.elanX ) -
	(dt/dy)*( pg[i*ny+(j+1)].F.elanX - pg[i*ny+j].F.elanX );
      pg[i*ny+j].U.elanY = pg[i*ny+j].U.elanY -
	(dt/dx)*( pg[(i+1)*ny+j].E.elanY - pg[i*ny+j].E.elanY ) -
	(dt/dy)*( pg[i*ny+(j+1)].F.elanY - pg[i*ny+j].F.elanY );
      pg[i*ny+j].U.nrg = pg[i*ny+j].U.nrg -
	(dt/dx)*( pg[(i+1)*ny+j].E.nrg - pg[i*ny+j].E.nrg ) -
	(dt/dy)*( pg[i*ny+(j+1)].F.nrg - pg[i*ny+j].F.nrg );
      
      decode(pg, i, j); //close the system with new predicted values
      
      //apply boundary conditions
      applyBCs(pg);

      //now calculate new E and F vectors with predicted values
      pg[i*ny+j].E.cont = pg[i*ny+j].U.elanX; //rho u
      pg[i*ny+j].E.elanX = pow(pg[i*ny+j].U.elanX, 2) + pg[i*ny+j].p - tauxx(pg, i , j); //rho u^2 +p - tauxx
      pg[i*ny+j].E.elanY = ( pg[i*ny+j].U.elanX*pg[i*ny+j].U.elanY/pg[i*ny+j].U.cont ) - tauxy(pg, i, j); //rho u v - tauxy
      pg[i*ny+j].E.nrg = ( (pg[i*ny+j].U.nrg + pg[i*ny+j].p)*pg[i*ny+j].U.elanX/pg[i*ny+j].U.cont ) - ( (pg[i*ny+j].U.elanX/pg[i*ny+j].U.cont)*tauxx(pg, i, j) ) - ( (pg[i*ny+j].U.elanY/pg[i*ny+j].U.cont)*tauxy(pg, i, j) ) + qx(pg, i , j); //(Et+p)u - u tauxx - v tauxy + qx
      
      pg[i*ny+j].F.cont = pg[i*ny+j].U.elanY; //rho v
      pg[i*ny+j].F.elanX = ( pg[i*ny+j].U.elanX*pg[i*ny+j].U.elanY/pg[i*ny+j].U.cont ) - tauxy(pg, i, j); //rho u v - tauxy
      pg[i*ny+j].F.elanY = pow(pg[i*ny+j].U.elanY, 2) + pg[i*ny+j].p - tauxx(pg, i , j); //rho v^2 +p - tauxx
      pg[i*ny+j].F.nrg = ( (pg[i*ny+j].U.nrg + pg[i*ny+j].p)*(pg[i*ny+j].U.elanY/pg[i*ny+j].U.cont) ) - ( (pg[i*ny+j].U.elanX/pg[i*ny+j].U.cont)*tauxy(pg, i, j) ) - ( (pg[i*ny+j].U.elanY/pg[i*ny+j].U.cont)*tauyy(pg, i, j)) + qy(pg, i , j); //(Et+p)v - u tauxy - v tauyy + qy

      applyBCs(pg);
    }
  }
  //entire flow should now be in predicted state at time t+deltaT

  //correct the flow field using opposite (rearward) diferencing
  // note that pg now holds the predicted values at time t+dt aka the "next" flow field and scratch holds the current flow field
  for(int i=0;i<nx;i++){
    for(int j=0;i<ny;j++){
      pg[i*ny+j].U.cont = 0.5*( scratch[i*ny+j].U.cont +
				pg[i*ny+j].U.cont -
				(dt/dx)*(pg[i*ny+j].E.cont - pg[(i-1)*ny+j].E.cont) -
				(dt/dy)*(pg[i*ny+j].F.cont - pg[i*ny+j-1].F.cont)
				);
      pg[i*ny+j].U.elanX = 0.5*( scratch[i*ny+j].U.elanX +
				pg[i*ny+j].U.elanX -
				(dt/dx)*(pg[i*ny+j].E.elanX - pg[(i-1)*ny+j].E.elanX) -
				(dt/dy)*(pg[i*ny+j].F.elanX - pg[i*ny+j-1].F.elanX)
				);
      pg[i*ny+j].U.elanY = 0.5*( scratch[i*ny+j].U.elanY +
				pg[i*ny+j].U.elanY -
				(dt/dx)*(pg[i*ny+j].E.elanY - pg[(i-1)*ny+j].E.elanY) -
				(dt/dy)*(pg[i*ny+j].F.elanY - pg[i*ny+j-1].F.elanY)
				);
      pg[i*ny+j].U.nrg = 0.5*( scratch[i*ny+j].U.nrg +
				pg[i*ny+j].U.nrg -
				(dt/dx)*(pg[i*ny+j].E.nrg - pg[(i-1)*ny+j].E.nrg) -
				(dt/dy)*(pg[i*ny+j].F.nrg - pg[i*ny+j-1].F.nrg)
				);

      decode(pg, i, j);
      
      //we _should_ be at a predicted-corrected t+dt flow field now

    }
  }

  //repeat until residuals drop below desired threshhold
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
 
  std::cout << "Preliminary setup... ";

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
  PointGrille* p_mesh = new PointGrille[nx*ny];
  PointGrille* p_scratch = new PointGrille[nx*ny]; //we've got shit tons of memory bro
  
  std::cout << " done." << std::endl;
  std::cout << "Beginning a series of max " << maxiterations << " iterations ..." << std:: endl;
  // start iterating
  for(iterations = 0; iterations >= maxiterations; iterations++) {
    if(iterations == 0){
      p_mesh->initICs();
    }
    //calculate time step deltaT
    deltaT = timeStep(p_mesh, deltaX, deltaY, nx, ny);
    mac(p_mesh, p_scratch, deltaT, deltaX, deltaY, nx, ny);
    
    if( checkConvergence() ) {
      std::cout << "Solution converged after " << iterations << " iterations!" << std::endl;
      break;
    }

  }

  if (checkContinuity() ) {
  }
  
  //write data to gnuplot somewhere here  

  return 0;
}
