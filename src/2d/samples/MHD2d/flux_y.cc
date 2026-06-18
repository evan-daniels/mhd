////////////////////////////////////////////////////////////////////////////////
//
// CentPack -- A generic numerical solver for hyperbolic conservation
//             laws and related time dependent problems
//
// Copyright (C) 2005 Jorge Balbas and Eitan Tadmor
//
// This program is distributed freely for research and instructional use only. 
// You may copy, modify,  and use this software for these non-commercial
// purposes, provided that the copyright notice and associated text is 
// reproduced on all copies
// 
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE.
//
// flux_y.cc -- Ideal MHD y-flux function
//
// Function called by Hy_flux, (model specific)
//
// Input (passed by reference):
//
// (1) u -- a doublearray1d type variable holding the conserved
//          quantities density, momentum (3 components), magnetic field
//          (3 components), and total energy at a point of the didcretized 
//          solution domain.
//
// (2) parameters -- a doublearray1d type variable holding additional parameters 
//					for the specific example
//
// Output (returned by reference):
//
// (1) g -- a doublearray1d type variable holding g(u), the flux in
//          the y direction.
//
////////////////////////////////////////////////////////////////////////////////

#include"centpack_2d_FD2.h"
#include"centpack_2d_SD2.h"

using namespace std;

void CENTPACK::flux_y(const doublearray1d& u, const doublearray1d& parameters, const doublearray1d& e_here, doublearray1d& g)
{
  double p, p_star;
  double gamma = parameters(0);

  double rho = u(0);
  double vx  = u(1)/rho;
  double vy  = u(2)/rho;
  double vz  = u(3)/rho;
  double Bx  = u(4);
  double By  = u(5);
  double Bz  = u(6);

  // src below
  // p = (gamma-1.0)*(u(7) - 0.5*((pow(u(1),2.0) + pow(u(2),2.0) + pow(u(3),2.0))/u(0)) - .5*(pow(u(4),2.0) + pow(u(5),2.0) + pow(u(6),2.0)));
  // p_star = p + 0.5*(pow(u(4),2.0) + pow(u(5),2.0) + pow(u(6),2.0));
  // src above

  p = (gamma-1.0)*(u(7) - 0.5*(u(1)*vx + u(2)*vy + u(3)*vz) - 0.5*(Bx*Bx + By*By + Bz*Bz));
  p_star = p + 0.5*(Bx*Bx + By*By + Bz*Bz);

  double e1 = e_here(0);  // Ex
  double e2 = e_here(1);  // Ey
  double e3 = e_here(2);  // Ez

  // ideal MHD momentum and density fluxes
  g(0) = u(2);                          // rho*vy
  g(1) = u(1)*vy - Bx*By;               // rho*vx*vy - Bx*By
  g(2) = u(2)*vy + p_star - By*By;      // rho*vy^2 + p* - By^2
  g(3) = u(3)*vy - By*Bz;               // rho*vz*vy - By*Bz

  // induction equation (dB/dt = -curl(E))
  g(4) =  e3;
  g(5) = 0.0;
  g(6) = -e1;

  // energy flux
  g(7) = (u(7) + p_star)*vy - By*(Bx*vx + By*vy + Bz*vz) + (e3*Bx - e1*Bz);

  // src below
  // g(0) = u(3);
  // g(1) = (u(1)*u(3)/u(0)) - u(4)*u(6);
  // g(2) = (u(2)*u(3)/u(0)) - u(5)*u(6);
  // g(3) = (pow(u(3),2))/u(0) + p_star - pow(u(6),2);
  // g(4) = (u(4)*u(3)/u(0)) - u(6)*u(1)/u(0);
  // g(5) = (u(5)*u(3)/u(0)) - u(6)*u(2)/u(0);
  // g(6) = 0.0;
  // g(7) = (u(7) + p_star)*(u(3)/u(0)) - u(6)*(u(4)*(u(1)/u(0)) + (u(2)/u(0))*u(5)+(u(3)/u(0))*u(6));
  // src above
}
