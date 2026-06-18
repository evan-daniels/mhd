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
// flux_x.cc -- Ideal MHD x-flux function
//
// Function called by Hx_flux, (model specific)
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
// (1) f -- a doublearray1d type variable holding f(u), the flux in
//          the x direction.
//
////////////////////////////////////////////////////////////////////////////////

#include"centpack_2d_FD2.h"
#include"centpack_2d_SD2.h"

using namespace std;

void CENTPACK::flux_x(const doublearray1d& u, const doublearray1d& parameters, const doublearray1d& e_here, doublearray1d& f)
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

  p = (gamma-1.0)*(u(7) - 0.5*(u(1)*vx + u(2)*vy + u(3)*vz) - 0.5*(Bx*Bx + By*By + Bz*Bz));
  p_star = p + 0.5*(Bx*Bx + By*By + Bz*Bz);

  // src below
  // p = (gamma-1.0)*(u(7) - 0.5*((pow(u(1),2.0) + pow(u(2),2.0) + pow(u(3),2.0))/u(0)) - .5*(pow(u(4),2.0) + pow(u(5),2.0) + pow(u(6),2.0)));
  // p_star = p + 0.5*(pow(u(4),2.0) + pow(u(5),2.0) + pow(u(6),2.0));
  // src above

  double e1 = e_here(0);  // Ex
  double e2 = e_here(1);  // Ey
  double e3 = e_here(2);  // Ez

  // ideal momentum and density fluxes
  f(0) = u(1);
  f(1) = u(1)*vx + p_star - Bx*Bx;
  f(2) = u(1)*vy - Bx*By;
  f(3) = u(1)*vz - Bx*Bz;
  
  // induction equation: dB/dt = -curl(E)
  f(4) = 0.0;
  f(5) = -e3;
  f(6) = e2;

  // energy flux
  f(7) = (u(7) + p_star)*vx - Bx*(Bx*vx + By*vy + Bz*vz) + (e2*Bz - e3*By);

  // SOURCE CODE BELOW
  // f(0) = u(1);
  // f(1) = (pow(u(1),2))/u(0) + p_star - pow(u(4),2);
  // f(2) = (u(1)*u(2)/u(0)) - u(4)*u(5);
  // f(3) = (u(1)*u(3)/u(0)) - u(4)*u(6);
  // f(4) = 0.0;
  // f(5) = (u(5)*u(1)/u(0)) - u(4)*u(2)/u(0);
  // f(6) = (u(6)*u(1)/u(0)) - u(4)*u(3)/u(0);
  // f(7) = (u(7) + p_star)*(u(1)/u(0)) - u(4)*(u(4)*(u(1)/u(0)) + (u(2)/u(0))*u(5)+(u(3)/u(0))*u(6));
  //SOURCE CODE ABOVE
}
