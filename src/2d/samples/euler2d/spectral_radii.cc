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
// spctral_radii.cc -- an auxiliary function of
//
// Function called by time_step_2d.cc, Hx_flux.cc and Hz_flux.cc (model 
// specific)
//
// This function calculates the spectral radii of the Jacobian matrices of
// f(u) and g(u) for ideal MHD equations
// 
// Input (passed by reference):
//
// (1) u -- a doublearray1d type variable holding the conserved
//          quantities density, momentum (3 components), magnetic field
//          (3 components), and total energy over the entire soltion domain
//          at a cell of the solution domain
//
// (2) gamma -- a double type variable holding the ratio of specific heats
//
// Output (returned by reference):
//
// (1) rx -- a double type variable holding the spectral radius of the 
//           Jacobian matrix of f(u)
//
// (2) ry -- a double type variable holding the spectral radius of the 
//           Jacobian matrix of g(u)
//
////////////////////////////////////////////////////////////////////////////////

#include "centpack_2d_FD2.h"
#include "centpack_2d_SD2.h"

using namespace std;

void CENTPACK::spectral_radii(const doublearray1d& u, const doublearray1d& parameters, double& rx, double& ry)
{
	double rho, vx, vy, p, H, c;
	double gamma = parameters(0);
	
	rho = u(0);
	vx = u(1)/rho;
	vy = u(2)/rho;
	p = (gamma - 1.0)*(u(3) - 0.5*rho*(pow(vx,2) + pow(vy,2)));
	c = sqrt(gamma*p/rho);
	rx = fabs(vx) + c;
	ry = fabs(vy) + c;
}
